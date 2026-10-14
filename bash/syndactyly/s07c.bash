#!/bin/bash
# Human Whole Genome Sequencing (WGS) — Preprocessing, Mapping, BQSR, Variant Calling & Functional Annotation
# Workflow:
#   Phase 1 : Reference Genome Preparation & Indexing (BWA index, Picard CreateSequenceDictionary, Samtools faidx)
#   Phase 2 : Read Quality Control & Preprocessing (FastQC raw, Trimmomatic PE, FastQC trimmed)
#   Phase 3 : Alignment, Coordinate Sorting & Deduplication (BWA-MEM, Samtools sort, Picard MarkDuplicates, Samtools index)
#   Phase 4 : Base Quality Score Recalibration (GATK BaseRecalibrator, GATK ApplyBQSR)
#   Phase 5 : Germline Variant Calling & Hard Filtering (GATK HaplotypeCaller, SelectVariants, VariantFiltration)
#   Phase 6 : Variant Functional Annotation (ANNOVAR database setup, convert2annovar, table_annovar)
# Building No. 3
set -euo pipefail

################
# GLOBAL SETUP #
################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="${SCRIPT_DIR}"
SCRIPTS_DIR="${PROJECT_DIR}/scripts"

# Base directory paths (with environment override & fallback handling)
REF_PATH="${REF_PATH:-/storage/student9/references}"
TOOLS_PATH="${TOOLS_PATH:-/storage/student9/tools}"
SAMPLE_PATH="${SAMPLE_PATH:-/storage/student9/projects/s07c}"
[ ! -d "${SAMPLE_PATH}" ] && [ -d "/storage/student9/s07c" ] && SAMPLE_PATH="/storage/student9/s07c"

# Tunable execution parameters
threads="${threads:-16}"
sample_id="${1:-S07C}" # Sample ID passed as first argument, default S07C

# Reference genome, sequence dictionary, adapters, and known polymorphic sites (dbSNP)
humanref="${REF_PATH}/hg38.fa"
humanref_dict="${REF_PATH}/hg38.dict"
truseq3="${REF_PATH}/TruSeq3-PE.fa"
dbSNP="${REF_PATH}/common_all_20180418.vcf"

# ANNOVAR setup (custom installation - called directly via perl)
annovar_path="${TOOLS_PATH}/annovar"
annovar_downdb="${annovar_path}/annotate_variation.pl"
annovar_convert="${annovar_path}/convert2annovar.pl"
annovar_table="${annovar_path}/table_annovar.pl"
humandb="${REF_PATH}/humandb"
annovar_buildver="hg38"
annovar_protocol="refGene,cytoBand,exac03,avsnp150,dbnsfp30a"
annovar_operation="g,r,f,f,f"

# Per-tool workflow output directories
RAW_WGS="${SAMPLE_PATH}/raw"
PREPROCESSING_WGS="${SAMPLE_PATH}/preprocessing_WGS"
FASTQC_RAW="${PREPROCESSING_WGS}/fastqc_raw"
FASTQC_PP="${PREPROCESSING_WGS}/fastqc_pp"
MAPPING_WGS="${SAMPLE_PATH}/mapping_WGS"
MAPPING_SAMPLE="${MAPPING_WGS}/${sample_id}"
CALLING_WGS="${SAMPLE_PATH}/variant_calling_WGS"
CALLING_SAMPLE="${CALLING_WGS}/${sample_id}"
ANNOTATION_WGS="${SAMPLE_PATH}/annotation_WGS"
ANNOTATION_SAMPLE="${ANNOTATION_WGS}/${sample_id}"
LOGS_DIR="${SAMPLE_PATH}/logs"

# Create all workflow directories up front
mkdir -p "${LOGS_DIR}"
mkdir -p "${PREPROCESSING_WGS}"
mkdir -p "${FASTQC_RAW}"
mkdir -p "${FASTQC_PP}"
mkdir -p "${MAPPING_SAMPLE}"
mkdir -p "${CALLING_SAMPLE}"
mkdir -p "${ANNOTATION_SAMPLE}"
mkdir -p "${humandb}"

# Dual-log setup: main run log and dedicated failure/skip tracking log
LOG="${LOGS_DIR}/WGS_${sample_id}.log"
FAIL_LOG="${LOGS_DIR}/WGS_${sample_id}.failed_skipped.log"
exec > >(tee -a "${LOG}") 2>&1

echo "======================================================" >> "${FAIL_LOG}"
echo " Failure & Skip Log — Started: $(date)" >> "${FAIL_LOG}"
echo "======================================================" >> "${FAIL_LOG}"

# Helper function to record failures/skips to FAIL_LOG with timestamp
log_failure() {
    local phase="$1"
    local sample="$2"
    local status="$3" # e.g. "SKIPPED_EXISTS", "INPUT_MISSING", "EXECUTION_FAILED", "OUTPUT_MISSING"
    local reason="$4"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [${phase}] [${sample}] [${status}] ${reason}" >> "${FAIL_LOG}"
}

echo "============================================================"
echo " Human WGS Variant Calling & Annotation: ${sample_id} - Started: $(date)"
echo " Reference genome   : ${humanref}"
echo " Threads allocated  : ${threads}"
echo " Full log           : ${LOG}"
echo " Failure & skip log : ${FAIL_LOG}"
echo "============================================================"

    #############################################
    # PHASE 1: REFERENCE GENOME PREPARATION     #
    #############################################

    echo -e "\e[33m ====================================================== \e[0m"
    echo -e "\e[33m PHASE 1: REFERENCE GENOME PREPARATION & INDEXING       \e[0m"
    echo -e "\e[33m ====================================================== \e[0m"

    # -- 1a. BWA index ---------------------------------------------------------
    echo -e "\e[31m ===================== \e[0m"
    echo -e "\e[31m BWA INDEX: REFERENCE  \e[0m"
    echo -e "\e[31m ===================== \e[0m"

    expected_bwt="${humanref}.bwt"
    if [ -s "${expected_bwt}" ]; then
        echo -e "\e[32m BWA index already present — skipping \e[0m"
        log_failure "Phase1_RefPrep" "hg38" "SKIPPED_EXISTS" "BWA index ${expected_bwt} already exists"
    else
        if [ ! -s "${humanref}" ]; then
            echo -e "\e[31m ERROR: Reference FASTA missing or empty: ${humanref} \e[0m"
            log_failure "Phase1_RefPrep" "hg38" "INPUT_MISSING" "File ${humanref} missing or 0 bytes"
            exit 1
        fi

        echo -e "\e[32m Building BWA index for ${humanref}... \e[0m"
        if ! conda run -n mapping bwa index "${humanref}"; then
            echo -e "\e[31m ERROR: BWA index execution failed \e[0m"
            log_failure "Phase1_RefPrep" "hg38" "EXECUTION_FAILED" "bwa index exited with non-zero status"
            exit 1
        fi

        if [ ! -s "${expected_bwt}" ]; then
            echo -e "\e[31m ERROR: Expected BWA index missing after run: ${expected_bwt} \e[0m"
            log_failure "Phase1_RefPrep" "hg38" "OUTPUT_MISSING" "File ${expected_bwt} not created"
            exit 1
        fi
        echo -e "\e[32m BWA index completed -> ${expected_bwt} \e[0m"
    fi

    # -- 1b. Picard sequence dictionary ---------------------------------------
    echo -e "\e[31m ================================= \e[0m"
    echo -e "\e[31m PICARD: CREATE SEQUENCE DICTIONARY\e[0m"
    echo -e "\e[31m ================================= \e[0m"

    if [ -s "${humanref_dict}" ]; then
        echo -e "\e[32m Sequence dictionary already present — skipping \e[0m"
        log_failure "Phase1_RefPrep" "hg38" "SKIPPED_EXISTS" "Sequence dictionary ${humanref_dict} already exists"
    else
        if [ ! -s "${humanref}" ]; then
            echo -e "\e[31m ERROR: Reference FASTA missing or empty: ${humanref} \e[0m"
            log_failure "Phase1_RefPrep" "hg38" "INPUT_MISSING" "File ${humanref} missing or 0 bytes"
            exit 1
        fi

        echo -e "\e[32m Creating sequence dictionary for ${humanref}... \e[0m"
        if ! conda run -n mapping picard CreateSequenceDictionary \
            REFERENCE="${humanref}" \
            OUTPUT="${humanref_dict}"; then
            echo -e "\e[31m ERROR: Picard CreateSequenceDictionary execution failed \e[0m"
            log_failure "Phase1_RefPrep" "hg38" "EXECUTION_FAILED" "Picard CreateSequenceDictionary failed"
            exit 1
        fi

        if [ ! -s "${humanref_dict}" ]; then
            echo -e "\e[31m ERROR: Sequence dictionary missing after run: ${humanref_dict} \e[0m"
            log_failure "Phase1_RefPrep" "hg38" "OUTPUT_MISSING" "File ${humanref_dict} not created"
            exit 1
        fi
        echo -e "\e[32m Sequence dictionary completed -> ${humanref_dict} \e[0m"
    fi

    # -- 1c. Samtools faidx ----------------------------------------------------
    echo -e "\e[31m ======================== \e[0m"
    echo -e "\e[31m SAMTOOLS: FASTA INDEXING \e[0m"
    echo -e "\e[31m ======================== \e[0m"

    expected_fai="${humanref}.fai"
    if [ -s "${expected_fai}" ]; then
        echo -e "\e[32m FASTA index already present — skipping \e[0m"
        log_failure "Phase1_RefPrep" "hg38" "SKIPPED_EXISTS" "FASTA index ${expected_fai} already exists"
    else
        if [ ! -s "${humanref}" ]; then
            echo -e "\e[31m ERROR: Reference FASTA missing or empty: ${humanref} \e[0m"
            log_failure "Phase1_RefPrep" "hg38" "INPUT_MISSING" "File ${humanref} missing or 0 bytes"
            exit 1
        fi

        echo -e "\e[32m Indexing ${humanref} with samtools faidx... \e[0m"
        if ! conda run -n mapping samtools faidx "${humanref}"; then
            echo -e "\e[31m ERROR: samtools faidx execution failed \e[0m"
            log_failure "Phase1_RefPrep" "hg38" "EXECUTION_FAILED" "samtools faidx failed"
            exit 1
        fi

        if [ ! -s "${expected_fai}" ]; then
            echo -e "\e[31m ERROR: FASTA index missing after run: ${expected_fai} \e[0m"
            log_failure "Phase1_RefPrep" "hg38" "OUTPUT_MISSING" "File ${expected_fai} not created"
            exit 1
        fi
        echo -e "\e[32m FASTA index completed -> ${expected_fai} \e[0m"
    fi

    #############################################
    # PHASE 2: QUALITY CONTROL & PREPROCESSING  #
    #############################################

    echo -e "\e[33m ====================================================== \e[0m"
    echo -e "\e[33m PHASE 2: QUALITY CONTROL & READ PREPROCESSING          \e[0m"
    echo -e "\e[33m ====================================================== \e[0m"

    # Identify input FASTQ files with flexible naming fallback
    read1="${RAW_WGS}/${sample_id}_1.fastq.gz"
    read2="${RAW_WGS}/${sample_id}_2.fastq.gz"
    [ ! -s "${read1}" ] && [ -s "${RAW_WGS}/${sample_id}_R1.fastq.gz" ] && read1="${RAW_WGS}/${sample_id}_R1.fastq.gz"
    [ ! -s "${read2}" ] && [ -s "${RAW_WGS}/${sample_id}_R2.fastq.gz" ] && read2="${RAW_WGS}/${sample_id}_R2.fastq.gz"
    [ ! -s "${read1}" ] && [ -s "${RAW_WGS}/${sample_id}_1.fq.gz" ] && read1="${RAW_WGS}/${sample_id}_1.fq.gz"
    [ ! -s "${read2}" ] && [ -s "${RAW_WGS}/${sample_id}_2.fq.gz" ] && read2="${RAW_WGS}/${sample_id}_2.fq.gz"

    # -- 2a. FastQC Raw Reads --------------------------------------------------
    echo -e "\e[31m ========================== \e[0m"
    echo -e "\e[31m FASTQC: ${sample_id} (RAW) \e[0m"
    echo -e "\e[31m ========================== \e[0m"

    r1_base=$(basename "${read1}" | sed -E 's/\.(fastq|fq)(\.gz)?$//')
    r2_base=$(basename "${read2}" | sed -E 's/\.(fastq|fq)(\.gz)?$//')
    fastqc_raw_html1="${FASTQC_RAW}/${r1_base}_fastqc.html"
    fastqc_raw_html2="${FASTQC_RAW}/${r2_base}_fastqc.html"

    if [ -s "${fastqc_raw_html1}" ] && [ -s "${fastqc_raw_html2}" ]; then
        echo -e "\e[32m Raw FastQC reports already completed — skipping \e[0m"
        log_failure "Phase2_FastQC_Raw" "${sample_id}" "SKIPPED_EXISTS" "FastQC raw reports exist"
    else
        if [ ! -s "${read1}" ] || [ ! -s "${read2}" ]; then
            echo -e "\e[31m ERROR: Raw FASTQ files missing or empty for ${sample_id}: ${read1} / ${read2} \e[0m"
            log_failure "Phase2_FastQC_Raw" "${sample_id}" "INPUT_MISSING" "Raw reads missing or empty"
            exit 1
        fi

        echo -e "\e[32m Running FastQC on raw reads for ${sample_id}... \e[0m"
        if ! conda run -n preprocessing fastqc \
            --threads "${threads}" \
            --outdir "${FASTQC_RAW}" \
            "${read1}" "${read2}"; then
            echo -e "\e[31m ERROR: FastQC raw execution failed \e[0m"
            log_failure "Phase2_FastQC_Raw" "${sample_id}" "EXECUTION_FAILED" "fastqc raw non-zero exit code"
            exit 1
        fi

        if [ ! -s "${fastqc_raw_html1}" ] || [ ! -s "${fastqc_raw_html2}" ]; then
            echo -e "\e[31m ERROR: FastQC raw reports missing after run \e[0m"
            log_failure "Phase2_FastQC_Raw" "${sample_id}" "OUTPUT_MISSING" "Raw FastQC html reports not created"
            exit 1
        fi
        echo -e "\e[32m Raw FastQC completed -> ${FASTQC_RAW} \e[0m"
    fi

    # -- 2b. Trimmomatic Adapter & Quality Trimming -----------------------------
    echo -e "\e[31m ========================= \e[0m"
    echo -e "\e[31m TRIMMOMATIC: ${sample_id} \e[0m"
    echo -e "\e[31m ========================= \e[0m"

    trim_output_1="${PREPROCESSING_WGS}/${sample_id}.R1"
    trim_output_2="${PREPROCESSING_WGS}/${sample_id}.R2"
    read1t="${trim_output_1}.paired.fastq.gz"
    read2t="${trim_output_2}.paired.fastq.gz"
    read1u="${trim_output_1}.unpaired.fastq.gz"
    read2u="${trim_output_2}.unpaired.fastq.gz"

    if [ -s "${read1t}" ] && [ -s "${read2t}" ]; then
        echo -e "\e[32m Trimmed reads already present — skipping Trimmomatic \e[0m"
        log_failure "Phase2_Trimmomatic" "${sample_id}" "SKIPPED_EXISTS" "Trimmed reads exist"
    else
        if [ ! -s "${read1}" ] || [ ! -s "${read2}" ]; then
            echo -e "\e[31m ERROR: Raw reads missing for trimming: ${read1} / ${read2} \e[0m"
            log_failure "Phase2_Trimmomatic" "${sample_id}" "INPUT_MISSING" "Input raw reads missing"
            exit 1
        fi

        # Check adapter file presence
        if [ ! -s "${truseq3}" ]; then
            echo -e "\e[31m WARNING: Adapter file ${truseq3} not found. Attempting conda shared location... \e[0m"
            if [ -s "/storage/student9/miniconda3/envs/preprocessing/share/trimmomatic/adapters/TruSeq3-PE.fa" ]; then
                truseq3="/storage/student9/miniconda3/envs/preprocessing/share/trimmomatic/adapters/TruSeq3-PE.fa"
            elif [ -s "/storage/student9/tools/trimmomatic/adapters/TruSeq3-PE.fa" ]; then
                truseq3="/storage/student9/tools/trimmomatic/adapters/TruSeq3-PE.fa"
            fi
        fi

        echo -e "\e[32m Running Trimmomatic PE for ${sample_id}... \e[0m"
        if ! conda run -n preprocessing trimmomatic PE \
            -threads "${threads}" \
            -phred33 \
            "${read1}" "${read2}" \
            "${read1t}" "${read1u}" \
            "${read2t}" "${read2u}" \
            ILLUMINACLIP:"${truseq3}":2:30:10:8:true \
            HEADCROP:3 TRAILING:10 MINLEN:50; then
            echo -e "\e[31m ERROR: Trimmomatic execution failed \e[0m"
            log_failure "Phase2_Trimmomatic" "${sample_id}" "EXECUTION_FAILED" "trimmomatic non-zero exit code"
            exit 1
        fi

        if [ ! -s "${read1t}" ] || [ ! -s "${read2t}" ]; then
            echo -e "\e[31m ERROR: Trimmed paired reads missing after run: ${read1t} \e[0m"
            log_failure "Phase2_Trimmomatic" "${sample_id}" "OUTPUT_MISSING" "Trimmed paired output missing"
            exit 1
        fi
        echo -e "\e[32m Trimmomatic completed -> ${read1t} , ${read2t} \e[0m"
    fi

    # -- 2c. FastQC Trimmed Reads ----------------------------------------------
    echo -e "\e[31m ========================= \e[0m"
    echo -e "\e[31m FASTQC: ${sample_id} (PP) \e[0m"
    echo -e "\e[31m ========================= \e[0m"

    r1t_base=$(basename "${read1t}" | sed -E 's/\.(fastq|fq)(\.gz)?$//')
    r2t_base=$(basename "${read2t}" | sed -E 's/\.(fastq|fq)(\.gz)?$//')
    fastqc_pp_html1="${FASTQC_PP}/${r1t_base}_fastqc.html"
    fastqc_pp_html2="${FASTQC_PP}/${r2t_base}_fastqc.html"

    if [ -s "${fastqc_pp_html1}" ] && [ -s "${fastqc_pp_html2}" ]; then
        echo -e "\e[32m Trimmed FastQC reports already completed — skipping \e[0m"
        log_failure "Phase2_FastQC_PP" "${sample_id}" "SKIPPED_EXISTS" "FastQC trimmed reports exist"
    else
        if [ ! -s "${read1t}" ] || [ ! -s "${read2t}" ]; then
            echo -e "\e[31m ERROR: Trimmed reads missing for post-QC: ${read1t} / ${read2t} \e[0m"
            log_failure "Phase2_FastQC_PP" "${sample_id}" "INPUT_MISSING" "Trimmed reads missing or empty"
            exit 1
        fi

        echo -e "\e[32m Running FastQC on trimmed reads for ${sample_id}... \e[0m"
        if ! conda run -n preprocessing fastqc \
            --threads "${threads}" \
            --outdir "${FASTQC_PP}" \
            "${read1t}" "${read2t}"; then
            echo -e "\e[31m ERROR: FastQC trimmed execution failed \e[0m"
            log_failure "Phase2_FastQC_PP" "${sample_id}" "EXECUTION_FAILED" "fastqc trimmed non-zero exit code"
            exit 1
        fi

        if [ ! -s "${fastqc_pp_html1}" ] || [ ! -s "${fastqc_pp_html2}" ]; then
            echo -e "\e[31m ERROR: FastQC trimmed reports missing after run \e[0m"
            log_failure "Phase2_FastQC_PP" "${sample_id}" "OUTPUT_MISSING" "Trimmed FastQC html reports not created"
            exit 1
        fi
        echo -e "\e[32m Trimmed FastQC completed -> ${FASTQC_PP} \e[0m"
    fi

    #######################################################
    # PHASE 3: READ MAPPING, SORTING & DEDUPLICATION      #
    #######################################################

    echo -e "\e[33m ====================================================== \e[0m"
    echo -e "\e[33m PHASE 3: READ ALIGNMENT, SORTING & DEDUPLICATION       \e[0m"
    echo -e "\e[33m ====================================================== \e[0m"

    sorted_bam="${MAPPING_SAMPLE}/${sample_id}.sorted_coord.bam"
    dedup_bam="${MAPPING_SAMPLE}/${sample_id}.dedup.bam"
    dedup_metrics="${MAPPING_SAMPLE}/${sample_id}.dedup_metrics.txt"

    # -- 3a. BWA-MEM Alignment & Coordinate Sort -------------------------------
    echo -e "\e[31m ===================================== \e[0m"
    echo -e "\e[31m BWA-MEM & SAMTOOLS SORT: ${sample_id} \e[0m"
    echo -e "\e[31m ===================================== \e[0m"

    if [ -s "${dedup_bam}" ]; then
        echo -e "\e[32m Deduplicated BAM already exists — skipping alignment & duplicate marking \e[0m"
        log_failure "Phase3_Mapping" "${sample_id}" "SKIPPED_EXISTS" "Final BAM ${dedup_bam} already exists"
    elif [ -s "${sorted_bam}" ]; then
        echo -e "\e[32m Sorted BAM already exists — skipping BWA-MEM alignment \e[0m"
        log_failure "Phase3_BWA" "${sample_id}" "SKIPPED_EXISTS" "Sorted BAM ${sorted_bam} already exists"
    else
        if [ ! -s "${read1t}" ] || [ ! -s "${read2t}" ]; then
            echo -e "\e[31m ERROR: Trimmed reads missing for alignment: ${read1t} \e[0m"
            log_failure "Phase3_BWA" "${sample_id}" "INPUT_MISSING" "Trimmed reads missing or empty"
            exit 1
        fi
        if [ ! -s "${humanref}.bwt" ]; then
            echo -e "\e[31m ERROR: BWA index missing for ${humanref} \e[0m"
            log_failure "Phase3_BWA" "${sample_id}" "INPUT_MISSING" "BWA index missing"
            exit 1
        fi

        echo -e "\e[32m Aligning reads with BWA-MEM and streaming to samtools sort... \e[0m"
        # Streaming BWA-MEM directly to samtools sort eliminates massive uncompressed SAM file disk footprint
        if ! conda run -n mapping bwa mem \
            -t "${threads}" \
            -R "@RG\tID:${sample_id}\tSM:${sample_id}\tPL:illumina\tLB:lib1\tPU:unit1" \
            "${humanref}" \
            "${read1t}" "${read2t}" | \
            conda run -n mapping samtools sort \
            -@ "${threads}" \
            -o "${sorted_bam}" -; then
            echo -e "\e[31m ERROR: Alignment / sort pipeline failed \e[0m"
            log_failure "Phase3_BWA" "${sample_id}" "EXECUTION_FAILED" "bwa mem | samtools sort failed"
            exit 1
        fi

        if [ ! -s "${sorted_bam}" ]; then
            echo -e "\e[31m ERROR: Coordinate-sorted BAM missing or empty: ${sorted_bam} \e[0m"
            log_failure "Phase3_BWA" "${sample_id}" "OUTPUT_MISSING" "Sorted BAM missing or 0 bytes"
            exit 1
        fi
        echo -e "\e[32m Alignment and coordinate sort completed -> ${sorted_bam} \e[0m"
    fi

    # -- 3b. Picard MarkDuplicates ---------------------------------------------
    echo -e "\e[31m =========================== \e[0m"
    echo -e "\e[31m PICARD MARKDUPLICATES: ${sample_id} \e[0m"
    echo -e "\e[31m =========================== \e[0m"

    if [ -s "${dedup_bam}" ]; then
        echo -e "\e[32m Deduplicated BAM already exists — skipping MarkDuplicates \e[0m"
        log_failure "Phase3_MarkDuplicates" "${sample_id}" "SKIPPED_EXISTS" "Deduplicated BAM exists"
    else
        if [ ! -s "${sorted_bam}" ]; then
            echo -e "\e[31m ERROR: Sorted BAM missing for deduplication: ${sorted_bam} \e[0m"
            log_failure "Phase3_MarkDuplicates" "${sample_id}" "INPUT_MISSING" "Sorted BAM missing"
            exit 1
        fi

        echo -e "\e[32m Marking duplicate reads with Picard MarkDuplicates... \e[0m"
        if ! conda run -n mapping picard MarkDuplicates \
            INPUT="${sorted_bam}" \
            OUTPUT="${dedup_bam}" \
            METRICS_FILE="${dedup_metrics}" \
            REMOVE_DUPLICATES=true \
            CREATE_INDEX=true \
            VALIDATION_STRINGENCY=SILENT; then
            echo -e "\e[31m ERROR: Picard MarkDuplicates execution failed \e[0m"
            log_failure "Phase3_MarkDuplicates" "${sample_id}" "EXECUTION_FAILED" "Picard MarkDuplicates failed"
            exit 1
        fi

        if [ ! -s "${dedup_bam}" ]; then
            echo -e "\e[31m ERROR: Deduplicated BAM missing after run: ${dedup_bam} \e[0m"
            log_failure "Phase3_MarkDuplicates" "${sample_id}" "OUTPUT_MISSING" "dedup.bam missing or 0 bytes"
            exit 1
        fi
        echo -e "\e[32m Duplicate marking completed -> ${dedup_bam} \e[0m"
    fi

    # -- 3c. Samtools Index BAM ------------------------------------------------
    echo -e "\e[31m ======================== \e[0m"
    echo -e "\e[31m SAMTOOLS INDEX: ${sample_id} \e[0m"
    echo -e "\e[31m ======================== \e[0m"

    if [ -s "${dedup_bam}.bai" ]; then
        echo -e "\e[32m BAM index already exists — skipping \e[0m"
        log_failure "Phase3_Samtools_Index" "${sample_id}" "SKIPPED_EXISTS" "dedup_bam.bai already exists"
    else
        echo -e "\e[32m Indexing deduplicated BAM with samtools... \e[0m"
        if ! conda run -n mapping samtools index \
            -@ "${threads}" \
            "${dedup_bam}"; then
            echo -e "\e[31m ERROR: samtools index execution failed \e[0m"
            log_failure "Phase3_Samtools_Index" "${sample_id}" "EXECUTION_FAILED" "samtools index failed"
            exit 1
        fi
        echo -e "\e[32m BAM indexing completed -> ${dedup_bam}.bai \e[0m"
    fi

    #################################################
    # PHASE 4: BASE QUALITY SCORE RECALIBRATION     #
    #################################################

    echo -e "\e[33m ====================================================== \e[0m"
    echo -e "\e[33m PHASE 4: BASE QUALITY SCORE RECALIBRATION (BQSR)       \e[0m"
    echo -e "\e[33m ====================================================== \e[0m"

    recal_table="${MAPPING_SAMPLE}/${sample_id}.recal_data.grp"
    recal_bam="${MAPPING_SAMPLE}/${sample_id}.recal.bam"

    # -- 4a. GATK BaseRecalibrator ---------------------------------------------
    echo -e "\e[31m =============================== \e[0m"
    echo -e "\e[31m GATK BASERECALIBRATOR: ${sample_id} \e[0m"
    echo -e "\e[31m =============================== \e[0m"

    if [ -s "${recal_table}" ]; then
        echo -e "\e[32m BQSR recalibration table already present — skipping \e[0m"
        log_failure "Phase4_BaseRecalibrator" "${sample_id}" "SKIPPED_EXISTS" "recal_data.grp exists"
    else
        if [ ! -s "${dedup_bam}" ]; then
            echo -e "\e[31m ERROR: Deduplicated BAM missing for BQSR: ${dedup_bam} \e[0m"
            log_failure "Phase4_BaseRecalibrator" "${sample_id}" "INPUT_MISSING" "dedup.bam missing"
            exit 1
        fi
        if [ ! -s "${dbSNP}" ]; then
            echo -e "\e[31m ERROR: Known sites dbSNP VCF missing: ${dbSNP} \e[0m"
            log_failure "Phase4_BaseRecalibrator" "${sample_id}" "INPUT_MISSING" "dbSNP VCF missing"
            exit 1
        fi

        echo -e "\e[32m Computing BQSR recalibration table with GATK BaseRecalibrator... \e[0m"
        if ! conda run -n Hcalling gatk BaseRecalibrator \
            -R "${humanref}" \
            -I "${dedup_bam}" \
            --known-sites "${dbSNP}" \
            -O "${recal_table}"; then
            echo -e "\e[31m ERROR: GATK BaseRecalibrator execution failed \e[0m"
            log_failure "Phase4_BaseRecalibrator" "${sample_id}" "EXECUTION_FAILED" "BaseRecalibrator non-zero exit code"
            exit 1
        fi

        if [ ! -s "${recal_table}" ]; then
            echo -e "\e[31m ERROR: Recalibration table missing after run: ${recal_table} \e[0m"
            log_failure "Phase4_BaseRecalibrator" "${sample_id}" "OUTPUT_MISSING" "recal_table missing or 0 bytes"
            exit 1
        fi
        echo -e "\e[32m BQSR table completed -> ${recal_table} \e[0m"
    fi

    # -- 4b. GATK ApplyBQSR ----------------------------------------------------
    echo -e "\e[31m ======================== \e[0m"
    echo -e "\e[31m GATK APPLYBQSR: ${sample_id} \e[0m"
    echo -e "\e[31m ======================== \e[0m"

    if [ -s "${recal_bam}" ]; then
        echo -e "\e[32m Recalibrated BAM already exists — skipping ApplyBQSR \e[0m"
        log_failure "Phase4_ApplyBQSR" "${sample_id}" "SKIPPED_EXISTS" "recal.bam already exists"
    else
        if [ ! -s "${dedup_bam}" ] || [ ! -s "${recal_table}" ]; then
            echo -e "\e[31m ERROR: Inputs missing for ApplyBQSR: ${dedup_bam} or ${recal_table} \e[0m"
            log_failure "Phase4_ApplyBQSR" "${sample_id}" "INPUT_MISSING" "dedup.bam or recal_table missing"
            exit 1
        fi

        echo -e "\e[32m Applying BQSR recalibration to generate recalibrated BAM... \e[0m"
        if ! conda run -n Hcalling gatk ApplyBQSR \
            -R "${humanref}" \
            -I "${dedup_bam}" \
            --bqsr-recal-file "${recal_table}" \
            --create-output-bam-index true \
            -O "${recal_bam}"; then
            echo -e "\e[31m ERROR: GATK ApplyBQSR execution failed \e[0m"
            log_failure "Phase4_ApplyBQSR" "${sample_id}" "EXECUTION_FAILED" "ApplyBQSR non-zero exit code"
            exit 1
        fi

        if [ ! -s "${recal_bam}" ]; then
            echo -e "\e[31m ERROR: Recalibrated BAM missing after run: ${recal_bam} \e[0m"
            log_failure "Phase4_ApplyBQSR" "${sample_id}" "OUTPUT_MISSING" "recal.bam missing or 0 bytes"
            exit 1
        fi

        # Defensively ensure BAM index is created for downstream HaplotypeCaller
        if [ ! -s "${recal_bam}.bai" ] && [ ! -s "${recal_bam%.bam}.bai" ]; then
            echo -e "\e[32m Building index for recalibrated BAM... \e[0m"
            conda run -n mapping samtools index -@ "${threads}" "${recal_bam}"
        fi
        echo -e "\e[32m Recalibrated BAM completed -> ${recal_bam} \e[0m"
    fi

    #######################################################
    # PHASE 5: VARIANT CALLING & HARD FILTRATION          #
    #######################################################

    echo -e "\e[33m ====================================================== \e[0m"
    echo -e "\e[33m PHASE 5: GERMLINE VARIANT CALLING & HARD FILTRATION    \e[0m"
    echo -e "\e[33m ====================================================== \e[0m"

    raw_vcf="${CALLING_SAMPLE}/${sample_id}.raw_variants.vcf"
    raw_snps_vcf="${CALLING_SAMPLE}/${sample_id}.raw_snps.vcf"
    raw_indels_vcf="${CALLING_SAMPLE}/${sample_id}.raw_indels.vcf"
    filtered_snps_vcf="${CALLING_SAMPLE}/${sample_id}.filtered_snps.vcf"
    filtered_indels_vcf="${CALLING_SAMPLE}/${sample_id}.filtered_indels.vcf"

    # -- 5a. GATK HaplotypeCaller ----------------------------------------------
    echo -e "\e[31m ==================================== \e[0m"
    echo -e "\e[31m GATK HAPLOTYPECALLER: ${sample_id}   \e[0m"
    echo -e "\e[31m ==================================== \e[0m"

    if [ -s "${raw_vcf}" ]; then
        echo -e "\e[32m Raw variants VCF already present — skipping HaplotypeCaller \e[0m"
        log_failure "Phase5_HaplotypeCaller" "${sample_id}" "SKIPPED_EXISTS" "raw_variants.vcf exists"
    else
        if [ ! -s "${recal_bam}" ]; then
            echo -e "\e[31m ERROR: Recalibrated BAM missing for variant calling: ${recal_bam} \e[0m"
            log_failure "Phase5_HaplotypeCaller" "${sample_id}" "INPUT_MISSING" "recal.bam missing"
            exit 1
        fi

        echo -e "\e[32m Calling germline variants with GATK HaplotypeCaller... \e[0m"
        if ! conda run -n Hcalling gatk HaplotypeCaller \
            -R "${humanref}" \
            -I "${recal_bam}" \
            --native-pair-hmm-threads "${threads}" \
            --output-mode EMIT_VARIANTS_ONLY \
            --standard-min-confidence-threshold-for-calling 30 \
            -O "${raw_vcf}"; then
            echo -e "\e[31m ERROR: GATK HaplotypeCaller execution failed \e[0m"
            log_failure "Phase5_HaplotypeCaller" "${sample_id}" "EXECUTION_FAILED" "HaplotypeCaller non-zero exit code"
            exit 1
        fi

        if [ ! -s "${raw_vcf}" ]; then
            echo -e "\e[31m ERROR: Raw variants VCF missing after run: ${raw_vcf} \e[0m"
            log_failure "Phase5_HaplotypeCaller" "${sample_id}" "OUTPUT_MISSING" "raw_vcf missing or 0 bytes"
            exit 1
        fi
        echo -e "\e[32m Raw variant calling completed -> ${raw_vcf} \e[0m"
    fi

    # -- 5b. GATK SelectVariants (SNPs) ----------------------------------------
    echo -e "\e[31m ==================================== \e[0m"
    echo -e "\e[31m GATK SELECTVARIANTS (SNPs): ${sample_id} \e[0m"
    echo -e "\e[31m ==================================== \e[0m"

    if [ -s "${raw_snps_vcf}" ]; then
        echo -e "\e[32m Raw SNPs VCF already exists — skipping \e[0m"
        log_failure "Phase5_SelectVariants_SNPs" "${sample_id}" "SKIPPED_EXISTS" "raw_snps.vcf exists"
    else
        if [ ! -s "${raw_vcf}" ]; then
            echo -e "\e[31m ERROR: Raw variants VCF missing: ${raw_vcf} \e[0m"
            log_failure "Phase5_SelectVariants_SNPs" "${sample_id}" "INPUT_MISSING" "raw_vcf missing"
            exit 1
        fi

        echo -e "\e[32m Extracting SNPs with SelectVariants... \e[0m"
        if ! conda run -n Hcalling gatk SelectVariants \
            -R "${humanref}" \
            -V "${raw_vcf}" \
            --select-type-to-include SNP \
            -O "${raw_snps_vcf}"; then
            echo -e "\e[31m ERROR: SelectVariants SNPs failed \e[0m"
            log_failure "Phase5_SelectVariants_SNPs" "${sample_id}" "EXECUTION_FAILED" "SelectVariants SNPs failed"
            exit 1
        fi

        if [ ! -s "${raw_snps_vcf}" ]; then
            echo -e "\e[31m ERROR: Raw SNPs VCF missing after run: ${raw_snps_vcf} \e[0m"
            log_failure "Phase5_SelectVariants_SNPs" "${sample_id}" "OUTPUT_MISSING" "raw_snps_vcf missing or 0 bytes"
            exit 1
        fi
        echo -e "\e[32m Raw SNPs extracted -> ${raw_snps_vcf} \e[0m"
    fi

    # -- 5c. GATK SelectVariants (INDELs) --------------------------------------
    echo -e "\e[31m ====================================== \e[0m"
    echo -e "\e[31m GATK SELECTVARIANTS (INDELs): ${sample_id} \e[0m"
    echo -e "\e[31m ====================================== \e[0m"

    if [ -s "${raw_indels_vcf}" ]; then
        echo -e "\e[32m Raw INDELs VCF already exists — skipping \e[0m"
        log_failure "Phase5_SelectVariants_INDELs" "${sample_id}" "SKIPPED_EXISTS" "raw_indels.vcf exists"
    else
        if [ ! -s "${raw_vcf}" ]; then
            echo -e "\e[31m ERROR: Raw variants VCF missing: ${raw_vcf} \e[0m"
            log_failure "Phase5_SelectVariants_INDELs" "${sample_id}" "INPUT_MISSING" "raw_vcf missing"
            exit 1
        fi

        echo -e "\e[32m Extracting INDELs with SelectVariants... \e[0m"
        if ! conda run -n Hcalling gatk SelectVariants \
            -R "${humanref}" \
            -V "${raw_vcf}" \
            --select-type-to-include INDEL \
            -O "${raw_indels_vcf}"; then
            echo -e "\e[31m ERROR: SelectVariants INDELs failed \e[0m"
            log_failure "Phase5_SelectVariants_INDELs" "${sample_id}" "EXECUTION_FAILED" "SelectVariants INDELs failed"
            exit 1
        fi

        if [ ! -s "${raw_indels_vcf}" ]; then
            echo -e "\e[31m ERROR: Raw INDELs VCF missing after run: ${raw_indels_vcf} \e[0m"
            log_failure "Phase5_SelectVariants_INDELs" "${sample_id}" "OUTPUT_MISSING" "raw_indels_vcf missing or 0 bytes"
            exit 1
        fi
        echo -e "\e[32m Raw INDELs extracted -> ${raw_indels_vcf} \e[0m"
    fi

    # -- 5d. GATK VariantFiltration (SNPs) -------------------------------------
    echo -e "\e[31m ===================================== \e[0m"
    echo -e "\e[31m GATK VARIANTFILTRATION (SNPs): ${sample_id}\e[0m"
    echo -e "\e[31m ===================================== \e[0m"

    if [ -s "${filtered_snps_vcf}" ]; then
        echo -e "\e[32m Filtered SNPs VCF already exists — skipping \e[0m"
        log_failure "Phase5_VariantFiltration_SNPs" "${sample_id}" "SKIPPED_EXISTS" "filtered_snps.vcf exists"
    else
        if [ ! -s "${raw_snps_vcf}" ]; then
            echo -e "\e[31m ERROR: Raw SNPs VCF missing for filtration: ${raw_snps_vcf} \e[0m"
            log_failure "Phase5_VariantFiltration_SNPs" "${sample_id}" "INPUT_MISSING" "raw_snps.vcf missing"
            exit 1
        fi

        echo -e "\e[32m Applying hard filters to SNPs... \e[0m"
        if ! conda run -n Hcalling gatk VariantFiltration \
            -R "${humanref}" \
            -V "${raw_snps_vcf}" \
            --filter-expression "QD < 2.0" --filter-name "QualByDepth" \
            --filter-expression "FS > 60.0" --filter-name "FisherStrand" \
            --filter-expression "MQ < 40.0" --filter-name "RMSMappingQuality" \
            --filter-expression "HaplotypeScore > 13.0" --filter-name "HaplotypeScore" \
            --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSumTest" \
            --filter-expression "QUAL < 30.0 || DP < 6 || DP > 5000 || HRun > 5" --filter-name "StandardFilters" \
            --filter-expression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" --filter-name "HARD_TO_VALIDATE" \
            -O "${filtered_snps_vcf}"; then
            echo -e "\e[31m ERROR: VariantFiltration SNPs failed \e[0m"
            log_failure "Phase5_VariantFiltration_SNPs" "${sample_id}" "EXECUTION_FAILED" "VariantFiltration SNPs non-zero exit"
            exit 1
        fi

        if [ ! -s "${filtered_snps_vcf}" ]; then
            echo -e "\e[31m ERROR: Filtered SNPs VCF missing after run: ${filtered_snps_vcf} \e[0m"
            log_failure "Phase5_VariantFiltration_SNPs" "${sample_id}" "OUTPUT_MISSING" "filtered_snps.vcf missing or 0 bytes"
            exit 1
        fi
        echo -e "\e[32m Filtered SNPs completed -> ${filtered_snps_vcf} \e[0m"
    fi

    # -- 5e. GATK VariantFiltration (INDELs) -----------------------------------
    echo -e "\e[31m ======================================= \e[0m"
    echo -e "\e[31m GATK VARIANTFILTRATION (INDELs): ${sample_id}\e[0m"
    echo -e "\e[31m ======================================= \e[0m"

    if [ -s "${filtered_indels_vcf}" ]; then
        echo -e "\e[32m Filtered INDELs VCF already exists — skipping \e[0m"
        log_failure "Phase5_VariantFiltration_INDELs" "${sample_id}" "SKIPPED_EXISTS" "filtered_indels.vcf exists"
    else
        if [ ! -s "${raw_indels_vcf}" ]; then
            echo -e "\e[31m ERROR: Raw INDELs VCF missing for filtration: ${raw_indels_vcf} \e[0m"
            log_failure "Phase5_VariantFiltration_INDELs" "${sample_id}" "INPUT_MISSING" "raw_indels.vcf missing"
            exit 1
        fi

        echo -e "\e[32m Applying hard filters to INDELs... \e[0m"
        if ! conda run -n Hcalling gatk VariantFiltration \
            -R "${humanref}" \
            -V "${raw_indels_vcf}" \
            --filter-expression "QD < 2.0" --filter-name "QualByDepth" \
            --filter-expression "FS > 200.0" --filter-name "FisherStrand" \
            --filter-expression "MQ < 40.0" --filter-name "RMSMappingQuality" \
            --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSumTest" \
            --filter-expression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" --filter-name "HARD_TO_VALIDATE" \
            --filter-expression "QUAL < 30.0 || DP < 6 || DP > 5000 || HRun > 5" --filter-name "QualFilter" \
            -O "${filtered_indels_vcf}"; then
            echo -e "\e[31m ERROR: VariantFiltration INDELs failed \e[0m"
            log_failure "Phase5_VariantFiltration_INDELs" "${sample_id}" "EXECUTION_FAILED" "VariantFiltration INDELs non-zero exit"
            exit 1
        fi

        if [ ! -s "${filtered_indels_vcf}" ]; then
            echo -e "\e[31m ERROR: Filtered INDELs VCF missing after run: ${filtered_indels_vcf} \e[0m"
            log_failure "Phase5_VariantFiltration_INDELs" "${sample_id}" "OUTPUT_MISSING" "filtered_indels.vcf missing or 0 bytes"
            exit 1
        fi
        echo -e "\e[32m Filtered INDELs completed -> ${filtered_indels_vcf} \e[0m"
    fi

    #############################################
    # PHASE 6: VARIANT FUNCTIONAL ANNOTATION    #
    #############################################

    echo -e "\e[33m ====================================================== \e[0m"
    echo -e "\e[33m PHASE 6: FUNCTIONAL ANNOTATION (ANNOVAR)              \e[0m"
    echo -e "\e[33m ====================================================== \e[0m"

    # -- 6a. ANNOVAR Database Preparation --------------------------------------
    echo -e "\e[31m =============================== \e[0m"
    echo -e "\e[31m ANNOVAR DATABASE PREPARATION    \e[0m"
    echo -e "\e[31m =============================== \e[0m"

    IFS=',' read -r -a annovar_db_array <<< "${annovar_protocol}"
    for db in "${annovar_db_array[@]}"; do
        db_marker="${humandb}/${annovar_buildver}_${db}.txt"

        if [ -s "${db_marker}" ]; then
            echo -e "\e[32m ANNOVAR DB [${db}] already downloaded — skipping \e[0m"
            log_failure "Phase6_ANNOVAR_DB" "${db}" "SKIPPED_EXISTS" "DB marker ${db_marker} exists"
            continue
        fi

        if [ ! -f "${annovar_downdb}" ]; then
            echo -e "\e[31m WARNING: ANNOVAR annotate_variation.pl missing at ${annovar_downdb} \e[0m"
            log_failure "Phase6_ANNOVAR_DB" "${db}" "INPUT_MISSING" "Script ${annovar_downdb} missing"
            continue
        fi

        echo -e "\e[32m Downloading ANNOVAR database: ${db}... \e[0m"
        if [ "${db}" == "cytoBand" ]; then
            if ! perl "${annovar_downdb}" \
                -buildver "${annovar_buildver}" \
                -downdb "${db}" \
                "${humandb}/"; then
                echo -e "\e[31m WARNING: Download failed for ${db} \e[0m"
                log_failure "Phase6_ANNOVAR_DB" "${db}" "EXECUTION_FAILED" "Download of ${db} failed"
            fi
        else
            if ! perl "${annovar_downdb}" \
                -buildver "${annovar_buildver}" \
                -downdb -webfrom annovar "${db}" \
                "${humandb}/"; then
                echo -e "\e[31m WARNING: Download failed for ${db} \e[0m"
                log_failure "Phase6_ANNOVAR_DB" "${db}" "EXECUTION_FAILED" "Download of ${db} failed"
            fi
        fi
    done

    # -- 6b. VCF to ANNOVAR avinput format conversion --------------------------
    echo -e "\e[31m ============================== \e[0m"
    echo -e "\e[31m CONVERT2ANNOVAR: ${sample_id}  \e[0m"
    echo -e "\e[31m ============================== \e[0m"

    snps_avinput="${ANNOTATION_SAMPLE}/${sample_id}.filtered_snps.avinput"
    indels_avinput="${ANNOTATION_SAMPLE}/${sample_id}.filtered_indels.avinput"

    # Conversion for SNPs
    if [ -s "${snps_avinput}" ]; then
        echo -e "\e[32m SNPs avinput already exists — skipping conversion \e[0m"
        log_failure "Phase6_Convert_SNPs" "${sample_id}" "SKIPPED_EXISTS" "snps.avinput exists"
    else
        if [ ! -s "${filtered_snps_vcf}" ]; then
            echo -e "\e[31m ERROR: Filtered SNPs VCF missing for ANNOVAR conversion: ${filtered_snps_vcf} \e[0m"
            log_failure "Phase6_Convert_SNPs" "${sample_id}" "INPUT_MISSING" "filtered_snps.vcf missing"
            exit 1
        fi

        echo -e "\e[32m Converting filtered SNPs VCF to ANNOVAR avinput... \e[0m"
        if ! perl "${annovar_convert}" \
            -format vcf4 --includeinfo \
            "${filtered_snps_vcf}" -outfile "${snps_avinput}"; then
            echo -e "\e[31m ERROR: convert2annovar failed on SNPs \e[0m"
            log_failure "Phase6_Convert_SNPs" "${sample_id}" "EXECUTION_FAILED" "convert2annovar SNPs failed"
            exit 1
        fi

        if [ ! -s "${snps_avinput}" ]; then
            echo -e "\e[31m ERROR: SNPs avinput missing after conversion: ${snps_avinput} \e[0m"
            log_failure "Phase6_Convert_SNPs" "${sample_id}" "OUTPUT_MISSING" "snps_avinput missing or 0 bytes"
            exit 1
        fi
        echo -e "\e[32m SNPs converted -> ${snps_avinput} \e[0m"
    fi

    # Conversion for INDELs
    if [ -s "${indels_avinput}" ]; then
        echo -e "\e[32m INDELs avinput already exists — skipping conversion \e[0m"
        log_failure "Phase6_Convert_INDELs" "${sample_id}" "SKIPPED_EXISTS" "indels.avinput exists"
    else
        if [ ! -s "${filtered_indels_vcf}" ]; then
            echo -e "\e[31m ERROR: Filtered INDELs VCF missing for ANNOVAR conversion: ${filtered_indels_vcf} \e[0m"
            log_failure "Phase6_Convert_INDELs" "${sample_id}" "INPUT_MISSING" "filtered_indels.vcf missing"
            exit 1
        fi

        echo -e "\e[32m Converting filtered INDELs VCF to ANNOVAR avinput... \e[0m"
        if ! perl "${annovar_convert}" \
            -format vcf4 --includeinfo \
            "${filtered_indels_vcf}" -outfile "${indels_avinput}"; then
            echo -e "\e[31m ERROR: convert2annovar failed on INDELs \e[0m"
            log_failure "Phase6_Convert_INDELs" "${sample_id}" "EXECUTION_FAILED" "convert2annovar INDELs failed"
            exit 1
        fi

        if [ ! -s "${indels_avinput}" ]; then
            echo -e "\e[31m ERROR: INDELs avinput missing after conversion: ${indels_avinput} \e[0m"
            log_failure "Phase6_Convert_INDELs" "${sample_id}" "OUTPUT_MISSING" "indels_avinput missing or 0 bytes"
            exit 1
        fi
        echo -e "\e[32m INDELs converted -> ${indels_avinput} \e[0m"
    fi

    # -- 6c. Functional Annotation via table_annovar ----------------------------
    echo -e "\e[31m ============================== \e[0m"
    echo -e "\e[31m TABLE_ANNOVAR: ${sample_id}    \e[0m"
    echo -e "\e[31m ============================== \e[0m"

    annotated_snps_csv="${ANNOTATION_SAMPLE}/${sample_id}.annotated_snps.${annovar_buildver}_multianno.csv"
    annotated_indels_csv="${ANNOTATION_SAMPLE}/${sample_id}.annotated_indels.${annovar_buildver}_multianno.csv"

    # Annotation of SNPs
    if [ -s "${annotated_snps_csv}" ]; then
        echo -e "\e[32m Annotated SNPs CSV already exists — skipping table_annovar \e[0m"
        log_failure "Phase6_TableAnnovar_SNPs" "${sample_id}" "SKIPPED_EXISTS" "annotated_snps multianno.csv exists"
    else
        if [ ! -s "${snps_avinput}" ]; then
            echo -e "\e[31m ERROR: SNPs avinput missing for table_annovar: ${snps_avinput} \e[0m"
            log_failure "Phase6_TableAnnovar_SNPs" "${sample_id}" "INPUT_MISSING" "snps.avinput missing"
            exit 1
        fi

        echo -e "\e[32m Annotating SNPs with table_annovar.pl... \e[0m"
        if ! perl "${annovar_table}" \
            "${snps_avinput}" \
            "${humandb}/" \
            -buildver "${annovar_buildver}" \
            -out "${ANNOTATION_SAMPLE}/${sample_id}.annotated_snps" \
            -remove \
            -protocol "${annovar_protocol}" \
            -operation "${annovar_operation}" \
            -nastring . \
            -csvout \
            -polish \
            --otherinfo; then
            echo -e "\e[31m ERROR: table_annovar execution failed on SNPs \e[0m"
            log_failure "Phase6_TableAnnovar_SNPs" "${sample_id}" "EXECUTION_FAILED" "table_annovar SNPs failed"
            exit 1
        fi

        if [ ! -s "${annotated_snps_csv}" ]; then
            echo -e "\e[31m ERROR: Annotated SNPs CSV missing after run: ${annotated_snps_csv} \e[0m"
            log_failure "Phase6_TableAnnovar_SNPs" "${sample_id}" "OUTPUT_MISSING" "annotated_snps CSV missing"
            exit 1
        fi
        echo -e "\e[32m SNPs annotation completed -> ${annotated_snps_csv} \e[0m"
    fi

    # Annotation of INDELs
    if [ -s "${annotated_indels_csv}" ]; then
        echo -e "\e[32m Annotated INDELs CSV already exists — skipping table_annovar \e[0m"
        log_failure "Phase6_TableAnnovar_INDELs" "${sample_id}" "SKIPPED_EXISTS" "annotated_indels multianno.csv exists"
    else
        if [ ! -s "${indels_avinput}" ]; then
            echo -e "\e[31m ERROR: INDELs avinput missing for table_annovar: ${indels_avinput} \e[0m"
            log_failure "Phase6_TableAnnovar_INDELs" "${sample_id}" "INPUT_MISSING" "indels.avinput missing"
            exit 1
        fi

        echo -e "\e[32m Annotating INDELs with table_annovar.pl... \e[0m"
        if ! perl "${annovar_table}" \
            "${indels_avinput}" \
            "${humandb}/" \
            -buildver "${annovar_buildver}" \
            -out "${ANNOTATION_SAMPLE}/${sample_id}.annotated_indels" \
            -remove \
            -protocol "${annovar_protocol}" \
            -operation "${annovar_operation}" \
            -nastring . \
            -csvout \
            -polish \
            --otherinfo; then
            echo -e "\e[31m ERROR: table_annovar execution failed on INDELs \e[0m"
            log_failure "Phase6_TableAnnovar_INDELs" "${sample_id}" "EXECUTION_FAILED" "table_annovar INDELs failed"
            exit 1
        fi

        if [ ! -s "${annotated_indels_csv}" ]; then
            echo -e "\e[31m ERROR: Annotated INDELs CSV missing after run: ${annotated_indels_csv} \e[0m"
            log_failure "Phase6_TableAnnovar_INDELs" "${sample_id}" "OUTPUT_MISSING" "annotated_indels CSV missing"
            exit 1
        fi
        echo -e "\e[32m INDELs annotation completed -> ${annotated_indels_csv} \e[0m"
    fi

########################################################
# END-OF-SCRIPT SUMMARY & CHECKPOINT REPORT            #
########################################################

echo ""
echo -e "\e[32m ======================================================== \e[0m"
echo -e "\e[32m PIPELINE COMPLETE — ${sample_id} — $(date)               \e[0m"
echo -e "\e[32m ======================================================== \e[0m"
echo ""
echo -e "\e[32m -- PREPROCESSING & QC ---------------------------------- \e[0m"
echo -e "\e[32m  Raw FastQC         : ${FASTQC_RAW}                      \e[0m"
echo -e "\e[32m  Trimmed reads (R1) : ${read1t}                          \e[0m"
echo -e "\e[32m  Trimmed reads (R2) : ${read2t}                          \e[0m"
echo -e "\e[32m  Trimmed FastQC     : ${FASTQC_PP}                       \e[0m"
echo ""
echo -e "\e[32m -- ALIGNMENT & BQSR ------------------------------------ \e[0m"
echo -e "\e[32m  Deduplicated BAM   : ${dedup_bam}                       \e[0m"
echo -e "\e[32m  BQSR recal table   : ${recal_table}                     \e[0m"
echo -e "\e[32m  Recalibrated BAM   : ${recal_bam}                       \e[0m"
echo ""
echo -e "\e[32m -- VARIANT CALLING ------------------------------------- \e[0m"
echo -e "\e[32m  Raw VCF            : ${raw_vcf}                         \e[0m"
echo -e "\e[32m  Filtered SNPs VCF  : ${filtered_snps_vcf}               \e[0m"
echo -e "\e[32m  Filtered INDELs VCF: ${filtered_indels_vcf}             \e[0m"
echo ""
echo -e "\e[32m -- FUNCTIONAL ANNOTATION ------------------------------- \e[0m"
echo -e "\e[32m  SNPs Multianno CSV : ${annotated_snps_csv}              \e[0m"
echo -e "\e[32m  INDELs Multianno CSV: ${annotated_indels_csv}           \e[0m"
echo ""
echo -e "\e[32m -- LOG FILES ------------------------------------------- \e[0m"
echo -e "\e[32m  Full execution log : ${LOG}                             \e[0m"

fail_count=$(awk '/\[FAILED\]|\[EXECUTION_FAILED\]|\[OUTPUT_MISSING/ {c++} END {print c+0}' "${FAIL_LOG}")
skip_count=$(awk '/\[SKIPPED/ {c++} END {print c+0}' "${FAIL_LOG}")

if [ "${fail_count}" -gt 0 ]; then
    echo -e "\e[31m  Failure/Issues log : ${FAIL_LOG} (${fail_count} failures detected!) \e[0m"
    echo -e "\e[31m  >>> Inspect ${FAIL_LOG} to see which steps failed and why. \e[0m"
else
    echo -e "\e[32m  Failure/Issues log : ${FAIL_LOG} (0 errors recorded)     \e[0m"
fi
echo -e "\e[32m  Skipped checkpoints: ${skip_count} records \e[0m"
echo ""
echo -e "\e[32m -- SUGGESTED NEXT STEPS -------------------------------- \e[0m"
echo -e "\e[32m  1. Inspect ANNOVAR CSV annotations in R/Python for high-impact mutations (frameshifts, stopgain, nonsynonymous). \e[0m"
echo -e "\e[32m  2. Cross-reference candidate variants with ClinVar or ACMG pathogenicity classifications. \e[0m"
echo -e "\e[32m  3. Verify coverage and read depth distributions across exons using MultiQC or Mosdepth. \e[0m"
echo -e "\e[32m  Full log: ${LOG}                                        \e[0m"
echo -e "\e[32m ======================================================== \e[0m"