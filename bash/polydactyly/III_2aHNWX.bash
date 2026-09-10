#!/bin/bash
#script for whole genome sequencing: preprocessing, mapping, BQSR, variant calling, and annotation
#Building No.2
set -euo pipefail
################
#GLOBAL SETTING#
################

REF_PATH="/storage/student9/references"
TOOLS_PATH="/storage/student9/tools"
SAMPLE_PATH="/storage/student9/projects/III_2aHNWX"
threads=16

#sample id can be passed as the first argument, e.g. ./III_2aHNWX.bash III_2aHNWX
sample_id="${1:-III_2aHNWX}"

#path to references and databases
humanref="${REF_PATH}/hg38.fa"
humanref_dict="${REF_PATH}/hg38.dict"
truseq3="${REF_PATH}/TruSeq3-PE.fa"
dbSNP="${REF_PATH}/common_all_20180418.vcf"

#ANNOVAR setup (not on bioconda - see conda_envs.txt - called directly via perl)
annovar_path="${TOOLS_PATH}/annovar"
annovar_downdb="${annovar_path}/annotate_variation.pl"
annovar_convert="${annovar_path}/convert2annovar.pl"
annovar_table="${annovar_path}/table_annovar.pl"
humandb="${REF_PATH}/humandb"
annovar_buildver="hg38"
annovar_protocol="refGene,cytoBand,exac03,avsnp150,dbnsfp30a"
annovar_operation="g,r,f,f,f"

#create working directories
RAW_WGS="${SAMPLE_PATH}/raw"
PREPROCESSING_WGS="${SAMPLE_PATH}/preprocessing_WGS"
MAPPING_WGS="${SAMPLE_PATH}/mapping_WGS"
CALLING_WGS="${SAMPLE_PATH}/variant_calling_WGS"
ANNOTATION_WGS="${SAMPLE_PATH}/annotation_WGS"

mkdir -p "${SAMPLE_PATH}/logs"
mkdir -p "${PREPROCESSING_WGS}/fastqc_raw"
mkdir -p "${PREPROCESSING_WGS}/fastqc_pp"
mkdir -p "${MAPPING_WGS}/${sample_id}"
mkdir -p "${CALLING_WGS}/${sample_id}"
mkdir -p "${ANNOTATION_WGS}/${sample_id}"
mkdir -p "${humandb}"

# Dual-log setup: main run log and dedicated failure/skip tracking log
LOG="${SAMPLE_PATH}/logs/WGS_${sample_id}.log"
FAIL_LOG="${SAMPLE_PATH}/logs/WGS_${sample_id}.failed_skipped.log"
exec > >(tee -a "${LOG}") 2>&1

echo "======================================================" >> "${FAIL_LOG}"
echo " Failure & Skip Log — Started: $(date)" >> "${FAIL_LOG}"
echo "======================================================" >> "${FAIL_LOG}"

log_failure() {
    local phase="$1"
    local sample="$2"
    local status="$3" # e.g. "SKIPPED_EXISTS", "INPUT_MISSING", "EXECUTION_FAILED", "OUTPUT_MISSING"
    local reason="$4"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [${phase}] [${sample}] [${status}] ${reason}" >> "${FAIL_LOG}"
}

echo "============================================================"
echo " WGS variant calling: ${sample_id} - Started: $(date)"
echo " Full log           : ${LOG}"
echo " Failure & skip log : ${FAIL_LOG}"
echo "============================================================"

	#############################################
	#ONE-TIME REFERENCE SETUP (skipped if present)#
	#############################################

	echo -e "\e[31m ======================= \e[0m"
	echo -e "\e[31m REFERENCE PREPARATION \e[0m"
	echo -e "\e[31m ======================= \e[0m"

	#BWA index - only needs to be built once for the reference genome
	if [ ! -f "${humanref}.bwt" ]; then
		echo -e "\e[32m Building BWA index for ${humanref}... \e[0m"
		conda run -n mapping bwa index "${humanref}"
	else
		echo -e "\e[32m BWA index already present - skipping \e[0m"
	fi

	#Picard sequence dictionary - only needs to be built once for the reference genome
	if [ ! -f "${humanref_dict}" ]; then
		echo -e "\e[32m Creating sequence dictionary for ${humanref}... \e[0m"
		conda run -n mapping picard CreateSequenceDictionary \
		REFERENCE="${humanref}" \
		OUTPUT="${humanref_dict}"
	else
		echo -e "\e[32m Sequence dictionary already present - skipping \e[0m"
	fi

	#Reference fasta index (required by GATK/samtools) - only needs to be built once
	if [ ! -f "${humanref}.fai" ]; then
		echo -e "\e[32m Indexing ${humanref} with samtools faidx... \e[0m"
		conda run -n mapping samtools faidx "${humanref}"
	else
		echo -e "\e[32m FASTA index already present - skipping \e[0m"
	fi

	#########################
	#SAMPLE DATA PREPARATION#
	#########################

	read1="${RAW_WGS}/${sample_id}_1.fastq.gz"
	read2="${RAW_WGS}/${sample_id}_2.fastq.gz"

	echo -e "\e[31m ========================== \e[0m"
	echo -e "\e[31m FASTQC: ${sample_id} (RAW) \e[0m"
	echo -e "\e[31m ========================== \e[0m"

	#review raw fastq
	conda run -n preprocessing fastqc \
	--threads ${threads} \
	--outdir "${PREPROCESSING_WGS}/fastqc_raw" \
	"${read1}" "${read2}"

	echo -e "\e[31m ========================= \e[0m"
	echo -e "\e[31m TRIMMOMATIC: ${sample_id} \e[0m"
	echo -e "\e[31m ========================= \e[0m"

	trim_output_1="${PREPROCESSING_WGS}/${sample_id}.R1"
	trim_output_2="${PREPROCESSING_WGS}/${sample_id}.R2"

	#fastq data manipulation
	conda run -n preprocessing trimmomatic PE \
	-threads ${threads} \
	-phred33 \
	"${read1}" "${read2}" \
	"${trim_output_1}.paired.fastq.gz" "${trim_output_1}.unpaired.fastq.gz" \
	"${trim_output_2}.paired.fastq.gz" "${trim_output_2}.unpaired.fastq.gz" \
	ILLUMINACLIP:"${truseq3}":2:30:10:8:true \
	HEADCROP:3 TRAILING:10 MINLEN:25

	read1t="${trim_output_1}.paired.fastq.gz"
	read2t="${trim_output_2}.paired.fastq.gz"

	echo -e "\e[31m ========================= \e[0m"
	echo -e "\e[31m FASTQC: ${sample_id} (PP) \e[0m"
	echo -e "\e[31m ========================= \e[0m"

	#review trimmed fastq
	conda run -n preprocessing fastqc \
	--threads ${threads} \
	--outdir "${PREPROCESSING_WGS}/fastqc_pp" \
	"${read1t}" "${read2t}"

	##################
	#MAPPING WITH BWA#
	##################

	MAPPING_SAMPLE="${MAPPING_WGS}/${sample_id}"
	sam_file="${MAPPING_SAMPLE}/${sample_id}.sam"
	sorted_bam="${MAPPING_SAMPLE}/${sample_id}.sorted_coord.bam"
	dedup_bam="${MAPPING_SAMPLE}/${sample_id}.dedup.bam"
	dedup_metrics="${MAPPING_SAMPLE}/${sample_id}.dedup_metrics.txt"

	echo -e "\e[31m ================= \e[0m"
	echo -e "\e[31m BWA MEM: ${sample_id} \e[0m"
	echo -e "\e[31m ================= \e[0m"

	#BUG FIX: original used a "-p" (interleaved-input) flag while passing two
	#separate paired FASTQ files, and the read group used "BL:lib1" instead of
	#the valid "LB:lib1" tag. Both are fixed below; RG is built per sample_id.
	conda run -n mapping bwa mem \
	-t ${threads} \
	-R "@RG\tID:${sample_id}\tSM:${sample_id}\tPL:illumina\tLB:lib1\tPU:unit1" \
	"${humanref}" \
	"${read1t}" "${read2t}" \
	> "${sam_file}"

	echo -e "\e[31m ================================= \e[0m"
	echo -e "\e[31m SAMTOOLS SORT & INDEX: ${sample_id} \e[0m"
	echo -e "\e[31m ================================= \e[0m"

	#conversion + coordinate sort in one step (replaces Picard SortSam)
	conda run -n mapping samtools sort \
	-@ ${threads} \
	-o "${sorted_bam}" \
	"${sam_file}"

	rm "${sam_file}"

	echo -e "\e[31m =========================== \e[0m"
	echo -e "\e[31m PICARD MARKDUPLICATES: ${sample_id} \e[0m"
	echo -e "\e[31m =========================== \e[0m"

	conda run -n mapping picard MarkDuplicates \
	INPUT="${sorted_bam}" \
	OUTPUT="${dedup_bam}" \
	METRICS_FILE="${dedup_metrics}" \
	REMOVE_DUPLICATES=true \
	CREATE_INDEX=true \
	VALIDATION_STRINGENCY=SILENT

	echo -e "\e[31m ======================== \e[0m"
	echo -e "\e[31m SAMTOOLS INDEX: ${sample_id} \e[0m"
	echo -e "\e[31m ======================== \e[0m"

	#BUG FIX: original indexed the .bam but BaseRecalibrator was later pointed
	#at the .bam.bai file instead of the .bam itself - CREATE_INDEX=true above
	#already produces dedup.bai, this confirms/re-indexes defensively
	conda run -n mapping samtools index \
	-@ ${threads} \
	"${dedup_bam}"

	#########################################
	#BASE QUALITY SCORE RECALIBRATION (BQSR)#
	#########################################

	recal_table="${MAPPING_SAMPLE}/${sample_id}.recal_data.grp"
	recal_bam="${MAPPING_SAMPLE}/${sample_id}.recal.bam"

	echo -e "\e[31m =============================== \e[0m"
	echo -e "\e[31m GATK BASERECALIBRATOR: ${sample_id} \e[0m"
	echo -e "\e[31m =============================== \e[0m"

	conda run -n Hcalling gatk BaseRecalibrator \
	-R "${humanref}" \
	-I "${dedup_bam}" \
	--known-sites "${dbSNP}" \
	-O "${recal_table}"

	echo -e "\e[31m ======================== \e[0m"
	echo -e "\e[31m GATK APPLYBQSR: ${sample_id} \e[0m"
	echo -e "\e[31m ======================== \e[0m"

	#BUG FIX: the original pipeline generated a recalibration table but never
	#applied it - HaplotypeCaller ran on the raw dedup.bam. ApplyBQSR is added
	#here so the recalibrated BAM is what actually gets called on below.
	conda run -n Hcalling gatk ApplyBQSR \
	-R "${humanref}" \
	-I "${dedup_bam}" \
	--bqsr-recal-file "${recal_table}" \
	-O "${recal_bam}"

	##################
	#CALLING VARIANTS#
	##################

	CALLING_SAMPLE="${CALLING_WGS}/${sample_id}"
	raw_vcf="${CALLING_SAMPLE}/${sample_id}.raw_variants.vcf"

	echo -e "\e[31m ==================================== \e[0m"
	echo -e "\e[31m GATK HAPLOTYPECALLER: ${sample_id} \e[0m"
	echo -e "\e[31m ==================================== \e[0m"

	conda run -n Hcalling gatk HaplotypeCaller \
	-R "${humanref}" \
	-I "${recal_bam}" \
	--output-mode EMIT_VARIANTS_ONLY \
	--standard-min-confidence-threshold-for-calling 30 \
	-O "${raw_vcf}"

	echo -e "\e[31m ============================ \e[0m"
	echo -e "\e[31m GATK SELECTVARIANTS: ${sample_id} \e[0m"
	echo -e "\e[31m ============================ \e[0m"

	raw_snps_vcf="${CALLING_SAMPLE}/${sample_id}.raw_snps.vcf"
	raw_indels_vcf="${CALLING_SAMPLE}/${sample_id}.raw_indels.vcf"

	conda run -n Hcalling gatk SelectVariants \
	-R "${humanref}" \
	-V "${raw_vcf}" \
	--select-type-to-include SNP \
	-O "${raw_snps_vcf}"

	conda run -n Hcalling gatk SelectVariants \
	-R "${humanref}" \
	-V "${raw_vcf}" \
	--select-type-to-include INDEL \
	-O "${raw_indels_vcf}"

	####################
	#VARIANT FILTRATION#
	####################

	echo -e "\e[31m ========================== \e[0m"
	echo -e "\e[31m GATK VARIANT FILTRATION: ${sample_id} \e[0m"
	echo -e "\e[31m ========================== \e[0m"

	filtered_snps_vcf="${CALLING_SAMPLE}/${sample_id}.filtered_snps.vcf"
	filtered_indels_vcf="${CALLING_SAMPLE}/${sample_id}.filtered_indels.vcf"

	#NOTE: HaplotypeScore and MQ0 are legacy GATK3 annotations. If GATK 4.6.x
	#does not emit them, the corresponding filter will simply never trigger
	#(harmless) rather than error - kept here to match the original protocol.
	conda run -n Hcalling gatk VariantFiltration \
	-R "${humanref}" \
	-V "${raw_snps_vcf}" \
	--filter-expression "QD < 2.0" --filter-name "QualByDepth" \
	--filter-expression "FS > 60.0" --filter-name "FisherStrand" \
	--filter-expression "MQ < 40.0" --filter-name "RMSMappingQuality" \
	--filter-expression "HaplotypeScore > 13.0" --filter-name "HaplotypeScore" \
	--filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSumTest" \
	--filter-expression "QUAL < 30.0 || DP < 6 || DP > 5000 || HRun > 5" --filter-name "StandardFilters" \
	--filter-expression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" --filter-name "HARD_TO_VALIDATE" \
	-O "${filtered_snps_vcf}"

	conda run -n Hcalling gatk VariantFiltration \
	-R "${humanref}" \
	-V "${raw_indels_vcf}" \
	--filter-expression "QD < 2.0" --filter-name "QualByDepth" \
	--filter-expression "FS > 200.0" --filter-name "FisherStrand" \
	--filter-expression "MQ < 40.0" --filter-name "RMSMappingQuality" \
	--filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSumTest" \
	--filter-expression "MQ0 >= 4 && ((MQ0 / (1.0*DP)) > 0.1)" --filter-name "HARD_TO_VALIDATE" \
	--filter-expression "QUAL < 30.0 || DP < 6 || DP > 5000 || HRun > 5" --filter-name "QualFilter" \
	-O "${filtered_indels_vcf}"

	echo -e "\e[32m Variant calling complete for ${sample_id} \e[0m"
	echo -e "\e[32m Filtered SNPs  : ${filtered_snps_vcf} \e[0m"
	echo -e "\e[32m Filtered indels: ${filtered_indels_vcf} \e[0m"

	############################################
	#ANNOVAR DATABASE SETUP (skipped if present)#
	############################################

	echo -e "\e[31m =============================== \e[0m"
	echo -e "\e[31m ANNOVAR DATABASE PREPARATION \e[0m"
	echo -e "\e[31m =============================== \e[0m"

	IFS=',' read -r -a annovar_db_array <<< "${annovar_protocol}"
for db in "${annovar_db_array[@]}"; do
	#refGene/exac03/avsnp150/dbnsfp30a are pulled from the annovar webserver;
	#cytoBand is pulled from the UCSC mirror - both use -downdb, only the
	#-webfrom flag differs, so cytoBand is special-cased here
	db_marker="${humandb}/${annovar_buildver}_${db}.txt"

	if [ -f "${db_marker}" ]; then
		echo -e "\e[32m ${db}: already downloaded - skipping \e[0m"
		continue
	fi

	echo -e "\e[32m Downloading ANNOVAR database: ${db}... \e[0m"

	if [ "${db}" == "cytoBand" ]; then
		perl "${annovar_downdb}" \
		-buildver "${annovar_buildver}" \
		-downdb "${db}" \
		"${humandb}/"
	else
		perl "${annovar_downdb}" \
		-buildver "${annovar_buildver}" \
		-downdb -webfrom annovar "${db}" \
		"${humandb}/"
	fi
done

	############
	#ANNOTATION#
	############

	ANNOTATION_SAMPLE="${ANNOTATION_WGS}/${sample_id}"
	snps_avinput="${ANNOTATION_SAMPLE}/${sample_id}.filtered_snps.avinput"
	indels_avinput="${ANNOTATION_SAMPLE}/${sample_id}.filtered_indels.avinput"

	echo -e "\e[31m ============================== \e[0m"
	echo -e "\e[31m PREPARE INPUT FILE FOR ANNOVAR: ${sample_id} \e[0m"
	echo -e "\e[31m ============================== \e[0m"

	#--includeinfo carries genotype/depth fields from the VCF through to the
	#annotated table - the original commands dropped this for both variant types
	perl "${annovar_convert}" \
	-format vcf4 --includeinfo \
	"${filtered_snps_vcf}" -outfile "${snps_avinput}"

	perl "${annovar_convert}" \
	-format vcf4 --includeinfo \
	"${filtered_indels_vcf}" -outfile "${indels_avinput}"

	echo -e "\e[31m ================== \e[0m"
	echo -e "\e[31m VARIANT ANNOTATION: ${sample_id} \e[0m"
	echo -e "\e[31m ================== \e[0m"

	perl "${annovar_table}" \
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
	--otherinfo

	perl "${annovar_table}" \
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
	--otherinfo

	echo -e "\e[32m Annotation complete for ${sample_id} \e[0m"

echo ""
echo "========================================================"
echo " PIPELINE COMPLETE - ${sample_id} - $(date)"
echo "========================================================"
echo "  Trimmed reads      : ${read1t} , ${read2t}"
echo "  Deduplicated BAM   : ${dedup_bam}"
echo "  BQSR recal table   : ${recal_table}"
echo "  Recalibrated BAM   : ${recal_bam}"
echo "  Raw VCF            : ${raw_vcf}"
echo "  Filtered SNPs VCF  : ${filtered_snps_vcf}"
echo "  Filtered indels VCF: ${filtered_indels_vcf}"
echo "  SNPs annotation    : ${ANNOTATION_SAMPLE}/${sample_id}.annotated_snps.${annovar_buildver}_multianno.csv"
echo "  Indels annotation  : ${ANNOTATION_SAMPLE}/${sample_id}.annotated_indels.${annovar_buildver}_multianno.csv"
echo "  Full execution log : ${LOG}"

fail_count=$(grep -c '\[FAILED\]\|\[EXECUTION_FAILED\]\|\[OUTPUT_MISSING' "${FAIL_LOG}" 2>/dev/null || echo 0)
skip_count=$(grep -c '\[SKIPPED' "${FAIL_LOG}" 2>/dev/null || echo 0)

if [ "${fail_count}" -gt 0 ]; then
	echo -e "\e[31m  Failure/Issues log : ${FAIL_LOG} (${fail_count} failures detected!) \e[0m"
	echo -e "\e[31m  >>> Inspect ${FAIL_LOG} to see which steps failed and why. \e[0m"
else
	echo -e "\e[32m  Failure/Issues log : ${FAIL_LOG} (0 errors recorded) \e[0m"
fi
echo -e "\e[32m  Skipped checkpoints: ${skip_count} records \e[0m"
echo "========================================================"