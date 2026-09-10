#!/bin/bash
#HMH: combined WES pipeline for 7 syndactyly patients
#Preprocessing -> short-variant calling -> structural-variant calling -> CNV calling
#Building No.5
set -euo pipefail
################
#GLOBAL SETTING#
################

REF_PATH="/storage/student9/references"
TOOLS_PATH="/storage/student9/tools"
SAMPLE_PATH="/storage/student9/projects/syndactyly_samples"
threads=16

#set to "false" to keep large intermediate BAMs (SAM/sort/fixmate/group/filter stages)
cleanup_intermediates="true"

#7 syndactyly patients, sample naming convention S1..S7
sample_ids=("1" "2" "3" "4" "5" "6" "7")

#path to references and databases
truseq3="${REF_PATH}/TruSeq3-PE.fa"
humanref="${REF_PATH}/hg38.fa"
dbSNP="${REF_PATH}/common_all_20180418.vcf"        #currently unused (no BQSR step in this pipeline)
hg38exon="${REF_PATH}/ucsc_exon_hg38.bed"
baits="${REF_PATH}/syndactyly_samples_baits.bed"
access="${REF_PATH}/access-10kb.hg38.bed"
hg38excl="${REF_PATH}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.r101.s501.blacklist"

#ANNOVAR (not on bioconda - see conda_envs.txt - called directly via perl)
annovar_path="${TOOLS_PATH}/annovar"
annovar_table="${annovar_path}/table_annovar.pl"
humandb="${REF_PATH}/humandb"
annovar_buildver="hg38"

#protocol/operation sets per variant type (match original per-script choices)
protocol_short="refGeneWithVer,cytoBand,gnomad211_exome,avsnp151,dbnsfp47a,clinvar_20240917"
operation_short="g,r,f,f,f,f"
protocol_sv="refGeneWithVer,cytoBand,gnomad211_exome,avsnp151,dbnsfp47a,dgvMerged"
operation_sv="g,r,f,f,f,r"
protocol_cnv="${protocol_sv}"
operation_cnv="${operation_sv}"
annovar_arg="'-hgvs',,,,,"

#NOTE: cnvkit / guess_baits.py are not documented in conda_envs.txt.
#This script assumes a "cnvkit_env" conda environment exists
#(e.g. `conda create -n cnvkit_env -c bioconda cnvkit`) - create it before running
#the CNV CALLING section if it does not already exist.

#create working directories
PREPROCESSING_SYN="${SAMPLE_PATH}/preprocessing_syndactyly"
MAPPING_SYN="${SAMPLE_PATH}/mapping_syndactyly"
SHORT_CALLING_SYN="${SAMPLE_PATH}/short_variant_calling_syndactyly"
SV_CALLING_SYN="${SAMPLE_PATH}/SV_calling_syndactyly"
CNV_CALLING_SYN="${SAMPLE_PATH}/CNV_calling_syndactyly"
ANNOTATION_SYN="${SAMPLE_PATH}/annotation_syndactyly"

mkdir -p "${SAMPLE_PATH}/logs"
mkdir -p "${PREPROCESSING_SYN}/fastqc_raw"
mkdir -p "${PREPROCESSING_SYN}/fastqc_pp"
mkdir -p "${SHORT_CALLING_SYN}/merged"
mkdir -p "${SV_CALLING_SYN}"
mkdir -p "${CNV_CALLING_SYN}"
mkdir -p "${ANNOTATION_SYN}/short_variants"
mkdir -p "${ANNOTATION_SYN}/SVs"
mkdir -p "${ANNOTATION_SYN}/CNVs"

for sample_id in "${sample_ids[@]}"; do
	mkdir -p "${MAPPING_SYN}/S${sample_id}"
	mkdir -p "${SHORT_CALLING_SYN}/S${sample_id}"
done

# Dual-log setup: main run log and dedicated failure/skip tracking log
LOG="${SAMPLE_PATH}/logs/syndactyly_WES_pipeline.log"
FAIL_LOG="${SAMPLE_PATH}/logs/syndactyly_WES_pipeline.failed_skipped.log"
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
echo " Syndactyly 7-patient WES pipeline - Started: $(date)"
echo " Full log           : ${LOG}"
echo " Failure & skip log : ${FAIL_LOG}"
echo "============================================================"

	################
	#PREPROCESSING #
	################

for sample_id in "${sample_ids[@]}"; do

	sample="S${sample_id}"
	MAPPING_SAMPLE="${MAPPING_SYN}/${sample}"

	echo -e "\e[31m ======================= \e[0m"
	echo -e "\e[31m SAMPLE DATA PREPARATION: ${sample} \e[0m"
	echo -e "\e[31m ======================= \e[0m"

	read1="${SAMPLE_PATH}/${sample}.R1.fastq.gz"
	read2="${SAMPLE_PATH}/${sample}.R2.fastq.gz"

	echo -e "\e[31m ========================== \e[0m"
	echo -e "\e[31m FASTQC: ${sample} (RAW) \e[0m"
	echo -e "\e[31m ========================== \e[0m"

	conda run -n preprocessing fastqc \
	--threads ${threads} \
	--outdir "${PREPROCESSING_SYN}/fastqc_raw" \
	"${read1}" "${read2}"

	echo -e "\e[31m ================ \e[0m"
	echo -e "\e[31m TRIMMOMATIC: ${sample} \e[0m"
	echo -e "\e[31m ================ \e[0m"

	trim_output_1="${PREPROCESSING_SYN}/${sample}.R1"
	trim_output_2="${PREPROCESSING_SYN}/${sample}.R2"

	#BUG FIX: a trailing "\" after MINLEN:50 in the original script left the
	#command open, so the blank line + "read1=" assignment below it were
	#swallowed into the trimmomatic call. Continuation now ends cleanly.
	conda run -n preprocessing trimmomatic PE \
	-threads ${threads} \
	${read1} ${read2} \
	${trim_output_1}.paired.fastq.gz ${trim_output_1}.unpaired.fastq.gz \
	${trim_output_2}.paired.fastq.gz ${trim_output_2}.unpaired.fastq.gz \
	ILLUMINACLIP:"${truseq3}":2:30:10:8:true HEADCROP:3 TRAILING:10 MINLEN:50

	read1t="${trim_output_1}.paired.fastq.gz"
	read2t="${trim_output_2}.paired.fastq.gz"

	echo -e "\e[31m ============== \e[0m"
	echo -e "\e[31m FASTQC: ${sample} (PP) \e[0m"
	echo -e "\e[31m ============== \e[0m"

	conda run -n preprocessing fastqc \
	--threads ${threads} \
	--outdir "${PREPROCESSING_SYN}/fastqc_pp" \
	"${read1t}" "${read2t}"

	echo -e "\e[31m =============== \e[0m"
	echo -e "\e[31m MAPPING WITH BWA: ${sample} \e[0m"
	echo -e "\e[31m =============== \e[0m"

	sam_file="${MAPPING_SAMPLE}/${sample}.bam"
	conda run -n mapping bwa mem \
	-t ${threads} \
	"${humanref}" "${read1t}" "${read2t}" \
	> "${sam_file}"

	echo -e "\e[31m ================ \e[0m"
	echo -e "\e[31m SORT WITH PICARD: ${sample} \e[0m"
	echo -e "\e[31m ================ \e[0m"

	sort_bam="${MAPPING_SAMPLE}/${sample}.sort.bam"
	conda run -n mapping picard SortSam \
	INPUT="${sam_file}" \
	OUTPUT="${sort_bam}" \
	SORT_ORDER=queryname \
	CREATE_INDEX=true \
	VALIDATION_STRINGENCY=SILENT

	echo -e "\e[31m =================== \e[0m"
	echo -e "\e[31m FIXMATE WITH PICARD: ${sample} \e[0m"
	echo -e "\e[31m =================== \e[0m"

	fixmate_bam="${MAPPING_SAMPLE}/${sample}.sort.fixmate.bam"
	conda run -n mapping picard FixMateInformation \
	INPUT="${sort_bam}" \
	OUTPUT="${fixmate_bam}" \
	SORT_ORDER=coordinate \
	CREATE_INDEX=true \
	VALIDATION_STRINGENCY=SILENT

	echo -e "\e[31m ===================== \e[0m"
	echo -e "\e[31m ADD GROUP WITH PICARD: ${sample} \e[0m"
	echo -e "\e[31m ===================== \e[0m"

	group_bam="${MAPPING_SAMPLE}/${sample}.sort.fixmate.group.bam"
	conda run -n mapping picard AddOrReplaceReadGroups \
	INPUT="${fixmate_bam}" \
	OUTPUT="${group_bam}" \
	SORT_ORDER=coordinate \
	RGID="${sample_id}" \
	RGSM="${sample}" \
	RGLB=lib1 \
	RGPL=illumina \
	RGPU=unit1 \
	CREATE_INDEX=true

	echo -e "\e[31m ============== \e[0m"
	echo -e "\e[31m FILTER BY GATK: ${sample} \e[0m"
	echo -e "\e[31m ============== \e[0m"

	#BUG FIX: trailing "\" after --read-filter PairedReadFilter left the
	#command open in the original script
	filter_bam="${MAPPING_SAMPLE}/${sample}.sort.fixmate.group.filter.bam"
	conda run -n Hcalling gatk PrintReads \
	-I "${group_bam}" \
	-O "${filter_bam}" \
	--read-filter MappedReadFilter \
	--read-filter PairedReadFilter

	echo -e "\e[31m ============================= \e[0m"
	echo -e "\e[31m REMOVE DUPLICATES WITH PICARD: ${sample} \e[0m"
	echo -e "\e[31m ============================= \e[0m"

	#BUG FIX: missing "\" between REMOVE_DUPLICATES=true and CREATE_INDEX=true
	#in the original script meant CREATE_INDEX=true was never passed to
	#MarkDuplicates (it ran as an unrelated shell assignment on its own line)
	rmdup_bam="${MAPPING_SAMPLE}/${sample}.sort.fixmate.group.filter.rmdup.bam"
	rmdup_metrics="${MAPPING_SAMPLE}/${sample}.sort.rmdp.txt"
	conda run -n mapping picard MarkDuplicates \
	INPUT="${filter_bam}" \
	OUTPUT="${rmdup_bam}" \
	METRICS_FILE="${rmdup_metrics}" \
	REMOVE_DUPLICATES=true \
	CREATE_INDEX=true

	echo -e "\e[31m ================================== \e[0m"
	echo -e "\e[31m FILTER LOW QUALITY MAPPING BY GATK: ${sample} \e[0m"
	echo -e "\e[31m ================================== \e[0m"

	#short-variant preprocessing (MQ >= 10)
	flqm_short_bam="${MAPPING_SAMPLE}/${sample}.sort.fixmate.group.filter.rmdup.FLQM.bam"
	conda run -n Hcalling gatk PrintReads \
	-I "${rmdup_bam}" \
	-O "${flqm_short_bam}" \
	--read-filter MappingQualityReadFilter \
	--minimum-mapping-quality 10

	#structural-variant preprocessing (MQ >= 1, more permissive for split/discordant reads)
	flqm_sv_bam="${MAPPING_SAMPLE}/${sample}_SVs.sort.fixmate.group.filter.rmdup.FLQM.bam"
	conda run -n Hcalling gatk PrintReads \
	-I "${rmdup_bam}" \
	-O "${flqm_sv_bam}" \
	--read-filter MappingQualityReadFilter \
	--minimum-mapping-quality 1

	###########################
	#REMOVE INTERMEDIATE FILES#
	###########################

	if [ "${cleanup_intermediates}" = "true" ]; then
		echo -e "\e[32m Removing intermediate files for ${sample} (disk space)... \e[0m"
		rm -f "${trim_output_1}.unpaired.fastq.gz" "${trim_output_2}.unpaired.fastq.gz"
		rm -f "${sam_file}"
		rm -f "${sort_bam}" "${MAPPING_SAMPLE}/${sample}.sort.bai"
		rm -f "${fixmate_bam}" "${MAPPING_SAMPLE}/${sample}.sort.fixmate.bai"
		rm -f "${group_bam}" "${MAPPING_SAMPLE}/${sample}.sort.fixmate.group.bai"
		rm -f "${filter_bam}" "${MAPPING_SAMPLE}/${sample}.sort.fixmate.group.filter.bai"
	else
		echo -e "\e[32m cleanup_intermediates=false - keeping intermediate BAMs for ${sample} \e[0m"
	fi

	echo -e "\e[32m Preprocessing complete: ${sample} \e[0m"

done

	echo -e "\e[32m PREPROCESSING COMPLETED FOR ALL SAMPLES \e[0m"

	#Build path arrays reused by the SV/CNV sections below
	flqm_short_bams=()
	flqm_sv_bams=()
for sample_id in "${sample_ids[@]}"; do
	sample="S${sample_id}"
	flqm_short_bams+=("${MAPPING_SYN}/${sample}/${sample}.sort.fixmate.group.filter.rmdup.FLQM.bam")
	flqm_sv_bams+=("${MAPPING_SYN}/${sample}/${sample}_SVs.sort.fixmate.group.filter.rmdup.FLQM.bam")
done

	##########################
	#SHORT VARIANT CALLING   #
	##########################

for sample_id in "${sample_ids[@]}"; do

	sample="S${sample_id}"
	CALLING_SAMPLE="${SHORT_CALLING_SYN}/${sample}"
	flqm_bam="${MAPPING_SYN}/${sample}/${sample}.sort.fixmate.group.filter.rmdup.FLQM.bam"

	echo -e "\e[31m ========================================== \e[0m"
	echo -e "\e[31m CALLING VARIANT BY HAPLOTYPECALLER IN GATK: ${sample} \e[0m"
	echo -e "\e[31m ========================================== \e[0m"

	raw_vcf="${CALLING_SAMPLE}/${sample}.raw.vcf"
	conda run -n Hcalling gatk HaplotypeCaller \
	--native-pair-hmm-threads ${threads} \
	-R "${humanref}" \
	-I "${flqm_bam}" \
	-O "${raw_vcf}"

	echo -e "\e[31m ==================================================== \e[0m"
	echo -e "\e[31m SEPARATE RAW DATA INTO SNPS AND INDELS FILES BY GATK: ${sample} \e[0m"
	echo -e "\e[31m ==================================================== \e[0m"

	raw_snps_vcf="${CALLING_SAMPLE}/${sample}.raw.SNPs.vcf"
	raw_indels_vcf="${CALLING_SAMPLE}/${sample}.raw.indels.vcf"

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

	conda run -n Hcalling gatk IndexFeatureFile --input "${raw_snps_vcf}"
	conda run -n Hcalling gatk IndexFeatureFile --input "${raw_indels_vcf}"

	echo -e "\e[31m ========================== \e[0m"
	echo -e "\e[31m VARIANT FILTRATION BY GATK: ${sample} \e[0m"
	echo -e "\e[31m ========================== \e[0m"

	snps_vcf="${CALLING_SAMPLE}/${sample}.SNPs.vcf"
	indels_vcf="${CALLING_SAMPLE}/${sample}.indels.vcf"

	#BUG FIX: a trailing "\" after the last filter of the SNPs block in the
	#original script merged the indels VariantFiltration call into this one
	#as extra (invalid) arguments - each block now ends cleanly
	conda run -n Hcalling gatk VariantFiltration \
	-R "${humanref}" \
	-V "${raw_snps_vcf}" \
	-O "${snps_vcf}" \
	--filter-expression "QD < 2.0" --filter-name "QualByDepth" \
	--filter-expression "FS > 60.0" --filter-name "FisherStrand" \
	--filter-expression "MQ < 40.0" --filter-name "RMSMappingQuality" \
	--filter-expression "HaplotypeScore > 13.0" --filter-name "HaplotypeScore" \
	--filter-expression "MappingQualityRankSum < -12.5" --filter-name "MappingQualityRankSum" \
	--filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSumTest" \
	--filter-expression "QUAL < 30.0 || DP < 6 || DP > 5000 || HRun > 5" --filter-name "StandardFilters" \
	--filter-expression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" --filter-name "HARD_TO_VALIDATE"

	conda run -n Hcalling gatk VariantFiltration \
	-R "${humanref}" \
	-V "${raw_indels_vcf}" \
	-O "${indels_vcf}" \
	--filter-expression "QD < 2.0" --filter-name "QualByDepth" \
	--filter-expression "FS > 200.0" --filter-name "FisherStrand" \
	--filter-expression "MQ < 40.0" --filter-name "RMSMappingQuality" \
	--filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSumTest" \
	--filter-expression "MQ0 >= 4 && ((MQ0 / (1.0*DP)) > 0.1)" --filter-name "HARD_TO_VALIDATE" \
	--filter-expression "QUAL < 30.0 || DP < 6 || DP > 5000 || HRun > 5" --filter-name "QualFilter"

	echo -e "\e[32m Short-variant calling complete: ${sample} \e[0m"

done

	####################
	#MERGE CALLED FILES#
	####################

	echo -e "\e[31m ====================================== \e[0m"
	echo -e "\e[31m MERGE CALLED FILES WITH BCFTOOLS MERGE \e[0m"
	echo -e "\e[31m ====================================== \e[0m"

	MERGED_SHORT="${SHORT_CALLING_SYN}/merged"
	snps_gz_list=()
	indels_gz_list=()

for sample_id in "${sample_ids[@]}"; do
	sample="S${sample_id}"
	CALLING_SAMPLE="${SHORT_CALLING_SYN}/${sample}"

	snps_gz="${CALLING_SAMPLE}/${sample}.SNPs.vcf.gz"
	indels_gz="${CALLING_SAMPLE}/${sample}.indels.vcf.gz"

	conda run -n Hcalling bgzip -c "${CALLING_SAMPLE}/${sample}.SNPs.vcf" > "${snps_gz}"
	conda run -n Hcalling bgzip -c "${CALLING_SAMPLE}/${sample}.indels.vcf" > "${indels_gz}"

	conda run -n Hcalling tabix -p vcf "${snps_gz}"
	conda run -n Hcalling tabix -p vcf "${indels_gz}"

	snps_gz_list+=("${snps_gz}")
	indels_gz_list+=("${indels_gz}")
done

	merged_snps="${MERGED_SHORT}/merged_samples.SNPs.vcf.gz"
	merged_indels="${MERGED_SHORT}/merged_samples.indels.vcf.gz"

	conda run -n Hcalling bcftools merge \
	"${snps_gz_list[@]}" \
	-o "${merged_snps}" -O z

	conda run -n Hcalling bcftools merge \
	"${indels_gz_list[@]}" \
	-o "${merged_indels}" -O z

	echo -e "\e[31m ================== \e[0m"
	echo -e "\e[31m VARIANT ANNOTATION \e[0m"
	echo -e "\e[31m ================== \e[0m"

	ANNOTATION_SHORT="${ANNOTATION_SYN}/short_variants"

	perl "${annovar_table}" \
	"${merged_snps}" \
	"${humandb}/" \
	-buildver "${annovar_buildver}" \
	-out "${ANNOTATION_SHORT}/merged_samples.SNPs" \
	-remove \
	-protocol "${protocol_short}" \
	-operation "${operation_short}" \
	-arg ${annovar_arg} \
	-nastring . \
	-vcfinput \
	-polish

	perl "${annovar_table}" \
	"${merged_indels}" \
	"${humandb}/" \
	-buildver "${annovar_buildver}" \
	-out "${ANNOTATION_SHORT}/merged_samples.indels" \
	-remove \
	-protocol "${protocol_short}" \
	-operation "${operation_short}" \
	-arg ${annovar_arg} \
	-nastring . \
	-vcfinput \
	-polish

	echo -e "\e[31m =========================== \e[0m"
	echo -e "\e[31m EXTRACT WITH BCFTOOLS QUERY \e[0m"
	echo -e "\e[31m =========================== \e[0m"

	#BUG FIX: "%CADD_phred\%tCLNALLELEID" in the original had the tab
	#character mistyped as "\%t" - fixed to "\t%CLNALLELEID" below. The
	#indels string also had a duplicated "%ALLELE_END" field - removed.
	conda run -n Hcalling bcftools query \
	--print-header \
	-f "%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%ANNOVAR_DATE\t%ALLELE_END\t%Gene.refGeneWithVer\t%Func.refGeneWithVer\t%AAChange.refGeneWithVer\t%cytoBand\t%AF\t%AF_popmax\t%AF_male\t%AF_female\t%AF_eas\t%AF_sas\t%avsnp151\t%SIFT_score\t%Polyphen2_HDIV_score\t%MutationTaster_score\t%PROVEAN_score\t%CADD_phred\t%CLNALLELEID\t%CLNDN[\t%GT:%DP:%GQ]\n" \
	"${ANNOTATION_SHORT}/merged_samples.SNPs.hg38_multianno.vcf" \
	> "${ANNOTATION_SHORT}/merged_samples.SNPs.txt"

	conda run -n Hcalling bcftools query \
	--print-header \
	-f "%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%ANNOVAR_DATE\t%ALLELE_END\t%Gene.refGeneWithVer\t%Func.refGeneWithVer\t%AAChange.refGeneWithVer\t%cytoBand\t%AF\t%AF_popmax\t%AF_male\t%AF_female\t%AF_eas\t%AF_sas\t%avsnp151\t%CADD_phred\t%CLNALLELEID\t%CLNDN\t%CLNSIG[\t%GT:%DP:%GQ]\n" \
	"${ANNOTATION_SHORT}/merged_samples.indels.hg38_multianno.vcf" \
	> "${ANNOTATION_SHORT}/merged_samples.indels.txt"

	echo -e "\e[32m Short-variant calling and annotation complete \e[0m"

	###############################
	#STRUCTURAL VARIANT CALLING   #
	###############################

	echo -e "\e[31m ======================================== \e[0m"
	echo -e "\e[31m CALLING STRUCTURAL VARIANT BY DELLY CALL \e[0m"
	echo -e "\e[31m ======================================== \e[0m"

	sv_bcfs=()
for sample_id in "${sample_ids[@]}"; do
	sample="S${sample_id}"
	#BUG FIX: the original SV script called delly on the short-variant FLQM
	#BAM (MQ>=10), not on the dedicated SV FLQM BAM (MQ>=1) that the
	#preprocessing script generated specifically for SV calling
	sv_bam="${MAPPING_SYN}/${sample}/${sample}_SVs.sort.fixmate.group.filter.rmdup.FLQM.bam"
	sv_bcf="${SV_CALLING_SYN}/${sample}.SVs.delly.bcf"

	conda run -n Hcalling delly call \
	-g "${humanref}" \
	-o "${sv_bcf}" \
	-x "${hg38excl}" \
	"${sv_bam}"

	sv_bcfs+=("${sv_bcf}")
done

	echo -e "\e[31m =============================================== \e[0m"
	echo -e "\e[31m MERGE SV SITES INTO UNIFIED LIST BY DELLY MERGE \e[0m"
	echo -e "\e[31m =============================================== \e[0m"

	merged_sites="${SV_CALLING_SYN}/merged_sites.SVs.delly.bcf"
	conda run -n Hcalling delly merge \
	-o "${merged_sites}" \
	"${sv_bcfs[@]}"

	echo -e "\e[31m ============================================= \e[0m"
	echo -e "\e[31m GENOTYPING MERGED SV SITES LIST BY DELLY CALL \e[0m"
	echo -e "\e[31m ============================================= \e[0m"

	sv_geno_bcfs=()
for sample_id in "${sample_ids[@]}"; do
	sample="S${sample_id}"
	sv_bam="${MAPPING_SYN}/${sample}/${sample}_SVs.sort.fixmate.group.filter.rmdup.FLQM.bam"
	sv_geno_bcf="${SV_CALLING_SYN}/${sample}.SVs.delly.geno.bcf"

	conda run -n Hcalling delly call \
	-g "${humanref}" \
	-v "${merged_sites}" \
	-o "${sv_geno_bcf}" \
	-x "${hg38excl}" \
	"${sv_bam}"

	#BUG FIX: bcftools merge requires each input to be indexed - the
	#original script merged the genotyped BCFs without indexing them first
	conda run -n Hcalling bcftools index "${sv_geno_bcf}"

	sv_geno_bcfs+=("${sv_geno_bcf}")
done

	echo -e "\e[31m =========================================== \e[0m"
	echo -e "\e[31m MERGE GENOTYPED SAMPLES WITH BCFTOOLS MERGE \e[0m"
	echo -e "\e[31m =========================================== \e[0m"

	merged_sv="${SV_CALLING_SYN}/merged_samples.SVs.delly.bcf"
	conda run -n Hcalling bcftools merge \
	-m id -O b -o "${merged_sv}" \
	"${sv_geno_bcfs[@]}"

	conda run -n Hcalling bcftools index "${merged_sv}"

	echo -e "\e[31m ============================== \e[0m"
	echo -e "\e[31m FILTER SAMPLES BY DELLY FILTER \e[0m"
	echo -e "\e[31m ============================== \e[0m"

	sv_filtered="${SV_CALLING_SYN}/SVs.delly.filtered.bcf"
	conda run -n Hcalling delly filter \
	-f germline \
	-o "${sv_filtered}" \
	"${merged_sv}"

	echo -e "\e[31m ============================== \e[0m"
	echo -e "\e[31m PREPARE INPUT FILE FOR ANNOVAR \e[0m"
	echo -e "\e[31m ============================== \e[0m"

	#BUG FIX: a trailing "\" after the bcftools convert call in the original
	#script merged the following table_annovar.pl call into this one
	sv_filtered_vcf="${SV_CALLING_SYN}/SVs.delly.filtered.vcf.gz"
	conda run -n Hcalling bcftools convert \
	--output-type z --output "${sv_filtered_vcf}" \
	"${sv_filtered}"

	ANNOTATION_SV="${ANNOTATION_SYN}/SVs"

	perl "${annovar_table}" \
	"${sv_filtered_vcf}" \
	"${humandb}/" \
	-buildver "${annovar_buildver}" \
	-out "${ANNOTATION_SV}/merged_samples.SVs.delly" \
	-remove \
	-protocol "${protocol_sv}" \
	-operation "${operation_sv}" \
	-arg ${annovar_arg} \
	-nastring . \
	-vcfinput \
	-polish

	echo -e "\e[31m =========================== \e[0m"
	echo -e "\e[31m EXTRACT WITH BCFTOOLS QUERY \e[0m"
	echo -e "\e[31m =========================== \e[0m"

	#BUG FIX: original queried "merged_samples.delly.SVs.hg38_multianno.vcf"
	#(word order swapped), which does not match the "-out" prefix above and
	#would never have been found
	conda run -n Hcalling bcftools query \
	--print-header \
	-f "%CHROM\t%SVTYPE\t%SVLEN\t%POS\t%END\t%CHR2\t%POS2\t%REF\t%ALT\t%FILTER\t%PE\t%SR\t%MAPQ\t%Gene.refGeneWithVer\t%Func.refGeneWithVer\t%AAChange.refGeneWithVer\t%cytoBand\t%CIPOS\t%CIEND\t%AF\t%AF_popmax\t%AF_male\t%AF_female\t%AF_eas\t%AF_sas\t%avsnp151\t%CADD_phred\t%dgvMerged[\t%GT:%GQ:%FT]\n" \
	"${ANNOTATION_SV}/merged_samples.SVs.delly.hg38_multianno.vcf" \
	> "${ANNOTATION_SV}/merged_samples.SVs.delly.txt"

	echo -e "\e[32m Structural-variant calling and annotation complete \e[0m"

	###################
	#CNV CALLING       #
	###################

	echo -e "\e[31m ================ \e[0m"
	echo -e "\e[31m BUILD BAITS FILE \e[0m"
	echo -e "\e[31m ================ \e[0m"

	#BUG FIX: the original guess_baits.py call had trailing "\ " (backslash
	#followed by a space) on several lines, which escapes the space instead
	#of continuing the line and silently truncated the argument list -
	#fixed below. The comment in the original script implied this should
	#only run if no baits .bed file already exists, so it is now guarded.
	if [ ! -f "${baits}" ]; then
		echo -e "\e[32m ${baits} not found - guessing baits from BAM coverage... \e[0m"
		conda run -n Hcalling guess_baits.py \
		"${flqm_short_bams[@]}" \
		-t "${hg38exon}" \
		-o "${baits}"
	else
		echo -e "\e[32m Baits file already present: ${baits} - skipping \e[0m"
	fi

	echo -e "\e[31m ============ \e[0m"
	echo -e "\e[31m CNVKIT BATCH \e[0m"
	echo -e "\e[31m ============ \e[0m"

	CNVKIT_BATCH="${CNV_CALLING_SYN}/cnvkit_batch"
	mkdir -p "${CNVKIT_BATCH}"

	conda run -n Hcalling cnvkit.py batch \
	"${flqm_short_bams[@]}" \
	-n -t "${baits}" -f "${humanref}" --access "${access}" \
	--method hybrid \
	-d "${CNVKIT_BATCH}"

	##################
	#FILTER CNV CALLS#
	##################

	echo -e "\e[31m ================= \e[0m"
	echo -e "\e[31m CNVKIT SEGMETRICS \e[0m"
	echo -e "\e[31m ================= \e[0m"

	cnv_vcfs=()
for sample_id in "${sample_ids[@]}"; do
	sample="S${sample_id}"
	cnr="${CNVKIT_BATCH}/${sample}.sort.fixmate.group.filter.rmdup.FLQM.cnr"
	cns="${CNVKIT_BATCH}/${sample}.sort.fixmate.group.filter.rmdup.FLQM.cns"
	segmetrics_cns="${CNV_CALLING_SYN}/${sample}.segmetrics.cns"

	conda run -n Hcalling cnvkit.py segmetrics \
	"${cnr}" \
	-s "${cns}" \
	--ci --mean \
	-o "${segmetrics_cns}"
done

	echo -e "\e[31m ============ \e[0m"
	echo -e "\e[31m CNVKIT CALLS \e[0m"
	echo -e "\e[31m ============ \e[0m"

for sample_id in "${sample_ids[@]}"; do
	sample="S${sample_id}"
	segmetrics_cns="${CNV_CALLING_SYN}/${sample}.segmetrics.cns"
	filtered_cns="${CNV_CALLING_SYN}/${sample}.filtered.cns"

	#BUG FIX: "--drop-loew-coverage" was a typo for "--drop-low-coverage"
	conda run -n Hcalling cnvkit.py call \
	"${segmetrics_cns}" \
	--filter ci --drop-low-coverage \
	-o "${filtered_cns}"
done

	echo -e "\e[31m =========================== \e[0m"
	echo -e "\e[31m CNVKIT EXPORT TO VCF FORMAT \e[0m"
	echo -e "\e[31m =========================== \e[0m"

for sample_id in "${sample_ids[@]}"; do
	sample="S${sample_id}"
	filtered_cns="${CNV_CALLING_SYN}/${sample}.filtered.cns"
	cnv_vcf="${CNV_CALLING_SYN}/${sample}.CNVs.cnvkit.vcf"
	cnv_gz="${CNV_CALLING_SYN}/${sample}.CNVs.cnvkit.vcf.gz"

	#BUG FIX: original used "S$sample_id}.CNVs.cnvkit.vcf" (missing opening
	#brace) which produced a literal stray "}" in the output filename
	conda run -n Hcalling cnvkit.py export vcf \
	"${filtered_cns}" \
	-i "${sample}" \
	-o "${cnv_vcf}"

	#BUG FIX: bcftools merge requires bgzipped + indexed inputs - the
	#original script merged the plain .vcf files directly
	conda run -n Hcalling bgzip -c "${cnv_vcf}" > "${cnv_gz}"
	conda run -n Hcalling tabix -p vcf "${cnv_gz}"

	cnv_vcfs+=("${cnv_gz}")
done

	#merge and index vcf files
	merged_cnv="${CNV_CALLING_SYN}/merged_samples.CNVs.cnvkit.vcf.gz"
	conda run -n Hcalling bcftools merge \
	-m id -O z -o "${merged_cnv}" \
	"${cnv_vcfs[@]}"

	conda run -n Hcalling bcftools index "${merged_cnv}"

	echo -e "\e[31m ============================== \e[0m"
	echo -e "\e[31m PREPARE INPUT FILE FOR ANNOVAR \e[0m"
	echo -e "\e[31m ============================== \e[0m"

	ANNOTATION_CNV="${ANNOTATION_SYN}/CNVs"

	perl "${annovar_table}" \
	"${merged_cnv}" \
	"${humandb}/" \
	-buildver "${annovar_buildver}" \
	-out "${ANNOTATION_CNV}/merged_samples.CNVs.cnvkit" \
	-remove \
	-protocol "${protocol_cnv}" \
	-operation "${operation_cnv}" \
	-arg ${annovar_arg} \
	-nastring . \
	-vcfinput \
	-polish

	echo -e "\e[31m =========================== \e[0m"
	echo -e "\e[31m EXTRACT WITH BCFTOOLS QUERY \e[0m"
	echo -e "\e[31m =========================== \e[0m"

	conda run -n Hcalling bcftools query \
	--print-header \
	-f "%CHROM\t%POS\t%END\t%REF\t%ALT\t%ID\t%FOLD_CHANGE\t%FOLD_CHANGE_LOG\t%Gene.refGeneWithVer\t%Func.refGeneWithVer\t%AAChange.refGeneWithVer\t%cytoBand\t%CIPOS\t%CIEND\t%AF\t%AF_popmax\t%AF_male\t%AF_female\t%AF_eas\t%AF_sas\t%avsnp151\t%CADD_phred\t%dgvMerged[\t%GT:%GQ:%CN:%CNQ]\n" \
	"${ANNOTATION_CNV}/merged_samples.CNVs.cnvkit.hg38_multianno.vcf" \
	> "${ANNOTATION_CNV}/merged_samples.CNVs.cnvkit.txt"

	echo -e "\e[32m CNV calling and annotation complete \e[0m"

echo ""
echo "========================================================"
echo " PIPELINE COMPLETE - $(date)"
echo "========================================================"
echo "  Preprocessing outputs : ${MAPPING_SYN}/"
echo "  -- SHORT VARIANTS ------------------------------------"
echo "  Per-sample calls      : ${SHORT_CALLING_SYN}/S<id>/"
echo "  Merged SNPs           : ${MERGED_SHORT}/merged_samples.SNPs.vcf.gz"
echo "  Merged indels         : ${MERGED_SHORT}/merged_samples.indels.vcf.gz"
echo "  SNPs annotation       : ${ANNOTATION_SHORT}/merged_samples.SNPs.txt"
echo "  Indels annotation     : ${ANNOTATION_SHORT}/merged_samples.indels.txt"
echo "  -- STRUCTURAL VARIANTS --------------------------------"
echo "  Filtered SVs (bcf)    : ${sv_filtered}"
echo "  SVs annotation        : ${ANNOTATION_SV}/merged_samples.SVs.delly.txt"
echo "  -- CNVs -------------------------------------------------"
echo "  Baits file             : ${baits}"
echo "  Merged CNVs             : ${merged_cnv}"
echo "  CNVs annotation          : ${ANNOTATION_CNV}/merged_samples.CNVs.cnvkit.txt"
echo "  Full execution log       : ${LOG}"

fail_count=$(grep -c '\[FAILED\]\|\[EXECUTION_FAILED\]\|\[OUTPUT_MISSING' "${FAIL_LOG}" 2>/dev/null || echo 0)
skip_count=$(grep -c '\[SKIPPED' "${FAIL_LOG}" 2>/dev/null || echo 0)

if [ "${fail_count}" -gt 0 ]; then
	echo -e "\e[31m  Failure/Issues log       : ${FAIL_LOG} (${fail_count} failures detected!) \e[0m"
	echo -e "\e[31m  >>> Inspect ${FAIL_LOG} to see which steps/samples failed and why. \e[0m"
else
	echo -e "\e[32m  Failure/Issues log       : ${FAIL_LOG} (0 errors recorded) \e[0m"
fi
echo -e "\e[32m  Skipped checkpoints      : ${skip_count} records \e[0m"
echo "========================================================"
