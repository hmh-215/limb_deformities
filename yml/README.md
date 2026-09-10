# Conda Environments for limb_deformities

This directory contains the conda environment configuration files required to run the pipelines in this repository.

## Environment Specifications

| Environment File | Environment Name | Tools & Key Packages | Purpose |
|---|---|---|---|
| [`preprocessing.yml`](file:///E:/huong/Projects/github/limb_deformities/yml/preprocessing.yml) | `preprocessing` | FastQC (0.12.1), Trimmomatic (0.38), MultiQC (1.33) | Raw read quality control and adapter trimming |
| [`mapping.yml`](file:///E:/huong/Projects/github/limb_deformities/yml/mapping.yml) | `mapping` | BWA (0.7.19), SAMtools (1.23.1), Picard (3.4.0) | Read alignment, sorting, deduplication |
| [`Hcalling.yml`](file:///E:/huong/Projects/github/limb_deformities/yml/Hcalling.yml) | `Hcalling` | GATK4 (4.6.2.0), Delly (0.7.6), BCFtools (1.23.1), HTSlib, CNVkit (0.9.10) | Germline variant calling (SNPs, indels, SVs, CNVs) |

## Creating Environments

Create each conda environment from its respective YAML file:
```bash
conda env create -f yml/preprocessing.yml
conda env create -f yml/mapping.yml
conda env create -f yml/Hcalling.yml
```

## Note on ANNOVAR
ANNOVAR is a standalone Perl script suite and is not distributed via Bioconda.
It should be downloaded directly from [ANNOVAR](https://annovar.openbioinformatics.org/) and extracted to `${TOOLS_PATH}/annovar`.
Ensure Perl 5 is available in your system path.
