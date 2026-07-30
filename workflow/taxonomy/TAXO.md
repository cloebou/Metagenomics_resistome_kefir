# Taxonomy Analysis Module

## Overview
This module performs **taxonomic classification and microbial diversity analysis** using a Nextflow workflow. It includes preprocessing, host removal, taxonomic assignment and diversity analysis to characterize microbial communities.

The taxonomy step integrates into the broader pipeline after FASTQ preprocessing and before downstream statistical and visualization analyses.

---

## Workflow Summary
1. Perform quality control on raw reads (FastQC)
2. Trim reads (Trimmomatic)
3. Remove host contamination (Bowtie2)
4. Classify reads with Kraken2
5. Estimate abundances with Bracken
6. Compute alpha and beta diversity
7. Generate visualizations (family and species level)

---

## Installation
Conda environment :
conda create --name nf-env bioconda::nextflow
conda activate nf-env
export NXF_HOME=$HOME/.nxfhome
mkdir -p $HOME/.nxfhome

---

## Database setup

Host databases :
cd host
wget https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip
unzip GRCh38_noalt_as.zip
wget https://s3.amazonaws.com/igenomes.illumina.com/Bos_taurus/Ensembl/Btau_4.0/Bos_taurus_Ensembl_Btau_4.0.tar.gz
tar -xvzf Bos_taurus_Ensembl_Btau_4.0.tar.gz

Kraken database :
cd kraken
wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_16_GB_20260226.tar.gz
tar -xvzf k2_pluspfp_16_GB_20260226.tar.gz

---

## Run pipeline
NXF_CONDA_ENABLED=true nextflow run main.nf -resume
Remove -resume to start from scratch

---

## Diversity analysis

Alpha diversity :
echo "sample,alpha,value" > alpha.csv
for f in ../results/kraken/*/*.bracken; do
  sample=$(basename "$f" .bracken)
  value=$(python stats/alpha_diversity.py -f "$f" -a Sh)
  echo "${sample},Sh,${value}" >> alpha.csv
done

Beta diversity :
mapfile -d '' bracken_files < <(find ../results/kraken -type f -name "\*.bracken" -print0)
python stats/beta_diversity.py -i "${bracken_files[@]}" --type bracken > braken.txt

---

## Visualization
Plot microbial composition (family and species level)

Rscript plot/visu_families.R
Rscript plot/visu_species.R

---

## Notes
- Adapt the path to FASTQ files in the main.nf (or symbolic links)
- Ensure enough memory for Kraken database (~16GB+)
- Keep consistent sample naming for merging and downstream analysis
