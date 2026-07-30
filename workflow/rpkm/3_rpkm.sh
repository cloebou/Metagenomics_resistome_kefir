#!/bin/bash
# Estimating the relative abundance 
# Normalization with RPKM 
# Reads Per Kilobase Million mapped
# This normalization procedure adjusts for 
# both sample sequencing depth and gene length
# Hypothesis for the estimation :  
# Coverage uniformity all along the contig
# Input : merged.tsv, with at least :
# coverage, read length, total reads in sample 

awk -F',' '
BEGIN {OFS=","}

# Header
FNR==1 {
    gsub(/\r/, "", $0)
    print $0",rpkm"
    next
}
# Calculating RPKM 
{
    gsub(/\r/, "", $0)

    cov=$10
    count=$11
    read_len=$12

    # Calcul RPKM
    rpkm = (cov * 1e9) / (read_len * count)

    print $0, rpkm
}
' merged.csv > rpkm.csv
