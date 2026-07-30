#!/bin/bash
# Joining the 2 tables containing coverage of the contig,
# Mean length of the reads and total number of reads per sample
# Input : 
# sample_nbreads_length.csv -> sample, number of total reads, mean length
# sample_contig_ARGMGE.csv -> sample, contig, ARG, MGE, classif, cov
awk -F',' '
BEGIN {OFS=","}

NR==FNR {
    gsub(/\r/, "", $0)
    stats[$1]=$2","$3   # count, mean_length
    next
}

# Header
FNR==1 {
    gsub(/\r/, "", $0)
    print $0",count,mean_length"
    next
}

# Merging
{
    gsub(/\r/, "", $0)

    sample=$1

    if (sample in stats) {
        print $0, stats[sample]
    }
}
'  sample_nbreads_length.csv sample_contig_ARGMGE.csv > merged.csv
