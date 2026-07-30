#!/bin/bash
# Extracting the contig's coverage
# Input : 'contig_arg_classif.csv' : 1 sample with 1 contig 
# and 1 ARG, all type of WHO classification

awk -F',' '
BEGIN {OFS=","}

NR==1 {
    gsub(/\r/, "", $0)
    print $0",coverage"
    next
}

{
    gsub(/\r/, "", $0)

    contig=$2
    match(contig, /cov_([0-9.]+)/, arr)
    cov = arr[1]

    print $0, cov
}
' contig_arg_classif.csv > sample_contig_ARGMGE.csv
