#!/bin/bash

awk '
BEGIN {FS=","; OFS=","}

# Lecture du fichier CSV (;)
NR==FNR {
    gsub(/\r/, "", $0)

    split($0, a, ";")

    sample = a[3]
    biop[sample] = a[2]

    next
}

# Header
FNR==1 {
    print $0, "bioproject"
    next
}

# Merge
{
    sample = $1

    if (sample in biop) {
        print $0, biop[sample]
    }
}
' sample_collection.csv rpkm.csv > rpkm_biop.csv
