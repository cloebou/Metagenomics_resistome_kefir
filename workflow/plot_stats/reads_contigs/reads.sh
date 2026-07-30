#!/bin/bash

> 'read.csv'  # écrase le fichier

for file in *_1.fastq
do
    sample=$(basename "$file" _1.fastq)
    count=$(grep -c '^@' "$file")
    echo "${sample},${count}" >> 'read.csv'
done

