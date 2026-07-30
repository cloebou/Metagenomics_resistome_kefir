#!/bin/bash
# To download SRA files. Names of the target shoud be in samples.txt.

INPUT="samples.txt"
OUTDIR="fastq_samples"
LOGDIR="logs"

mkdir -p "$OUTDIR" "$LOGDIR"

while read -r SAMPLE; do
    echo "Processing $SAMPLE ..."
    fasterq-dump "$SAMPLE" -O "$OUTDIR" &> "$LOGDIR/${SAMPLE}.log"

    if [[ $? -eq 0 ]]; then
        echo "$SAMPLE end"
    else
        echo " Error in $SAMPLE, go see $LOGDIR/${SAMPLE}.log"
    fi
done < "$INPUT"
