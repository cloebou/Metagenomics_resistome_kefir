#!/bin/bash
# Retriving contigs with at least 1 ARG.

set -euo pipefail
shopt -s nullglob

for f in *_contig_coloc.fasta; do
    sample=${f%_contig_coloc.fasta}

    echo "$f"

    fastq_in="${sample}_assembled.fastq"
    fasta_out="${sample}_contig_coloc_seq.fasta"

    [[ -f "$fastq_in" ]] || {
        echo " $fastq_in cannot be found, sample $sample ignoré" >&2
        continue
    }

    > "$fasta_out"

    while read -r head; do
        [[ $head == ">"* ]] || continue

        contig="${head#>}"

        seq=$(awk -v c="$contig" '
            $0 == "@"c { getline; print; exit }
        ' "$fastq_in")

        if [[ -n "$seq" ]]; then
            {
                echo ">$contig"
                echo "$seq"
            } >> "$fasta_out"
        fi

        echo "File $fasta_out with the contig : $contig"
    done < "$f"

done
