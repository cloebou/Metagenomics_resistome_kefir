# Merging for each sample and each contig the ARGs info and the species related (found with Blast)

import csv
import glob
from pathlib import Path

for f in glob.glob("*_contig_blast"):
    sample = f.replace("_contig_blast","")
    coloc_file = f"{sample}_assembled.fastq_colocalizations.csv"
    blast_file = f"{sample}_contig_blast"
    output_file = f"{sample}_merged.csv"

# Reading blast file -> dictionnary
    blast_info = {}

    with open(blast_file) as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)  # query subject evalue identity

        for row in reader:
            query, subject, evalue, identity = row

            # Extraction of the contig before"|ARG"
            contig = query.split("|ARG")[0]

            blast_info[contig] = {
                "blast_subject": subject,
                "blast_evalue": evalue,
                "blast_identity": identity
            }

# Reading CSV with colocalization and merging
    with open(coloc_file) as coloc, open(output_file, "w", newline="") as out:
    
        coloc_reader = csv.DictReader(coloc)
        fieldnames = (
            coloc_reader.fieldnames
            + ["blast_subject", "blast_evalue", "blast_identity"]
        )
    
        writer = csv.DictWriter(out, fieldnames=fieldnames)
        writer.writeheader()

        for row in coloc_reader:
            contig = row["read"]
    
            if contig in blast_info:
                row.update(blast_info[contig])
            else:
                # If no BLAST for this contig : 
                row.update({
                    "blast_subject": "NA",
                    "blast_evalue": "NA",
                    "blast_identity": "NA"
                })

            writer.writerow(row)

    print(f"Joined file : {output_file}")
