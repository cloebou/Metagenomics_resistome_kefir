# Retriving sequences with the first ARG found in each contigs. 

import csv
from Bio import SeqIO
from pathlib import Path

for csv_file in Path(".").glob("*_assembled.fastq_colocalizations.csv"):
    
    sample = csv_file.name.replace("_assembled.fastq_colocalizations.csv", "")
    fasta_in = f"{sample}_contig_coloc_seq.fasta"
    fasta_out = f"{sample}_seq_coloc.fasta"
    
    seq_dict = SeqIO.to_dict(SeqIO.parse(fasta_in, "fasta"))
    extracted_records = []
    
    with open(csv_file, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            contig = row["read"]
            first_pos = row["ARG positions"].split(";")[0]
            start, end = map(int, first_pos.split(":"))
            start_idx = start 
            end_idx = end
            if contig not in seq_dict:
                print(f"contig {contig} non trouvé dans le FASTA")
                continue
            seq = seq_dict[contig].seq[start_idx:end_idx]
            record = seq_dict[contig][start_idx:end_idx]
            record.id = f"{contig}|ARG_{start}_{end}"
            record.description = ""
            extracted_records.append(record)
    SeqIO.write(extracted_records, fasta_out, "fasta")
    print(f" Extracted sequences → {fasta_out}\n")
