# BLAST each sequences which contains ARGs. The best hit is kept. 

import os
import glob
import time
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

E_VALUE_THRESH = 1e-20

for fasta_file in glob.glob("*_seq_coloc.fasta"):
    sample = fasta_file.replace("_seq_coloc.fasta", "")
    output_file = f"{sample}_contig_blast"
    xml_file = f"{sample}.xml"
    print(f"BLAST : {sample}")
    
    if os.path.exists(output_file):
        print(f"{output_file} exist, sample ignored\n")
        continue
    if os.stat(fasta_file).st_size == 0:
        print(f"FASTA empty(0 byte), sample ignored : {sample}\n")
        continue
    records = list(SeqIO.parse(fasta_file, "fasta"))
    if not records:
        print(f"FASTA without sequence, sample ignored : {sample}\n")
        continue

    # Reading FASTA
    with open(fasta_file) as f:
        sequence_data = f.read()
    time.sleep(20)
    
    # Lauching BLAST
    result_handle = NCBIWWW.qblast(
        program = "blastn",
        database="nt",
        sequence = sequence_data,
        format_type="XML"
    )
    # Saving XML
    with open(xml_file, "w") as save_file:
        save_file.write(result_handle.read())
        
    # Parsing XML and keeping best hit
    with open(xml_file, "rb") as handle, open(output_file, "w") as out:
        out.write("query\tsubject\tevalue\tidentity\n")
        for record in NCBIXML.parse(handle):   
            if not record.alignments:
                continue
            best_align = record.alignments[0]
            best_hsp = best_align.hsps[0]
            if best_hsp.expect < E_VALUE_THRESH:
                identity = best_hsp.identities / best_hsp.align_length * 100
                out.write(
                        f"{record.query}\t"
                        f"{best_align.title}\t"
                        f"{best_hsp.expect}\t"
                        f"{identity:.2f}\n"     
                        )
