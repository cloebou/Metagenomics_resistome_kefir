# BLAST Processing Module

## Overview

This module performs **sequence similarity searches using BLAST** on contigs generated from the TELComP workflow. It enables the identification and annotation of sequences, followed by merging and aggregation steps across samples.

The BLAST step integrates into the broader pipeline after contig generation and colocalization analysis, and before downstream statistical and visualization analyses.

---

## Workflow Summary
Steps:

1. Extract contigs with at least one ARG 
2. Extract the sequences from contigs with the first ARG detected
3. Perform BLAST search  
4. Merge results per sample  
5. Aggregate all samples into a global dataset  

---

## Package requirement
Python : 
- Bio
- pathlib
- csv
- os
- glob
- time
- Bio.Blast

---

## Notes & Recommendations

- BLAST is not required to be installed since the online version is beeing used. Due to this, the jobs parallelization is limited.  
- Consider filtering BLAST hits (identity, coverage) depending on analysis goals.  
- For the input either : add the path to the TELCoMB output or create symbolic links

## Example Usage

```bash
bash 1_contig.sh
python 2_seq.py
python 3_blast.py
python 4_merged.py
python 5_global_merged.py
