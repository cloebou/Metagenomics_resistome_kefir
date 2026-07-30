
# RPKM Quantification Module

## Overview
This module computes an estimation of the **RPKM (Reads Per Kilobase per Million reads)** values for ARGs across samples. It integrates coverage and count data to normalize gene abundance and enable comparisons between samples.

The RPKM step is performed after ARG detection and counting, and before downstream statistical analyses.
Hypothesis :  Uniform coverage along contig

---

## Workflow Summary
Steps:

1. Add coverage information to ARG data  
2. Merge all samples  
3. Compute RPKM values  
4. Add bioproject info
