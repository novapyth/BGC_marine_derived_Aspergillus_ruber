# Biosynthetic gene clusters (BGCs) of marine derived fungus (MDF) Aspergillus ruber

## Herein you will find details of reproducible pipeline for Biosynthetic gene clusters (BGCs) mining for Aspergillus ruber genome. It can be


This document contains the complete files you requested:
Basic Bash Script:

tools.yaml — conda environment for CLI tools (BLAST, HMMER, FastTree, COBALT, ncbi-datasets, Snakemake/Nextflow support, etc.)
python.yaml — conda environment for Python tools (PANNZER2, Biopython, etc.)
config.yaml — pipeline configuration used by Snakemake and Nextflow
Snakefile — Snakemake workflow implementing download → BLAST → antiSMASH → HMMER → PANNZER2 → reannotation → alignment → tree
nextflow.config — Nextflow configuration
main.nf — Nextflow pipeline equivalent to the Snakemake workflow
