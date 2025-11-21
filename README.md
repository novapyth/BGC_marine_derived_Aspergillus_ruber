# Biosynthetic gene clusters (BGCs) of marine derived fungus (MDF) Aspergillus ruber

## Herein you will find details of reproducible pipeline for Biosynthetic gene clusters (BGCs) mining for Aspergillus ruber genome. It can be very easily applied to other microbial genomes. 

# Biosynthetic Gene Cluster (BGC) Mining Pipeline  
### Marine-Derived Fungus *Aspergillus ruber*

[![Snakemake](https://img.shields.io/badge/Workflow-Snakemake-5C4EE5?logo=snakemake)]()
[![Nextflow](https://img.shields.io/badge/Workflow-Nextflow-3ac382?logo=nextflow)]()
[![Conda](https://img.shields.io/badge/Environment-conda-44A833?logo=anaconda)]()
![License](https://img.shields.io/badge/License-MIT-blue.svg)

This repository contains a fully reproducible pipeline for **Biosynthetic Gene Cluster (BGC) mining** in the marine-derived fungus **_Aspergillus ruber_**.  
It supports **both Snakemake and Nextflow**, enabling flexible, scalable, and reproducible genomics analysis.

---

## 📁 Repository Structure

```text
BGC-Pipeline/
├── BGC_mining_pipeline.sh      # Main Bash script
├── tools.yaml                  # Conda environment for CLI tools
├── python.yaml                 # Conda environment for Python tools
├── config.yaml                 # Common configuration for both workflows
│
├── snakemake/
│   └── Snakefile               # Snakemake pipeline
│
└── nextflow/
    ├── main.nf                 # Nextflow pipeline
    └── nextflow.config         # Nextflow configuration


Summary of Files

File / Directory
Description
BGC_mining_pipeline.sh
Core Bash script orchestrating key BGC mining steps.
tools.yaml
Conda environment (BLAST, HMMER, FastTree, COBALT, ncbi-datasets, antiSMASH, Snakemake/Nextflow support).
python.yaml
Environment for Python tools: PANNZER2, Biopython, utility scripts.
config.yaml
Pipeline configuration shared by Snakemake and Nextflow.
Snakefile
Full Snakemake workflow: download → BLAST → antiSMASH → HMMER → PANNZER2 → reannotation → alignment → tree.
nextflow.config
Nextflow pipeline configuration.
main.nf
Nextflow implementation of the BGC mining workflow.
