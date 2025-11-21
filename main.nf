#!/usr/bin/env nextflow

params.assembly_accession = params.assembly_accession ?: 'GCA_000600275.1'
params.workdir = params.workdir ?: 'A_ruber_pipeline'
params.threads = params.threads ?: 8

workflow {
    download_genome()
    extract_genome()
    make_blast_db()
    blastp()
    run_antismash()
    hmmer_scan()
    pannzer_run()
    core_reannotation()
    align_and_tree()
}

process download_genome {
    tag "datasets"
    output:
    path "${params.workdir}/genome/ncbi_dataset.zip"
    script:
    """
    mkdir -p ${params.workdir}/genome
    cd ${params.workdir}/genome
    datasets download genome accession ${params.assembly_accession} --include genome,protein,gff3,seq-report --filename ncbi_dataset.zip
    """
}

process extract_genome {
    input:
    path zipfile from file("${params.workdir}/genome/ncbi_dataset.zip").ifExists()
    output:
    path "${params.workdir}/genome/A_ruber_genome.fna"
    path "${params.workdir}/genome/A_ruber_proteins.faa"
    path "${params.workdir}/genome/A_ruber.gff3"
    script:
    """
    unzip -o ${zipfile} -d ${params.workdir}/genome/dataset
    cp ${params.workdir}/genome/dataset/ncbi_dataset/data/${params.assembly_accession}/**/*.fna ${params.workdir}/genome/A_ruber_genome.fna || true
    cp ${params.workdir}/genome/dataset/ncbi_dataset/data/${params.assembly_accession}/**/*.faa ${params.workdir}/genome/A_ruber_proteins.faa || true
    cp ${params.workdir}/genome/dataset/ncbi_dataset/data/${params.assembly_accession}/**/*.gff* ${params.workdir}/genome/A_ruber.gff3 || true
    """
}

process make_blast_db {
    input:
    path proteins from file("${params.workdir}/genome/A_ruber_proteins.faa").ifExists()
    output:
    path "${params.workdir}/blast/A_ruber_prot_db*"
    script:
    """
    mkdir -p ${params.workdir}/blast
    makeblastdb -in ${proteins} -dbtype prot -out ${params.workdir}/blast/A_ruber
