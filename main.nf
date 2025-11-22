#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.assembly_accession = params.assembly_accession ?: 'GCA_000600275.1'
params.workdir            = params.workdir ?: 'A_ruber_pipeline'
params.threads            = params.threads ?: 8

workflow {
    download_genome()
    extract_genome()
    run_busco()
    run_blobtoolkit()
    make_blast_db()
    blastp()
    run_antismash()
    hmmer_scan()
    pannzer_run()
    core_reannotation()
    align_and_tree()
}

process download_genome {
    tag "download-genome"
    output:
        path "${params.workdir}/genome/ncbi_dataset.zip"
    script:
    """
    mkdir -p ${params.workdir}/genome
    cd ${params.workdir}/genome
    datasets download genome accession ${params.assembly_accession} \
        --include genome,protein,gff3,seq-report \
        --filename ncbi_dataset.zip
    """
}

process extract_genome {
    tag "extract-genome"
    input:
        path zipfile from file("${params.workdir}/genome/ncbi_dataset.zip").ifExists()
    output:
        path "${params.workdir}/genome/A_ruber_genome.fna"
        path "${params.workdir}/genome/A_ruber_proteins.faa"
        path "${params.workdir}/genome/A_ruber.gff3"
    script:
    """
    mkdir -p ${params.workdir}/genome/extracted
    unzip -o ${zipfile} -d ${params.workdir}/genome/extracted

    find ${params.workdir}/genome/extracted -name "*.fna" | head -1 | xargs cp -t ${params.workdir}/genome/ || true
    find ${params.workdir}/genome/extracted -name "*.faa" | head -1 | xargs cp -t ${params.workdir}/genome/ || true
    find ${params.workdir}/genome/extracted -name "*.gff3" | head -1 | xargs cp -t ${params.workdir}/genome/ || true

    mv ${params.workdir}/genome/*.fna  ${params.workdir}/genome/A_ruber_genome.fna
    mv ${params.workdir}/genome/*.faa  ${params.workdir}/genome/A_ruber_proteins.faa
    mv ${params.workdir}/genome/*.gff3 ${params.workdir}/genome/A_ruber.gff3
    """
}

process run_busco {
    tag "busco-quality"
    input:
        path genome from file("${params.workdir}/genome/A_ruber_genome.fna")
    output:
        path "${params.workdir}/qc/busco/"
    script:
    """
    mkdir -p ${params.workdir}/qc/busco
    busco -i ${genome} \
          -l fungi_odb10 \
          -m genome \
          -c ${params.threads} \
          -o A_ruber_busco \
          --out_path ${params.workdir}/qc/busco
    """
}

process run_blobtoolkit {
    tag "blobtoolkit-taxonomy"
    input:
        path genome from file("${params.workdir}/genome/A_ruber_genome.fna")
    output:
        path "${params.workdir}/qc/blobtoolkit/"
    script:
    """
    mkdir -p ${params.workdir}/qc/blobtoolkit
    blobtools create --fasta ${genome} --out ${params.workdir}/qc/blobtoolkit/A_ruber_blobdb

    # Alignment against nt or custom database recommended
    # Example using Diamond NR
    diamond blastx -d /db/nr.dmnd \
        -q ${genome} \
        -o ${params.workdir}/qc/blobtoolkit/A_ruber.diamond.tsv \
        --threads ${params.threads}

    blobtools taxify --hits ${params.workdir}/qc/blobtoolkit/A_ruber.diamond.tsv \
                     --db /db/names.dmp \
                     --nodes /db/nodes.dmp \
                     --out ${params.workdir}/qc/blobtoolkit/taxified_hits

    blobtools view --out ${params.workdir}/qc/blobtoolkit/view \
                   ${params.workdir}/qc/blobtoolkit/A_ruber_blobdb

    blobtools plot --out ${params.workdir}/qc/blobtoolkit/plots \
                   ${params.workdir}/qc/blobtoolkit/A_ruber_blobdb
    """
}
