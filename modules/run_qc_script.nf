// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.

process run_qc_script {
    /*
    *               Compute QC Metrics

    This process performs quality control (QC) on sequence data using
    a custom Python QC script (available at `bin/` dir). It uses a
    set of input files (BAM, FASTA, and reference files) along with
    the results from samtools mpileup to generate a QC report in CSV
    format. The key steps include calculating statistics on the BAM
    file, parsing mpileup output, and running a custom QC script.
    This process generates a CSV file with QC metrics and a PNG image
    visualizing sequencing depth, which are essential for assessing
    the quality and reliability of the sequencing data. The script
    is a modified version of [the `qc.py` of the ARTIC pipeline approach]
    (https://gitlab.internal.sanger.ac.uk/malariagen1/ncov2019-artic-nf/-/blob/main/bin/qc.py?ref_type=heads).

    NOTE: the gitlab for the ARTIC pipeline is currently only availabe
    at SANGER internal gitlab.
    * ---------------------------------------------------------------
    * Input
        - `meta`: Metadata associated with the sample, assumes the following keys  are available: 
            - `id`: internal pipeline ID
            - `sample_id`: sample ID
            - `taxid`: taxonomic ID
        - `bam`: Path to the BAM file containing aligned sequencing reads.
        - `fasta`: Path to the consensus FASTA file for the sample.
        - `ref`: Path to the reference file against which the reads are aligned.
        - `samtools_mpileup`: Path to the output file from samtools mpileup that contains coverage information.

    * Output
        - A CSV file containing QC metrics (`${meta.id}.qc.csv`).
        - Standard output (`stdout`) prints the first row of the QC CSV file for quick verification.

    * ---------------------------------------------------------------
    */
    tag {meta.id}

    label "qc"
    label "mem_250M"
    label "cpu_1"

    input:
    tuple val(meta), path(bam), path(bam_index), path(fasta), path(ivar_variants)

    output:
    tuple val(meta), path(bam), path( bam_index), path(fasta), path(ivar_variants), path("${meta.id}.qc.json")

    script:
    samtools_flagstat="${meta.id}.flagstat.txt"
    samtools_depth="${meta.id}.depths.txt"
    samtools_bam_header="${meta.id}.sam_header.txt"
    """
    # Generate required samtools flagstat file
    samtools flagstat ${bam} > ${samtools_flagstat}
    samtools depth -a ${bam} > ${samtools_depth}
    samtools view -H ${bam} > ${samtools_bam_header}

    # Run QC script
    qc.py \
        --outfile ${meta.id}.qc.json \
        --fasta_file ${fasta} \
        --samtools_depth_file ${samtools_depth} \
        --samtools_flagstat_file ${samtools_flagstat} \
        --samtools_bam_header_file ${samtools_bam_header} \
        --ivar_variants_file ${ivar_variants}
    """

}

