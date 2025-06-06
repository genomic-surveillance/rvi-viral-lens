// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.

include {run_qc_script} from '../modules/run_qc_script.nf'

workflow COMPUTE_QC_METRICS {
    /*
    -----------------------------------------------------------------
    Obtain QC Metrics

    The `COMPUTE_QC_METRICS` workflow is designed to compute quality
    control (QC) metrics for consensus sequences generated from 
    sequencing data. The workflow processes each sample's data to
    evaluate the quality and coverage of the generated consensus 
    sequences. The QC metrics include the percentage of bases 
    covered, the percentage of N bases, the longest segment without
    N bases, and read alignment statistics (total reads aligned,
    unmapped and mapped).

    -----------------------------------------------------------------
    # Inputs
        A channel containing tuples of metadata and the consensus
    FASTA file. Metadata (`meta`) must include the following keys:
        - `id`: Unique identifier combining sample ID and taxonomic
        ID.
        - `taxid`: Taxonomic ID of the sample.
        - `sample_id`: Sample identifier.
        - `bam_file`: Sorted Bam file
        - `ref_files`: Paths to reference genome files.
        - `mpileup_file`: Path to output file from `samtools mpileup`

    -----------------------------------------------------------------
    # Key Processes
        - **QC Metrics Calculation**: Uses a custom QC script (check
        QC script documentation) to compute various quality metrics
        based on the aligned BAM file and the consensus FASTA 
        sequence.

    -----------------------------------------------------------------
    # Outputs
        - `qc_Out_ch`: A channel containing tuples with updated
        metadata that includes computed QC metrics.
    */

    take:
        qc_metrics_In_ch // [meta, fasta_file]

    main:
        qc_metrics_In_ch
            | map {meta, fasta_file ->
                // store fasta_files at meta
                meta.consensus_fa = fasta_file
                tuple(meta, meta.bam_file, fasta_file, meta.ref_files[0], meta.mpileup_file)
            }
            | set{qc_script_In_ch}

        // compute percentage of bases covered
        run_qc_script(qc_script_In_ch)
        
        // populate meta with the qc values
        run_qc_script.out
            | map {meta, qc_csv, stdout_str -> 
                def tokens_lst = stdout_str.tokenize(",")
                meta.percentage_of_N_bases = tokens_lst[1] // Proportion of ambiguous bases in the consensus sequence.
                meta.percentage_genome_coverage = tokens_lst[2] // Proportion of the genome covered by aligned reads.
                meta.longest_no_N_segment = tokens_lst[3] // The longest continuous segment without N bases.
                meta.total_aligned_reads = tokens_lst[4] // Number of reads aligned to the reference.
                meta.total_unmapped_reads = tokens_lst[10].replace("\n","") // Number of reads that did not align.
                meta.total_mapped_reads = tokens_lst[8] // Number of reads that mapped successfully.

                tuple(meta, meta.bam_file)
            }
            | set {qc_Out_ch}
    emit:
        qc_Out_ch // (meta, bam_file)

//-------------------------------------------------------------------
// Note: The QC script's output is expected to be in a specific 
// format; any changes in this format may require adjustments in the
// parsing logic within the workflow.
}


workflow {
    manifest_channel = Channel.fromPath(params.manifest_file)
    | splitCsv(header: true, sep: ',')
    | map { row ->
        def meta = [id:row.id,
            taxid:row.taxid,
            sample_id:row.sample_id,
            bam_file:row.bam_file,
            ref_files:[row.ref_files],
            mpileup_file:row.mpileup_file
        ]
        [meta, row.fasta_file]
    }
    COMPUTE_QC_METRICS(manifest_channel)
}

//-------------------------------------------------------------------
// TODO: We should consider to make reference and bam files
//           explicitly on the input channels.
