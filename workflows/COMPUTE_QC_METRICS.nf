// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.

import groovy.json.JsonSlurper

include {run_qc_script} from '../modules/run_qc_script.nf'

workflow COMPUTE_QC_METRICS {
    /*
    -----------------------------------------------------------------
    Compute QC Metrics

    The `COMPUTE_QC_METRICS` workflow is designed to compute quality
    control (QC) metrics for consensus sequences generated from 
    sequencing data. The workflow processes each sample's data to
    evaluate the quality and coverage of the generated consensus 
    sequences. 

    -----------------------------------------------------------------
    # Key Processes
        - **QC Metrics Calculation**: Uses a custom QC script (check
        QC script documentation) to compute various quality metrics
        based on the aligned BAM file and the consensus FASTA 
        sequence.

    -----------------------------------------------------------------
    */

    take:
        qc_metrics_In_ch // [meta, bam, bam_idx, fasta_file, ivar_variants_file]

    main:

        run_qc_script(qc_metrics_In_ch)

        run_qc_script.out
            | map {meta, bam, bam_idx, consensus, variants, qc_json ->
                def json_map = new JsonSlurper().parse(new File(qc_json.toString()))
                // create a new meta with QC metrics
                def new_meta = meta.plus(json_map)
                tuple(new_meta, bam, bam_idx, consensus, variants, qc_json) }
            | set {qc_out_ch}

    emit:
        qc_out_ch // (meta, bams, fasta, variants)
}

workflow {
    manifest_channel = Channel.fromPath(params.manifest_file)
    | splitCsv(header: true, sep: ',')
    | map { row ->
        def meta = [id:row.id,
            taxid:row.taxid,
            sample_id:row.sample_id]
        [meta, row.bam_file, row.fasta_file, row.variants_file, row.qc_json_file]
    }
    COMPUTE_QC_METRICS(manifest_channel)
}
