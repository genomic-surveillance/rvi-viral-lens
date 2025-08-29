// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.

include {run_pangolin} from '../modules/run_pangolin.nf'

workflow SCOV2_SUBTYPING {

    /*
    -----------------------------------------------------------------
    Subtype SCOV2 sequences

    The SCOV2_SUBTYPING workflow is designed to determine the 
    SARS-CoV-2 lineage (subtype) of consensus sequences using the
    PANGOLIN tool. This workflow takes in a channel of consensus
    sequences along with their metadata, runs the PANGOLIN lineage
    classification, and outputs updated metadata with the assigned
    lineage.

    -----------------------------------------------------------------
    # Inputs

    - **Consensus Sequence Channel**: A channel containing tuples of
    metadata and consensus sequences. Each tuple should include:
        - `meta`: Metadata for each sample, which must include the
        following keys:
        - `id`: Unique identifier for the sample.
        - `taxid`: Taxonomic ID of the sample.
        - `sample_id`: Sample identifier.
    - `consensus_seq`: The consensus sequence (FASTA file) of the
    sample.

    ```groovy
    // Example of creating the consensus_seq_ch channel
    consensus_seq_ch = Channel.of(
        [ [id: 'sample_001', taxid: '2697049', 
            sample_id: 'sample_001'], 'consensus_sample_001.fasta' ],
        [ [id: 'sample_002', taxid: '2697049',
            sample_id: 'sample_002', 'consensus_sample_002.fasta' ]]
    )
    ```

    -----------------------------------------------------------------
    # Key Processes
    - **SARS-CoV-2 Lineage Classification**: The workflow uses the
    PANGOLIN tool to classify each consensus sequence into a 
    SARS-CoV-2 lineage.

    - **Metadata Update**: The assigned lineage is added to the 
    sample's metadata.

    -----------------------------------------------------------------
    # Outputs
    - **SCOV2 Subtype Channel**: A channel emitting a tuple with updated 
    metadata for each sample, now including the assigned SARS-CoV-2 lineage
    (virus subtype), and the consensus sequence (propagated from input)

    */

    take:
        consensus_seq_ch // tuple (meta, consensus_seq)

    main:
        // get SCOV2 lineage classification
        run_pangolin(consensus_seq_ch)
        run_pangolin.out
            .map { meta, consensus_seq, pangolin_csv -> 
                // Pangolin CSV file will contain one row (single sample), in addition to header
                def rows = pangolin_csv.splitCsv(header:true)
                def new_meta = meta.plus(virus_subtype: rows[0].lineage)
                tuple(new_meta, consensus_seq)
            }
            .set {scov2_subtype_out_ch}

    emit:
        scov2_subtype_out_ch // tuple (meta, consensus_seq)
}
