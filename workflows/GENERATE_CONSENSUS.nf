// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.

import groovy.json.JsonSlurper

include {bwa_alignment_and_post_processing} from '../modules/bwa_alignment.nf'
include {run_ivar} from '../modules/run_ivar.nf'

workflow GENERATE_CONSENSUS {
    /*
    -----------------------------------------------------------------
    Obtain Consensus Sequences

    The GENERATE_CONSENSUS workflow performs read alignment and 
    consensus sequence generation for sequencing data. It processes
    paired-end reads by aligning them to reference genomes using BWA,
    followed by consensus calling with iVar. This workflow is designed
    to take in sequencing data for different samples and taxonomic 
    IDs, process them, and produce consensus sequences.

    -----------------------------------------------------------------
    # Inputs
    - **Sample Taxid Channel **: A channel containing tuples of 
    metadata and paired-end FASTQ files. Metadata (`meta`) must
    include the following keys:
        - `id`: Unique identifier combining sample ID and taxonomic
        ID.
        - `taxid`: Taxonomic ID of the sample.
        - `sample_id`: Sample identifier.

    -----------------------------------------------------------------
    # Key Processes
    - **BWA Alignment**: Aligns sequencing reads to the provided
    reference genome.
    - **Consensus Calling with iVar**: Generates consensus sequences
    from the aligned reads.

    -----------------------------------------------------------------
    # Outputs
    - `run_ivar.out`: A channel containing tuples of metadata and the generated consensus FASTA file.

    */

    take:
        sample_taxid_ch // tuple (meta, reads, ref_files)
        
    main:
        // align reads to reference
        bwa_alignment_and_post_processing (sample_taxid_ch)
        bams_ch = bwa_alignment_and_post_processing.out 

        // set ivar input channel
        bams_ch
            | map {meta, fastq, ref_fa, ref_indices, bam, bam_idx ->
                tuple(meta, bam, bam_idx, ref_fa)
            }
            | set {ivar_in_ch} 

        run_ivar(ivar_in_ch)
        
    emit:
        run_ivar.out

}

def check_generate_consensus_params(){
    // < placeholder >
    def errors = 0
    return errors
}
