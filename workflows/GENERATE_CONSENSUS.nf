// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.

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
        - `ref_files`: Paths to reference genome files.

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

        sample_taxid_ch // tuple (meta, reads)
        
    main:
        // prepare bwa input channel
        sample_taxid_ch
            | map {meta, reads -> tuple(meta,reads,meta.ref_files)}
            | set {bwa_input_ch} // tuple(meta, reads, ref_genome_paths)

        // align reads to reference
        bwa_alignment_and_post_processing (bwa_input_ch)
        bams_ch = bwa_alignment_and_post_processing.out // tuple (meta, [sorted_bam, bai])

        // set ivar input channel
        bams_ch
            | map {meta, bams ->
                // store bam file on meta (check TODO)
                meta.bam_file = bams[0]
                tuple(meta, bams, meta.ref_files[0])}
            | set {ivar_in_ch} // tuple (meta, bam, ref_fasta)

        // generate consensus using ivar
        run_ivar(ivar_in_ch)

        // add mpileup output file to meta
        run_ivar.out // tuple (meta, fasta_file, mpileup_file)
            | map {meta, fasta_file, mpileup_file, stdout -> 
                mut_tokens_lst = stdout.tokenize("---")[-1].tokenize("\n")

                meta.total_mutations = mut_tokens_lst[1].tokenize(":")[-1]
                meta.n_insertions = mut_tokens_lst[2].tokenize(":")[-1]
                meta.n_deletions = mut_tokens_lst[3].tokenize(":")[-1]
                meta.n_snps = mut_tokens_lst[4].tokenize(":")[-1]
                meta.n_ti = mut_tokens_lst[5].tokenize(":")[-1]
                meta.n_tv = mut_tokens_lst[6].tokenize(":")[-1]
                meta.ti_tv_ratio = mut_tokens_lst[7].tokenize(":")[-1]

                meta.mpileup_file = mpileup_file
                tuple(meta, fasta_file)}
            | set {out_ch}

    emit:
        out_ch // tuple (meta, fasta_file)

//-------------------------------------------------------------------
// TODO: We should consider output bam files explicitly instead 
//           of implictily stored on meta.
//-------------------------------------------------------------------
}



def parse_consensus_mnf_meta(consensus_mnf) {
    // consensus_mnf <Channel.fromPath()>
    def mnf_ch =  Channel.fromPath(consensus_mnf)
                        | splitCsv(header: true, sep: ',')
                        | map {row -> 
                            // set meta
                            meta = [sample_id: row.sample_id,
                                    taxid: row.taxid,
                                    ref_files: row.ref_files.split(";").collect()]

                            meta.id = "${row.sample_id}.${row.taxid}"

                            // set files
                            reads = [row.reads_1, row.reads_2]

                            // declare channel shape
                            tuple(meta, reads)
                        }
    return mnf_ch // tuple(index, [fastq_pairs])
}

def check_generate_consensus_params(){
    // < placeholder >
    def errors = 0
    return errors
}
