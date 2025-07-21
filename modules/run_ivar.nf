// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.
params.ivar_min_depth = 10
params.ivar_freq_threshold = 0.75 // Minimum frequency threshold(0 - 1) to call variants (Ivar Default: 0.03)
params.ivar_min_quality_threshold = 0.20 // Minimum quality score threshold to count base (Ivar Default: 0.20)
params.qc_minimum_depth = 10

process run_ivar {
  /*
  *       Obtain Consensus Sequences using ivar
  
   The run_ivar process in this Nextflow pipeline is designed to 
   reproduce the ARTIC pipeline approach for generating consensus
   sequences from BAM files using the samtools and ivar tools.
   Its main objective is create consensus sequences based on
   specified depth and frequency thresholds. The resulting files 
   are organized by sample and taxonomic ID in the specified
   results directory.
  
  The ARTIC Pipeline can be found at:
  https://gitlab.internal.sanger.ac.uk/malariagen1/ncov2019-artic-nf/-/tree/main?ref_type=heads

  * ----------------------------------------------------------------
  * Input:
   - `meta`: Metadata associated with the sample, the process 
     assumes it contains:
       - pipeline internal ID (`meta.id`)
       - taxonomic ID (`meta.taxid`)
       - sample ID (`meta.sample_id`)
   - `bams`: Path to the BAM file(s) to be processed.
  
  * Output:
   - `meta`: same metadata provided as the input
   - `${meta.id}.consensus.fa`: the path to the generated consensus 
       FASTA file (`.fa`) for the given sample.
   - `mpileup_output`: output file from `samtools mpileup`.
  
  * Parameters:
   - `ivar_min_depth`: The minimum depth required to make a base call
        in the consensus sequence. If the depth at a given position is
        below this threshold, an `'N'` will be used instead. The 
        default value is `10`.
   - `ivar_freq_threshold`: Minimum frequency threshold (0 - 1) to call
        consensus. Bases with a frequency below this threshold are not 
        called. The default value is `0.75`.
   - `ivar_min_quality_threshold`: Minimum quality score threshold to 
        count base (Ivar Default: 0.20) 
  
  * -----------------------------------------------------------------
  */

  tag "${meta.id}"
  label "sample_output"
  label "ivar"

  input:
    tuple val(meta), path(bams), path(reference_fasta)

  output:
    tuple val(meta), path(bams), path("${meta.id}.consensus.fa"), path("${meta.id}_mutations.tsv")

  script:
    sorted_bam = "${meta.id}.sorted.bam"
    mpileup_output="${meta.id}.mpileup.txt"
    flagstat_output="${meta.id}.flagstat.txt"
    depth_output="${meta.id}.depths.txt"
    sam_header_output="${meta.id}.sam_header.txt"
    """
    set -e
    set -o pipefail

    samtools mpileup -aa -A -B -d 0 -Q0 ${sorted_bam} > ${mpileup_output}
    cat ${mpileup_output} | ivar consensus -t ${params.ivar_freq_threshold} -m ${params.ivar_min_depth} -n N -p ${meta.id}.consensus
    cat ${mpileup_output} | ivar variants -t ${params.ivar_freq_threshold} -q ${params.ivar_min_quality_threshold} -r ${reference_fasta} -p ${meta.id}_mutations
    """
}

/*
# Script Breakdown

- `set -e`: This command ensures that the script exits immediately if 
    a command exits with a non-zero status, providing error handling 
    during execution.

- `set -o pipefail`: This option causes the pipeline to return the 
    exit status of the last command in the pipeline that failed, which
    helps in debugging errors.

- `samtools mpileup`: This command generates a pileup of reads from the
    BAM file. The options used are:

    - `-aa`: Output absolutely all positions, including unused 
      reference sequences and zero depth.

    - `-A`: Do not skip anomalous read pairs in variant calling. 
      - Anomalous read pairs are those marked in the FLAG field as 
        paired in sequencing but without the properly-paired flag set.

    - `-B`: Disable BAQ 
      - ([Base Alignment Quality](https://academic.oup.com/bioinformatics/article/27/8/1157/227268)).

    - `-d 0`: Set the maximum depth per file to unlimited. 
      - At a position, read maximally INT reads per input file. Setting
        this limit reduces the amount of memory and time needed to process
        regions with very high coverage. Passing zero for this option sets
        it to the highest possible value, effectively removing the depth 
        limit.

    - `-Q0`: Set the minimum base quality to 0.
      - Minimum base quality for a base to be considered. Note base-quality
        0 is used as a filtering mechanism for overlap removal which marks 
        bases as having quality zero and lets the base quality filter remove
        them. Hence using `-Q0` will make the overlapping bases
        reappear, albeit with quality zero.

    - The output is stored on a text file (`$mpileup_output`)

- `ivar consensus`: This command creates a consensus sequence from the pileup data:
    - `-t ${params.ivar_freq_threshold}`: Sets the threshold for base frequency.
    - `-m ${params.ivar_min_depth}`: Minimum depth to call consensus.
    - `-n N`: Set `N` as the character to print in regions with less than minimum coverage.
    - `-p ${meta.id}.consensus`: Specifies the prefix for the output consensus sequence file.

- `ivar variants`:  call variants - single nucleotide variants(SNVs) and indels.
    - `-t ${params.ivar_freq_threshold}`: Minimum frequency threshold(0 - 1) to call variants (Default: 0.03)
    - `-q ${params.ivar_min_quality_treshold}`: Minimum quality score threshold to count base (Default: 20)
    - `-r ${reference_fasta}`: Reference file used for alignment
    - `-p ${meta.id}`: Prefix for the output tsv variant file

Information about parameters were obtained from 
  - [samtools mpileup documentation](http://www.htslib.org/doc/samtools-mpileup.html)
  - [ivar documentation](https://andersen-lab.github.io/ivar/html/manualpage.html#**autotoc_md19)

*/
