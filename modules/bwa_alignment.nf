// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.
process bwa_alignment_and_post_processing {
    /*
    *              Map reads to reference
    
    The process is responsible for mapping sequencing reads
    to a reference genome using BWA (Burrows-Wheeler Aligner),
    followed by post-processing steps including conversion to
    BAM format, sorting, and indexing. At current version of
    this pipeline, this process generates high-quality, sorted
    BAM files that are essential for consensus sequence 
    generation.
    
    * --------------------------------------------------------------
    * Input Tuple:
        - `meta`: Metadata associated with the sample. The process
            assumes there is:
            - `taxid` : taxonomic ID
            - `sample_id`: sample ID
        - `fastq`: Path to the FASTQ file(s) containing sequencing
            reads to be aligned.
        - `ref_files`: Path to the reference files required for 
            alignment, including the reference FASTA 
            (`${meta.taxid}.fa`).
    
    * Output Tuple:
        - Input tuple, plus an additional list comprising:
          - Sorted BAM file (`*.sorted.bam`).
          - Index file for the sorted BAM (`*.sorted.bam.bai`).
    
    * ---------------------------------------------------------------
    * 
    */
    
    input:
        tuple val(meta), path(fastq), path(ref_fa), path(ref_indices)

    output:
        tuple val(meta), path(fastq), path(ref_fa), path(ref_indices), path("${meta.id}.sorted.bam"), path("${meta.id}.sorted.bam.bai")

    script:
        """
        set -e
        set -o pipefail

        #for each fasta file, get simple file name and run bwa mem
        bwa mem ${ref_fa} ${fastq} > ${meta.id}.sam

        # convert sam to bam
        samtools view -S -b ${meta.id}.sam -o ${meta.id}.bam

        # sort alignment by leftmost coordinates
        samtools sort ${meta.id}.bam -o ${meta.id}.sorted.bam

        # generate indexes
        samtools index ${meta.id}.sorted.bam
        """
}

/*

# Script Breakdown

1. **Error Handling**:

- `set -e`: Exits the script if any command returns a non-zero status,
    ensuring that errors are caught early.
- `set -o pipefail`: Ensures that errors in piped commands are 
    correctly propagated.

2. **BWA Alignment**:

- Command: `bwa mem ${ref_fa} ${fastq} > ${meta.id}.sam`
  - Aligns the input sequencing reads (`${fastq}`) to the specified
    reference genome (`${ref_fa}`), producing a SAM file 
    (`${meta.id}.sam`) containing the raw alignments.

3. **Convert SAM to BAM**:

- Command: `samtools view -S -b ${meta.id}.sam -o ${meta.id}.bam`
  - Converts the SAM file to BAM format, which is a compressed binary 
  version of the alignment file, making it more efficient for storage 
  and processing.

4. **Sort BAM File**:

- Command: `samtools sort ${meta.id}.bam -o ${meta.id}.sorted.bam`
  - Sorts the BAM file by the leftmost coordinates of the alignments,
  a necessary step for many downstream analyses.

5. **Generate BAM Index**:

- Command: `samtools index ${meta.id}.sorted.bam`
  - Indexes the sorted BAM file, creating an index file 
  (`*.sorted.bam.bai`) that allows for fast access to specific regions
  of the alignments.
*/