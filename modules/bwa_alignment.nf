process bwa_alignment_and_post_processing {
    /*
    * Map reads to reference
    */

    publishDir "${params.results_dir}/", overwrite: true, mode: "copy"

    input:
        tuple val(file_id), path(fastq), val(reference_fasta)

    output:
        tuple val(file_id), path(bam_file), path("${bam_file}.bai")

    script:
        bam_file="${file_id}.bam"

        """
        set -e
        set -o pipefail
        bwa mem -p -Y -K 100000000 -t 1 \
            "${reference_fasta}" \
            "${fastq}" | \
            samtools view -bu - | \
            samtools sort -n - | \
            samtools fixmate - - | \
            samtools sort -o ${bam_file}.sorted.bam 
        
        samtools index ${bam_file}.sorted.bam

        mv ${bam_file}.sorted.bam ${bam_file}
        mv ${bam_file}.sorted.bam.bai ${bam_file}.bai

        """
}

// note we may need to provide the index for the ref genome
// better than generate at run time for every single sample

