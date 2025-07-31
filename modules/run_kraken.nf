// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.
process run_kraken {
    /*
    *                 Assign Reads to Taxids
    
    The `run_kraken` process is designed to perform taxonomic 
    classification of sequencing reads using Kraken2, a widely-used
    tool for rapid classification of metagenomic sequences. This
    process classifies reads into taxonomic categories, and outputs
    classified, unclassified reads, and a detailed report, organizing
    the results by sample ID.

    * ---------------------------------------------------------------
    * Input
    - `meta`: Metadata associated with the sample, it assumes it 
        includes
        - `ID`: identifier
    - `fastqs`: Paths to the paired-end FASTQ files that will be
        analyzed.
    - `db_path`: Absolute path to the Kraken2 database used for
        classification.

    * Output
    - Classification results `(*.kraken.output)`.
    - FASTQ files with classified reads `(*.class_seqs*)`.
    - FASTQ files with unclassified reads `(*.unclass_seqs*)`.
    - A summary report `(*.report.txt)`.
    * ---------------------------------------------------------------

    > TODO: We should consider add multithreading options 
        (`--threads ${params.kraken2_threads}`). This add control for
        the user to set the amount of cpus to be allocated.
    * ---------------------------------------------------------------
    */

    tag "${meta.id} - c=${task.cpus} - m=${task.memory}"
    label "kraken"
    // baseline for kraken2 will be 4 CPUs and 2 GB; but both of these are overridable 
    // with command-line options (see config)
    label 'mem_2'
    label "cpu_4"

    input:
        tuple val(meta), path(fastqs) // tuple(sample_id, [fastq_pairs])
        val(db_path) // (absolute) path to kraken DB

    output:
        tuple val(meta), path("*.kraken.output"), path("*.class_seqs*"), path("*.unclass_seqs*"), path("*.report.txt")

    script:
        """
        #!/bin/bash

        set -e
        set -u

        kraken2 \
        --db ${db_path} \
        --output ${meta.id}.kraken.output \
        --paired \
        --classified-out ${meta.id}.class_seqs#.fq \
        --unclassified-out ${meta.id}.unclass_seqs#.fq \
        --report ${meta.id}.report.txt \
        --threads ${task.cpus} \
        ${fastqs}
        """
}

/*
# Command Breakdown

- `set -e`: This command ensures that the script exits immediately if 
    any command exits with a non-zero status, which is useful for
    handling errors during execution.
- `set -u`: This command treats unset variables as an error and exits
    immediately, adding an additional layer of error checking.
- `kraken2`: This command runs Kraken2 to classify the sequencing 
    reads. Options used:
  - `--db ${db_path}`: Specifies the path to the Kraken2 database used
    for classification.
  - `--output ${meta.id}.kraken.output`: Outputs the classification 
    results to a file named with the sample ID and .kraken.output
    extension.
  - `--paired`: Indicates that the input reads are paired-end.
  - `--classified-out ${meta.id}.class_seqs#.fq`: Outputs classified 
    reads to files prefixed with the sample ID and .class_seqs
    extension.
  - `--unclassified-out ${meta.id}.unclass_seqs#.fq`: Outputs 
    unclassified reads to files prefixed with the sample ID and 
    .unclass_seqs extension.
  - `--report ${meta.id}.report.txt`: Generates a detailed report of
    the classification results, including the number of reads 
    classified to each taxonomic level.

For more details, check:
 [Kraken2 classification command documentation](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#classification).

*/
