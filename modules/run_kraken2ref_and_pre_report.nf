// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.
params.k2r_polling_mode = "max" // [kmeans or max]

process run_k2r_sort_reads {
    /*
    * Parse kraken report and write tax_to_reads_json

    This process runs Kraken2Ref to parse Kraken reports and sort 
    reads by taxonomic classification. It generates JSON files that
    map taxonomic IDs to read IDs and performs sorting based on the 
    reference-selection output JSON if available.

    * --------------------------------------------------------------
    * Input
        - `meta`: Metadata including identifiers like sample ID "id".
        - `kraken_output`: Path to the Kraken output file containing
            classification results.
        - `kraken_report`: Path to the Kraken report file summarizing
            classifications.
    * Output
        - JSON files mapping taxonomic IDs to read IDs
            (`${meta.id}_tax_to_reads.json`)
        - decomposed JSON file (`${meta.id}_decomposed.json`). The 
            reference-selection output JSON is optional and only
            generated if applicable.

    * --------------------------------------------------------------
    */

    tag "${meta.id} - ${task.attempt} - ${task.memory}"
    cache 'lenient'
    label 'kraken2ref'
    label 'mem_k2r_escalate'

    input:
        tuple val(meta), path(kraken_output), path(kraken_report)

    output:
        tuple val(meta), path("${meta.id}_tax_to_reads.json"), path("${meta.id}_decomposed.json"), optional: true, emit: json_files
    script:


    """
    kraken2ref -s ${meta.id} parse_report -i ${kraken_report} -o ./ \
               -t ${params.min_reads_for_taxid} -m ${params.k2r_polling_mode}
    
    # if empty file, no decomposed json file will be generated
    if [ -e "${meta.id}_decomposed.json" ]; then
        kraken2ref -s ${meta.id} sort_reads -k ${kraken_output} -r ./${meta.id}_decomposed.json -m tree -u
    else
        echo "Warning: JSON file does not exist. Skipping downstream commands."
    fi

    """
/*

# Script Breakdown

1. **Parsing Kraken Report**: 
    ```
    kraken2ref -s ${meta.id} parse_report -i ${kraken_report} -o ./
     -t ${params.min_reads_for_taxid}
    ```

    Parses the Kraken report to generate JSON mapping of taxonomic IDs
    to read IDs. The `-t` option sets the minimum number of reads for
    a taxonomic ID to be included.

2. **Sorting Reads**: 
    ```
    kraken2ref -s ${meta.id} sort_reads -k ${kraken_output} 
    -r ./${meta.id}_decomposed.json -m tree -u
    ```

    Sorts reads using the reference-selection output JSON  if available,
    based on a taxonomy tree structure. If the decomposed JSON file is
    missing, sorting is skipped with a warning message.

> **TODO**: we check for decomposed json file, at this point 
there should not be empty fastq files, so I think we should let the
pipeline brake if this is the case.
*/
}



process run_k2r_dump_fastqs_and_pre_report {
    /*
    * --------------------------------------------------------------
                Write Individual Fastq Files and the Pre Report

    The run_k2r_dump_fastqs_and_pre_report process uses Kraken2Ref
    to extract classified reads into FASTQ files and generate a
    preliminary report based on the taxonomic classification data.
    It processes classified reads and produces a detailed report for
    further analysis.

    * --------------------------------------------------------------
    * Input
        - `classified_fqs`: Paths to paired FASTQ files containing
            classified reads.
        - `json_tax_to_readsid`: JSON file mapping taxonomic IDs to
            read IDs.
        - `decomposed_json`: The reference-selection output JSON
            file representing detailed taxonomic hierarchy.
        - `kraken_report`: Kraken report summarizing classification
            data.

    * Output
        - FASTQ files with extracted reads by taxonomic 
            classification.
        - Preliminary TSV report file containing classification 
            summary data.
    
    * Parameters
        - `k2r_dump_fq_mem`: the amount of memory to be requested
            for this process 

    * --------------------------------------------------------------
    */

    tag "${meta.id} - ${task.attempt} - ${task.memory}"
    cache 'lenient'
    label 'kraken2ref'
    memory params.k2r_dump_fq_mem//'mem_k2r_escalate'

    input:
        tuple val(meta), path(classified_fqs), path(json_tax_to_readsid), path(decomposed_json), path(kraken_report)

    output:
        tuple val(meta), path("*_R{1,2}.fq"), optional: true, emit:fq_files // tuple(meta, [id.tax_id.extracted_{1,2}.fq])
        path("${meta.id}_pre_report.tsv"), optional: true, emit: report_file // pre_report.tsv

    shell:
    fq_1 = classified_fqs[0]
    fq_2 = classified_fqs[1]
    
    '''
    if [ "!{meta.splitted}" = "true" ]; then
        part=$(echo !{fq_1}| awk -F'[.]' '{print $(NF-1)}')
        prefix="${part}-!{meta.id}"
    else
        prefix="!{meta.id}"
    fi

    kraken2ref -s ${prefix} dump_fastqs \
            -fq1 !{fq_1} -fq2 !{fq_2} \
            --tax_to_readsid_path !{json_tax_to_readsid} \
            -o ./ --fq_load_mode !{params.k2r_fq_load_mode} \
            -r !{decomposed_json}

    # write pre_report
    k2r_report.py -i !{decomposed_json} -r !{kraken_report} -out_suffix _pre_report.tsv
    '''
/*
# Script Breakdown

This script is designed to handle the processing of classified FASTQ files using Kraken2Ref. It dynamically sets a prefix for the output files based on whether the input FASTQ files are split parts and then runs Kraken2Ref to extract reads into separate FASTQ files for further analysis.

1. **Determine Prefix Based on File Splitting**
Purpose: Sets the prefix for output files based on whether the FASTQ splitted files are splitted or not.
   1. if `meta.splitted` variable is set to "true", indicating that the FASTQ files are split into parts. then, it extracts the part identifier from the filename of `fq_1` and the prefix is set as `${part}-!{meta.id}`.
   2. If not split, the prefix is set as the ID: `!{meta.id}`.

2. **Run Kraken2Ref to Extract FASTQs**
Purpose: Runs Kraken2Ref to extract reads based on taxonomic classification into separate FASTQ files.
Parameters:
    - `-s ${prefix}`: Sets the prefix for the output files, determined in the previous step.
    - `-fq1 !{fq_1} -fq2 !{fq_2}`: Specifies the input paired FASTQ files.
    - `--tax_to_readsid_path !{json_tax_to_readsid}`: Provides the path to the JSON file that maps taxonomic IDs to read IDs, guiding which reads to extract.
    - `-o ./`: Sets the output directory to the current directory (`./`).
    - `--fq_load_mode !{params.k2r_fq_load_mode}`: Specifies the mode for loading FASTQ files, defined by the parameter `params.k2r_fq_load_mode`.
    - `-r !{decomposed_json}`: Uses the reference-selection output JSON file to provide detailed taxonomic structure for sorting reads.

3. **Generate preliminary report**
Purpose: Produce a preliminary report in TSV format summarizing the classification data.
*/
}

process concatenate_fqs_parts {
    /*
    *              Concatenate splitted fastqs

    This process concatenates FASTQ files from multiple parts into
    final combined FASTQ files for each taxonomic classification. 
    This process ensures that all parts corresponding to the same
    taxonomic ID are merged into single files.

    * --------------------------------------------------------------
    * Input
        - `id`: Identifier for the sample or taxonomic
            classification.
        - `fq_parts`: Paths to FASTQ file parts that need to be
            concatenated.

    * Output
        -Concatenated FASTQ files for each taxonomic classification.
    * --------------------------------------------------------------
    */

    cache 'lenient'

    input:
        tuple val(id), path(fq_parts)

    output:
       tuple val(id), path("*_R{1,2}.fq")

    shell:

    '''
    #!/bin/bash

    # Loop through all fastq files in the current directory
    for file in *_R[12].fq; do
        # Extract taxid
        taxid=$(echo "$file"| awk -F'[_]' '{print $(NF-1)}')

        # Define the output filename
        fq_1_out="!{id}_${taxid}_R1.fq"
        fq_2_out="!{id}_${taxid}_R2.fq"

        if [ -f "$fq_1_out" ]; then
            # concatenation confirmation
            echo "Concatenation for !{id}-${taxid} already done"
        else
            # Find and concatenate all matching files into the output file
            echo "Assembing $output_filename"
            cat *-!{id}_${taxid}_R1.fq >> "$fq_1_out"
            cat *-!{id}_${taxid}_R2.fq >> "$fq_2_out"
        fi
    done
    '''
/*
# Script Breakdown

1. **File Loop**: Iterates over FASTQ files with suffix `_R1.fq` and
    `_R2.fq`, assuming paired-end sequencing data.

2. **Output Filename Construction**: Constructs output filenames 
    based on sample or taxonomic ID and checks if concatenation has
    already been completed to avoid redundancy.

3. **Concatenation**: Concatenates parts into the designated output
    files using cat, appending parts that match the specified naming
    convention for each taxonomic classification.
*/
}