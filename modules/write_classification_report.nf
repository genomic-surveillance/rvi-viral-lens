// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.
process write_classification_report {
    /*
    * ---------------------------------------------------------------
    * Write classification report

    The process in this Nextflow pipeline generates a classification
    report in CSV format by writing provided lines of data into the
    report file. This report summarizes various metrics and
    classifications related to the sequencing analysis, including
    information on sample identification, taxonomic classification,
    genome coverage, and more.

    * ---------------------------------------------------------------
    * Input
        - A list of lines that contain data to be written into the
        classification report. Each line corresponds to a data entry
        related to the classification results of samples.

    * Output
        - Outputs the path to the generated classification report file
        (`classification_report.csv`).

    * ---------------------------------------------------------------
    */
    label "run_output"
    //publishDir "${params.outdir}/", mode: 'copy'

    input:
        val(list_of_report_lines)
        val(report_filename)
    output:
        path(output_report_file)

    script:
        output_report_file = "${report_filename}.csv"
        // Replace " with ' to prevent issues with writing lines to file
        sample_headers_1="Sample_ID,Virus_Taxon_ID,Virus,Species"
        sample_headers_2="Reference_Taxon_ID,Selected_Reference,Flu_Segment,Reference_Subtype,Sample_Subtype"
        qc_headers="Percentage_of_Genome_Covered,Total_Mapped_Reads,Longest_non_N_segment,Percentage_of_N_bases"
        mut_info_headers = "total_mutations,n_insertions,n_deletions,n_snps,ti_tv_ratio"
        report_lines = list_of_report_lines.join("").replaceAll(/"/, "'")

        """
        # Write header to output report file
        echo "${sample_headers_1},${sample_headers_2},${qc_headers},${mut_info_headers}" > ${output_report_file}_pre

        # Write data lines to report file
        echo "${report_lines}" >> ${output_report_file}_pre
        sed -e "s/\r//g" ${output_report_file}_pre > ${output_report_file}
        # NOTE: the sed expression is there to remove "^M" added characteres
        """
/*
    * ---------------------------------------------------------------
# Script Breakdown

1. **Output File Name**:
The output report file is named classification_report.csv.

2. **Replace Double Quotes**: 
    - `report_lines = list_of_report_lines.join("").replaceAll(/"/, "'")`
    - This line of Groovy script joins all the report lines into a
    single string and replaces double quotes (") with single quotes
    ('). This substitution helps prevent potential formatting issues
    when writing to the CSV file.

3. **Write Header**:

    - Command: 
    ```
    echo "Sample_ID,...,Percentage_of_N_bases" > ${output_report_file}_pre`
    ```
    - Writes the header line to the pre-output file 
    (${output_report_file}_pre). The header defines the columns in the
    CSV file, which include various identifiers and metrics related to
    the samples and their classifications.

4. **Write Data Lines**:

    - Command: `echo "${report_lines}" >> ${output_report_file}_pre`
    - Appends the processed data lines (${report_lines}) to the pre-output file.

5. **Remove Carriage Return Characters**:

    - Command: 
    ```
    sed -e "s/\r//g" ${output_report_file}_pre > ${output_report_file}`
    ```
    - This command uses sed to remove any carriage return (\r) characters
    that may be present in the file. These characters can be introduced 
    when handling files across different operating systems (e.g., 
    Windows vs. Unix-based systems), and their removal ensures 
    consistent file formatting.

*/
}
