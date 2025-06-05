// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.

include {write_classification_report} from '../modules/write_classification_report.nf'
include {write_classification_report as write_invalid_classification_report} from '../modules/write_classification_report.nf'
workflow GENERATE_CLASSIFICATION_REPORT {
    /*
    -----------------------------------------------------------------
    Write Classification Report

    The GENERATE_CLASSIFICATION_REPORT workflow generates a 
    classification report based on metadata associated with 
    sequencing samples. This workflow collects metadata from each
    sample, formats the data into a report line, and aggregates these
    lines into a final classification report file.

    -----------------------------------------------------------------
    # Inputs

    - **Metadata Channel**: A channel containing metadata for each
    sample. The metadata must include the following keys:
        - `sample_id`: Unique identifier for the sample.
        - `taxid`: Taxonomic ID of the sample.
        - `ref_selected`: Reference sequence used for the 
        sample, formatted with commas but replaced by pipes (`|`).
        - `virus_name`: Name of the virus detected.
        - `virus_subtype`: Subtype of the virus, if applicable
        (defaults to `'None'` if `null`).
        - `flu_segment`: Influenza segment information, if 
        applicable (defaults to `'None'` if `null`).
        - `percentage_genome_coverage`: Percentage of the genome
        covered by aligned reads.
        - `total_mapped_reads`: Number of reads that mapped 
        successfully.
        - `longest_no_N_segment`: Length of the longest segment
        without ambiguous bases ('N').
        - `percentage_of_N_bases`: Proportion of ambiguous bases
        in the consensus sequence.

    -----------------------------------------------------------------
    # Key Processes:
        - **Report Line Generation**: Each sample's metadata is
        processed to generate a line of text summarizing the
        classification data.
        - **Report Aggregation**: All report lines are aggregated
        into a single report file.

    -----------------------------------------------------------------
    # Outputs
        - Classification Report Channel: A text file containing
        aggregated classification data for all samples, formatted as
        a CSV.
    */

    take:
        meta_ch // meta

    main:
    
        // Create a report line for every sample and then aggregate them
        report_lines_ch = meta_ch
        .branch{ it ->
            sample_cov = it[0].percentage_genome_coverage as Float
            min_cov = params.min_coverage_percent as Float

            valid: sample_cov > min_cov
            not_valid: sample_cov <= min_cov
        }

        report_lines_ch.valid.map{it ->
            // convert null values for type and segments to None strings
            if (it[0].virus_subtype == null){
                virus_subtype='None'
            } else {
                virus_subtype = it[0].virus_subtype
            }

            if (it[0].flu_segment==null){
                flu_segment='None'
            } else {
                flu_segment = it[0].flu_segment
            }
            // sample_id, virus, report_name, virus_name, taxid, reference_selected, flu_segment, 
            // virus_subtype, sample_subtype, percentage_genome_coverage, total_mapped_reads,
            // longest_no_N_segment, percentage_of_N_bases
            def sample_info_1 = "${it[0].sample_id},${it[0].virus},${it[0].report_name},${it[0].virus_name}"
            def sample_info_2 = "${it[0].taxid},${it[0].ref_selected.replace(",","|")},${flu_segment},${virus_subtype},${it[0].sample_subtype}"
            def qc_info_v = "${it[0].percentage_genome_coverage},${it[0].total_mapped_reads.replace("^M", "")},${it[0].longest_no_N_segment},${it[0].percentage_of_N_bases}"
            def mut_info_v = "${it[0].total_mutations},${it[0].n_insertions},${it[0].n_deletions},${it[0].n_snps},${it[0].ti_tv_ratio}"
            "${sample_info_1},${sample_info_2},${qc_info_v},${mut_info_v}\n"
        }.collect().set{valid_lines_ch}

        // Write all of the per-sample report lines to a report file
        write_classification_report(valid_lines_ch, "classification_report")

        // write non valid report
        report_lines_ch.not_valid.map{it ->
            // convert null values for type and segments to None strings
            if (it[0].virus_subtype == null){
                def virus_subtype='None'
            } else {
                def virus_subtype = it[0].virus_subtype
            }

            if (it[0].flu_segment==null){
                def flu_segment='None'
            } else {
                def flu_segment = it[0].flu_segment
            }
            // sample_id, virus, report_name, virus_name, taxid, reference_selected, flu_segment, 
            // virus_subtype, sample_subtype, percentage_genome_coverage, total_mapped_reads,
            // longest_no_N_segment, percentage_of_N_bases
            def sample_info_1 = "${it[0].sample_id},${it[0].virus},${it[0].report_name},${it[0].virus_name}"
            def sample_info_2 = "${it[0].taxid},${it[0].ref_selected.replace(",","|")},${flu_segment},${virus_subtype},${it[0].sample_subtype}"
            def qc_info_v = "${it[0].percentage_genome_coverage},${it[0].total_mapped_reads.replace("^M", "")},${it[0].longest_no_N_segment},${it[0].percentage_of_N_bases}"
            def mut_info_v = "${it[0].total_mutations},${it[0].n_insertions},${it[0].n_deletions},${it[0].n_snps},${it[0].ti_tv_ratio}"
            "${sample_info_1},${sample_info_2},${qc_info_v},${mut_info_v}\n"
        }.collect().set{not_valid_lines_ch}

        write_invalid_classification_report(not_valid_lines_ch, "classification_report_invalid")

    emit:
        write_classification_report.out // report file
/*
---------------------------------------------------------------------
# Example Manifest File

The manifest file should be in CSV format and include headers 
matching the required metadata keys. Example:

```csv
sample_id,taxid,ref_selected,virus_name,virus_subtype,flu_segment,percentage_genome_coverage,total_mapped_reads,longest_no_N_segment,percentage_of_N_bases
sample_001,12345,ref1.fasta,Influenza A,H1N1,Segment 1,95.5,10000,500,2.5
sample_002,67890,ref2.fasta,Influenza B,,Segment 2,90.0,9500,450,5.0
```
---------------------------------------------------------------------
*/

}

workflow {
    manifest_channel = Channel.fromPath(params.manifest_file)
    | splitCsv(header: true, sep: ',')
    | map { row ->
        meta = [[id:row.sample_id,
            taxid:row.taxid,
            ref_selected:row.ref_selected,
            virus_name:row.virus_name,
            virus_subtype:row.virus_subtype,
            flu_segment:row.flu_segment,
            percentage_genome_coverage:row.percentage_genome_coverage,
            total_mapped_reads:row.total_mapped_reads,
            longest_no_N_segment:row.longest_no_N_segment,
            percentage_of_N_bases:row.percentage_of_N_bases
        ]]
    }

    GENERATE_CLASSIFICATION_REPORT(manifest_channel)
}


def check_classification_report_params(){
    /*
    -----------------------------------------------------------------
    Checks for necessary parameters and validates paths to ensure 
    they exist. Logs errors if any required parameters are missing.
    -----------------------------------------------------------------

    - **Output**: Number of errors encountered during the checks.

    -----------------------------------------------------------------

    */
    def errors = 0
    // Was the min_coverage_percent provided?
    if (params.min_coverage_percent == null){
        log.error("No min_coverage_percent set")
        errors +=1
    }

    // if provided, check if it is float
    min_cov_value = params.min_coverage_percent as Float

    // if float, it should be between [0.0,100.0]
    if (min_cov_value < 0.0 || min_cov_value > 100.0){
            log.error("min_coverage_percent value set ($min_cov_value) must be >=0.0 and <=100.0")
            errors +=1
    }

    return errors
}
