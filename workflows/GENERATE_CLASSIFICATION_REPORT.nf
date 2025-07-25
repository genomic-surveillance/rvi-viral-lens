// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.

include {write_classification_report_csv} from '../modules/write_classification_report.nf'
include { write_single_properties_json } from '../modules/write_properties_json.nf'
include { write_collated_properties_json } from '../modules/write_properties_json.nf'


workflow GENERATE_CLASSIFICATION_REPORT {
    /*
    -----------------------------------------------------------------
    Write Classification Report

    The GENERATE_CLASSIFICATION_REPORT workflow generates some reports
    based on metadata collected during the pipeline:
    - A JSON properties file for each consensus sequence
    - A collated JSON properties files for all consensus sequennce 
    - A collated TSV "classification report" file for all consensii

    -----------------------------------------------------------------
    # Inputs

    - **Metadata Channel**: A channel containing metadata for each
    sample. 

    # Outputs
        - Per-consensus properties.json file channel
        - Collated properties.json file channel
        - Classification report CSV channel (see code below for columns)
    */

    take:
        report_in_ch // meta

    main:
    
        report_in_ch.map{ meta ->
            //
            // Clean up the map to restrict it to simple types that can be written into reports; 
            // This shouldn't be necessary, but since the JSON dumping will fail completely if
            // if map contains any non-simple types, better to be safe than sorry
            //
            def clean 
            clean = { value ->
                if (value instanceof Map) {
                    return value.collectEntries { k, v -> [(k): clean(v)] }.findAll { k, v -> v != null }
                }
                if (value instanceof List) {
                    return value.collect { clean(it) }
                }
                if (value instanceof CharSequence || value instanceof Number || value instanceof Boolean) {
                    return value
                }
                return null
            }
            def cleaned = clean(meta)
            return cleaned
        }.set{ consensus_meta_ch }

        // write individual JSONs, and a collated JSON
        consensus_properties_ch = write_single_properties_json( consensus_meta_ch )
 
        collated_properties_ch = write_collated_properties_json( consensus_meta_ch.collect(), "consensus_sequence_properties" )

        consensus_meta_ch.map { meta ->
            // convert any empty or null values to "None"
            def new_meta = meta.collectEntries { key, value ->
                [(key): (value in [null, ''] ? 'None' : value)]
            }
            get_report_line( new_meta )
        }.collect().set{ report_lines_ch } 

        // Write all of the per-sample report lines to a report file
        classification_report_ch = write_classification_report_csv( 
            get_header_line(), 
            report_lines_ch , 
            "classification_report" 
        )

    emit:
        consensus_properties_ch
        collated_properties_ch
        classification_report_ch 
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


def get_header_line() {

    def sample_headers_1 = "Sample_ID,Virus_Taxon_ID,Virus,Species"
    def sample_headers_2="Reference_Taxon_ID,Selected_Reference,Flu_Segment,Reference_Subtype,Sample_Subtype"
    def qc_headers="Percentage_of_Genome_Covered,Total_Mapped_Reads,Total_Mapped_Bases,Longest_non_N_segment,Percentage_non_N_bases"
    def mut_info_headers = "total_mutations,n_insertions,n_deletions,n_snps,ti_tv_ratio"

    return "${sample_headers_1},${sample_headers_2},${qc_headers},${mut_info_headers}"
}

def get_report_line(meta) {
    /*
    -----------------------------------------------------------------
    Concatenates a list of report lines into a single string.

    - **Input**: A list of strings, where each string is a line of
    the report.
    - **Output**: A single string containing all lines concatenated
    together, separated by newlines.

    -----------------------------------------------------------------

    */
    // HEADERS:
    // sample_info_1:
    // sample_id, virus, report_name, virus_name, 
    // sample_info_2:
    // taxid, reference_selected, flu_segment, virus_subtype, sample_subtype, 
    // qc_info_v:
    // percentage_genome_coverage, total_mapped_reads, longest_no_N_segment, percentage_of_N_bases
    // mut_info_v:
    // total_mutations,n_insertions,n_deletions,n_snps,ti_tv_ratio
    
    def sample_info_1 = "${meta.sample_id},${meta.virus},${meta.report_name},${meta.virus_name}"
    def sample_info_2 = "${meta.taxid},${meta.ref_selected.replace(",","|")},${meta.flu_segment},${meta.virus_subtype},${meta.sample_subtype}"
    def qc_info_v = "${meta.percent_positions_exceeding_depth_10},${meta.reads_mapped},${meta.bases_mapped},${meta.longest_non_n_subsequence},${meta.percent_non_n_bases}"
    def mut_info_v = "${meta.mutations},${meta.insertions},${meta.deletions},${meta.snps},${meta.ti_tv_ratio}"
    return "${sample_info_1},${sample_info_2},${qc_info_v},${mut_info_v}\n"
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
    return errors
}


workflow {
    manifest_channel = Channel.fromPath(params.manifest_file)
    | splitCsv(header: true, sep: ',')
    | map { row ->
        def meta = [[id:row.sample_id,
            taxid:row.taxid,
            ref_selected:row.ref_selected,
            virus_name:row.virus_name,
            virus_subtype:row.virus_subtype,
            flu_segment:row.flu_segment,
            percent_positions_exceeding_depth_10:row.percent_positions_exceeding_depth_10,
            reads_mapped:row.reads_mapped,
            bases_mapped:row.bases_mapped,
            longest_non_n_subsequence:row.longest_non_n_subsequence,
            percent_non_n_bases:row.percent_non_n_bases,
            qc_status:row.qc_status,
            mutations:row.mutations,
            insertions:row.insertions,
            deletions:row.deletions,
            snps:row.snps,
            ti_tv_ratio:row.ti_tv_ratio
        ]]
    }

    GENERATE_CLASSIFICATION_REPORT(manifest_channel)
}

