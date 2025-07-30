#!/usr/bin/env python3
# Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.
import re
import json
import pandas as pd
import argparse

def get_args():
    parser = argparse.ArgumentParser("Script to convert kraken2ref JSON to TSV report.")
    parser.add_argument("-i", "--input_json", type = str, required = True, help = "The JSON file produced by kraken2ref.")
    parser.add_argument("-r", "--report", type = str, required = True, help = "The kraken2 taxonomic report.")
    parser.add_argument("-out_suffix", type = str, required=False, default = ".viral_pipe.report.tsv", 
        help="report filname suffix (default: '.viral_pipe.report.tsv'")

    return parser

def regex_subtyping(regex, search_string):
    """
    Searches for the first occurrence of the given regular expression in the provided string.

    Args:
        regex (str): The regular expression to be searched.
        search_string (str): The string to be searched.

    Returns:
        str: The entire matched string if found, otherwise None.

    ---------------------------------------------------------------------------
    TODOs:
      - either replace this with a function that returns subtype or remove it
      and replace the code that calls it with a single regex that uses
      capture groups to parse subtype and segment data from isolate name. This
      can be taken from kraken-flu utils for example.
    """

    pattern = re.compile(regex)
    found = pattern.search(search_string)
    if found:
        return found.group()
    else:
        return None


def get_report(in_file, report_file, out_suffix=".viral_pipe.report.tsv"):

    """
    Generate a viral report by processing data from a Kraken report and a reference JSON file.

    This function reads a Kraken classification report and a JSON file containing metadata and taxonomic 
    information. It processes the data to extract information about selected reference taxa, including their 
    species, virus names, subtypes (for influenza), and the number of reads per taxon. The processed data is 
    output as a tab-separated values (TSV) file, which contains detailed information about each selected 
    reference taxon.

    Args:
        in_file (str): Path to the input JSON file containing reference taxonomic metadata and results.
        report_file (str): Path to the Kraken report file in TSV format.
        out_suffix (str, optional): The suffix for the output report file. Defaults to ".viral_pipe.report.tsv".

    Returns:
        None: The function writes the output report to a TSV file in the current working directory. The file name 
              is generated using the sample ID from the input JSON file and the `out_suffix`.

    Output Structure:
        The generated report includes the following fields for each selected reference:
        - sample_id: The ID of the sample being processed.
        - virus: The source taxonomic ID of the virus species.
        - virus_name: Human-readable name of the virus species.
        - selected_taxid: Taxonomic ID of the selected reference.
        - ref_selected: Name of the selected reference.
        - sample_subtype: Influenza subtype if identified from the sample (specific to Flu Segments 4 and 6).
        - flu_segment: Segment number for influenza viruses if subtype data is available.
        - virus_subtype: General influenza subtype (if applicable).
        - parent_selected: Taxonomic ID of the parent taxon selected in the reference JSON.
        - num_reads: Number of reads assigned to the selected taxon.
        - report_name: General name for the virus (e.g., "Influenza A Virus" for Alphainfluenzavirus).

    The function handles influenza-specific logic to extract subtype information based on virus segment 
    identifiers (for H and N subtypes). The report will include influenza subtype information if it is relevant 
    to the selected reference taxa.

    ---------------------------------------------------------------------------
    TODOs:
        - This function has an unfortunate name. It would be far clearer what it
          does if it would be called write_report. I would expect a get_report
          function to return data but it doesn't.
        - we should return data rather than write straight to file whenever
          possible. The function is hard to test if it writes to file because
          the test would need to read and parse the file instead of working 
          with data.
    """

    ## read in kraken report and kraken2ref JSON
    report_df = pd.read_csv(report_file, sep = "\t", header = None)
    ref_json = json.load(open(in_file))

    ## collect tax_ids for selected reference taxa
    selected_ref_taxa = ref_json["metadata"]["selected"]

    ## collect tax_ids for Species associated with each chosen reference tax_id
    virus_ids = [int(ref_json["outputs"][k]["source_taxid"]) for k in ref_json["outputs"].keys()]

    ## collect the human-readable names associated with the viru_ids and selected_ref_taxa
    ## from kraken report
    virus_names_df = report_df[report_df[4].isin(virus_ids)]
    virus_desc_names_dict = dict(zip(virus_names_df[4], [i.strip() for i in virus_names_df[5]]))

    subset_df = report_df[report_df[4].isin(selected_ref_taxa)]
    desc_names_dict = dict(zip(subset_df[4], [i.strip() for i in subset_df[5]]))

    ## get the sample_id
    sample_id = ref_json["metadata"]["sample"]

    ## initialise output data structure
    report_output = {}
    for idx, selected_ref in enumerate(selected_ref_taxa):
        ## get data specific to the selected_ref
        selected_data_dict = ref_json["outputs"][str(selected_ref)]

        ## identify its Species and get Species name
        virus = selected_data_dict["source_taxid"]
        virus_name = None
        virus_name = virus_desc_names_dict[virus]

        ## get the name of the actual reference that was chosen
        ref_name = desc_names_dict[selected_ref]

        ## if selected_ref is ANY flu segment, collect its subtype info
        ## THIS IS THE SUBTYPE OF THE CHOSEN REFERENCE ONLY
        ## THIS MAY NOT BE THE SAME AS THE SUBTYOE OF THE SAMPLE ITSELF
        subtype = None
        segment = None
        generic_subtype = None

        # if flu B, get segment number, but not subtype
        # NOTE:
        #  This part of the code assumes flu B refname will ALWAYS
        #  have the the term "segment <int>" in its header 

        if ref_name.startswith("B/"):
            segment = regex_subtyping("(?<=segment )[0-9]", ref_name)

        else:
            # if flu A, generic subtypy will not be None, get segment and subtype
            ## if selected_ref is Flu Segment 4 or 6, collect subtype info
            ## THIS WILL BE THE SUBTYPE OF THE SAMPLE
            generic_subtype = regex_subtyping("H[0-9]+N[0-9]+", ref_name)

            if generic_subtype is not None:
                segment = regex_subtyping("(?<=segment )[0-9]", ref_name)
            if segment in [4, "4"]:
                subtype = regex_subtyping("H[0-9]+", ref_name)
            if segment in [6, "6"]:
                subtype = regex_subtyping("N[0-9]+", ref_name)

        ## collect number of reads written to each fastq filepair
        if "per_taxon" in ref_json["metadata"]["summary"].keys():
            num_reads = ref_json["metadata"]["summary"]["per_taxon"][str(selected_ref)]
        else:
            num_reads = None

        if "Alphainfluenzavirus" in virus_name:
            report_name = "Influenza A Virus"
        elif "Betainfluenzavirus" in virus_name:
            report_name = "Influenza B Virus"
        else:
            report_name = virus_name

        ## populate output data dict
        report_output[idx] = {
                                "sample_id": sample_id,
                                "virus": virus,
                                "virus_name": virus_name,
                                "selected_taxid": selected_ref,
                                "ref_selected": ref_name,
                                "sample_subtype": subtype,
                                "flu_segment": segment,
                                "virus_subtype": generic_subtype,
                                "parent_selected": selected_data_dict["parent_selected"],
                                "num_reads": num_reads,
                                "report_name": report_name
                            }

    ## once iteration over all chosen refs is done
    ## convert dict of dicts to tabular format and write to tsv
    output_df = pd.DataFrame.from_dict(report_output.values())
    # if not empty, sort dataframe and write csv
    if len(output_df) > 0:
        output_df.sort_values("selected_taxid").to_csv(f"{sample_id}{out_suffix}", sep = "\t", header=True, index = False)
    # if empty, write a csv file
    else:
        output_df.to_csv(f"{sample_id}{out_suffix}", sep = "\t", header=True, index = False)

def main():
    args = get_args().parse_args()
    get_report(args.input_json, args.report, out_suffix=args.out_suffix)

if __name__ == "__main__":
    main()
