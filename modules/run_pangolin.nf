process run_pangolin {
    /*
    *                Subtype SCOV2 sequence

    The run_pangolin process in this Nextflow pipeline is designed to
    classify SARS-CoV-2 sequences into lineages using the Pangolin
    tool. Pangolin (Phylogenetic Assignment of Named Global Outbreak
    LINeages) is a tool widely used in the analysis of SARS-CoV-2 
    genomes to determine the likely lineage of the virus based on the
    consensus sequence.
    
    Pangolin documentation:
    https://cov-lineages.org/resources/pangolin.html

    * ---------------------------------------------------------------
    * Input
        - `meta`: Metadata associated with the sample, including 
            identifiers like sample ID and taxonomic ID.
        - `mapped_fasta`: The path to the consensus FASTA file
            generated for the sample.

    * Output
        - `mapped_fasta`: The path to the input consensus FASTA file.
        - `env(lineage)`: Environment variables capturing the lineage
            assignment information, which are extracted from the
            Pangolin output.

    * ---------------------------------------------------------------
    */

    publishDir "${params.results_dir}/${meta.sample_id}/${meta.taxid}/", mode: "copy", pattern: "*.csv"
    label "pangolin"
    tag "${meta.id}"

    input:
        tuple val(meta), path(mapped_fasta)

    output:
        tuple val(meta), path(mapped_fasta), env(lineage)

    shell:
        lineage_report = "${meta.id}_lineage.csv"
        consensus_fasta = mapped_fasta
        '''
        # run pangolin
        pangolin !{consensus_fasta} --outfile !{lineage_report}

        # get sample lineage assignment
        # stores every value of every column as an env variable
        # PS: assumes only one row at lineage_report
        IFS=, read -r taxon lineage conflict ambiguity_score scorpio_call \
        scorpio_support scorpio_conflict scorpio_notes version \
        pangolin_version scorpio_version constellation_version \
        is_designated qc_status qc_notes note < <(tail -n +2 !{lineage_report})

        echo "Taxon: $taxon"
        echo "Lineage: $lineage"
        echo "Conflict: $conflict"
        '''
/*
# Script Breakdown

1. **Run Pangolin**:

    - Command:
    ```
    pangolin "${consensus_fasta}" --outfile "${lineage_report}"`
    ```

    - This command runs Pangolin on the input consensus FASTA file 
    and writes the results to a CSV file named with the sample ID
    and `_lineage.csv` suffix.

2. **Extract Lineage Information**:

    - The script uses the IFS (Internal Field Separator) command to 
    parse the CSV output file from Pangolin, extracting values from
    each column into environment variables.

    - The script assumes that there is only one row of data 
    (excluding the header) in the lineage report and reads this data
    into variables such as taxon, lineage, conflict, ambiguity_score,
    and others.

3. **Log Lineage Details**:

    - Lineage details such as taxon, lineage, and conflict are echoed 
    to in order to be captured as strings on the output channel.

    - Curretly, only `$lineage` is used for downstream analysis.

---------------------------------------------------------------------

> TODO: 
    -The `Taxon` and `Conflict` values are echoed in order to
have those value within easy access on the `.command.log` files, it
was usefull for debugging during development, but may be not necessary
to keep it anymore. We should remove it.
    - If for any reason the consensus fasta have more than one sequence,
we would be able to get results for the first sequence only. At current
implementation, this should never happen because ivar only output one 
consensus sequences per taxid sample combination. If this ever change,
we will need to adapt this process.
*/
}