# Nextflow Process Documentation: `run_qc_script`

## Overview

This process performs quality control (QC) on sequence data using a custom Python QC script (available at `bin/` dir).
It uses a set of input files (BAM, FASTA, and reference files) along with the results from samtools mpileup to generate a QC report in CSV format. The key steps include calculating statistics on the BAM file, parsing mpileup output, and running a custom QC script.
This process generates a CSV file with QC metrics and a PNG image visualizing sequencing depth, which are essential for assessing the quality and reliability of the sequencing data. The script is a modified version of [the `qc.py` of the ARTIC pipeline approach](https://gitlab.internal.sanger.ac.uk/malariagen1/ncov2019-artic-nf/-/blob/main/bin/qc.py?ref_type=heads).

## Process

This process runs a QC analysis on the input BAM, FASTA, and reference files, outputting QC metrics and a sequencing depth plot. Below is a detailed breakdown of the process components:

### Tags and Labels

- Tag: `${meta.id}` – The tag identifies the process run by the sample ID, facilitating tracking and logging within the pipeline.

### Input

- Input Tuple: `tuple val(meta), path(bam), path(fasta), path(ref)`
  - `meta`: Metadata associated with the sample, assumes the following keys  are available: 
    - `id`: internal pipeline ID
    - `sample_id`: sample ID
    - `taxid`: taxonomic ID
  - `bam`: Path to the BAM file containing aligned sequencing reads.
  - `fasta`: Path to the consensus FASTA file for the sample.
  - `ref`: Path to the reference file against which the reads are aligned.
  - `samtools_mpileup`: Path to the output file from samtools mpileup that contains coverage information.

### Output

- Output Tuple: `tuple val(meta), path("${meta.id}.qc.csv"), path("${meta.id}.depth.png"), stdout`
  - A CSV file containing QC metrics (`${meta.id}.qc.csv`).
  - A PNG file visualizing the depth of sequencing (`${meta.id}.depth.png`).
  - Standard output (`stdout`) prints the first row of the QC CSV file for quick verification.

### Publish Directory

- Publish Directory: `${params.results_dir}/${meta.sample_id}/${meta.taxid}/`
  - The results are published in a structured directory organized by sample ID and taxonomic ID, making it easy to locate and review QC results.

### Script

```bash
# Generate required samtools flagstat file
samtools flagstat ${bam} > ${samtools_flagstat}

# Parse out the 3 key columns we need from mpileup output file
cat ${samtools_mpileup} | cut -f 1,2,4 > ${mpileup_depths}

# Run QC script
qc.py \
  --outfile ${meta.id}.qc.csv \
  --sample ${meta.id} \
  --ref ${ref} \
  --bam ${bam} \
  --fasta ${fasta} \
  --depths_file ${mpileup_depths} \
  --flagstat_file ${samtools_flagstat} \
  --minimum_depth ${params.qc_minimum_depth} \
  --ivar_md ${params.ivar_min_depth}

# Print first row of output file to stdout
sed -n "2p" ${meta.id}.qc.csv
```

### Script Breakdown

The script performs the following steps:

1. **Generate flagstat report**: Runs samtools flagstat on the BAM file to generate summary statistics, outputting the result to a file (`.flagstat.out`).

2. **Parse mpileup output**: Extracts key columns (reference sequence name, position, and depth) from the samtools mpileup file and saves them to a new file (`.depths.out`).

3. **Run the custom QC script**: Executes the `qc.py` script, passing in the BAM, FASTA, reference, and depths files, along with sample metadata. The script generates the final QC report as a CSV file (`.qc.csv`). Key parameters include:

   - `--outfile`: Path to the output CSV file.
   - `--sample`: Sample ID.
   - `--ref`: Reference file.
   - `--bam`: BAM file.
   - `--fasta`: FASTA file.
   - `--depths_file`: Parsed mpileup depths file.
   - `--flagstat_file`: Samtools flagstat output.
   - `--minimum_depth`: Minimum depth required for quality evaluation (`params.qc_minimum_depth`).
   - `--ivar_md`: Minimum depth threshold for ivar trimming (`params.ivar_min_depth`).

4. **Print output**: The second row of the generated CSV file is printed to standard output using sed for logging purposes.

   - Command: `sed -n "2p" ${meta.id}.qc.csv`
     - This command prints the second line of the QC CSV file, the first row of data (excluding headers). It provides a quick check of the QC output for monitoring and verification purposes.

For more details, check [QC script documentation](../bin/qc.md)
