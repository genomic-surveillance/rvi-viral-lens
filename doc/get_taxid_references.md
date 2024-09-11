# Nextflow Process Documentation: `get_taxid_reference_files`

## Overview

The `get_taxid_reference_files` process is designed to extract reference sequences corresponding to a specified taxonomic ID (taxid) from a larger FASTA file. This process retrieves the relevant sequence, writes it to an output file, and indexes it using BWA (Burrows-Wheeler Aligner). This step is essential for downstream analysis that requires taxon-specific reference sequences.  In the context of the pipeline, it is used to extract all the reference files from the Kraken database which were observed on that input batch.

## Process: `get_taxid_reference_files`

The `get_taxid_reference_files` process retrieves sequences for a given taxid from a source FASTA file and indexes them for further analysis. Below is a detailed breakdown of the process components:

### Tags and Labels

- Tag: `${taxid}` – This tag is used to identify the process run by the taxonomic ID, allowing easy tracking and logging within the Nextflow pipeline.

### Input

- Input Values:
  - `val(taxid)`: The taxonomic ID for which reference sequences are being retrieved.
  - `val(kraken_db_library_path)`: The path to the source FASTA file containing sequences for multiple taxa.

### Output

- Output Tuple: `tuple val(taxid), path("${taxid}.fa*"), optional: true`. The outputs include:

  - A FASTA file (`${taxid}.fa`) containing the sequence corresponding to the specified taxid.

  - BWA index files generated from the FASTA file `(${taxid}.fa.*)`.

Output is optional, as sequences may not be found for all specified taxids.

> NOTE: we need to either remove the optional or, at least, raise a warning when that happens.

### Publish Directory

- Publish Directory: `${params.results_dir}/reference_files/${taxid}/`.
  - The results are published to a directory organized by taxid, facilitating easy access to reference files specific to each taxonomic group.

### Script

The script section includes the Python script executed within the process:

```python
#!/usr/bin/python3

import subprocess 

taxid_to_find = "${taxid}"
output_file = "${taxid}.fa"
source_fna_path = "${kraken_db_library_path}"

with open(source_fna_path, "r") as source_file:
    header = ""
    seq = ""
    found = False
    nfound = 0
    for line in source_file:
        line = line.strip()
        if line.startswith(">"):
            if found:
                break  # Exit loop after finding the first sequence matching taxid
            header = line
            if f"|{taxid_to_find}|" in header:
                nfound += 1
                found = True
        elif found:
            seq += line

# Write output only if a matching sequence is found
print(f"{nfound} sequences for {taxid_to_find}")
if found:
    with open(output_file, "w") as output:
        output.write(header + "\n" + seq + "\n")

    # Run bwa index on the output file
    subprocess.run(["bwa", "index", output_file])
```

### Script Breakdown

Python Script: The script is written in Python and is responsible for extracting reference sequences for a given taxid from a source FASTA file and indexing them with BWA.

- Reading the Source File:
  - The script reads through the `source_fna_path`, which contains sequences for multiple taxa.
  - It searches for headers (`>`) that contain the specified taxid (`|${taxid}|`).
  - Once a matching sequence is found, it reads the sequence data until the next header or end of the file.

- Output Generation:
  - If a matching sequence is found, it is written to a FASTA file named ${taxid}.fa.
  - The script prints the number of sequences found for the given taxid for logging and verification purposes.

- BWA Indexing:
  - If a sequence is found and written to the output file, BWA is used to index the FASTA file (`bwa index ${taxid}.fa`), preparing it for alignment and further analysis.
