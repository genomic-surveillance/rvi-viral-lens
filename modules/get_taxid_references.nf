// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.
process get_taxid_reference_files{
    /*
    *         Fetch Fasta Sequence Files for a Given Taxid
    
     The get_taxid_reference_files process is designed to extract 
     reference sequences corresponding to a specified taxonomic ID
     (taxid) from a larger FASTA file. This process retrieves the 
     relevant sequence, writes it to an output file, and indexes 
     it using BWA (Burrows-Wheeler Aligner). This step is essential
     for downstream analysis that requires taxon-specific reference
     sequences. In the context of the pipeline, it is used to 
     extract all the reference files from the Kraken database which
     were observed on that input batch.
    
    * --------------------------------------------------------------
     
    * Input:
       - taxid: The taxonomic ID for which reference sequences are
           being retrieved.
       - kraken_db_library_path: The path to the source FASTA file
           containing sequences for multiple taxa.
    
    * Output:
       - A FASTA file (`${taxid}.fa`) containing the sequence
         corresponding to the specified taxid.
       - BWA index files generated from the FASTA file `(${taxid}.fa.*)`
    
    * NOTE: Output is optional, as sequences may not be found for all 
    *      specified taxids.
    
    * --------------------------------------------------------------
    * > TODO: we need to either remove the optional or, at least,
    *            raise a warning when that happens.
    *
    * --------------------------------------------------------------
    */

    tag "${taxid}"


    input:
        val(taxid)
        val(kraken_db_library_path)
    output:
        tuple val(taxid), path("${taxid}.fa*"), optional : true

    script:
$/
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
                nfound +=1
                found = True
        elif found:
            seq += line

# Write output only if a matching sequence is found
print(f"{nfound} sequences for {taxid_to_find}")
if found:
    ## replace open_paren with hyphen and remove close_paren
    ## this is to avoid those chars being an issue when
    ## loading the subsequently generated bam files into IGV
    header_to_write = header.replace("(", "-").replace(")", "")
    with open(output_file, "w") as output:
        output.write(header_to_write + "\n" + seq + "\n")

    # Run bwa index on the output file
    subprocess.run(["bwa", "index", output_file])
/$
}

/*
# Script Breakdown

Python Script: The script is written in Python and is responsible for
extracting reference sequences for a given taxid from a source FASTA
file and indexing them with BWA.

- Reading the Source File:
  - The script reads through the `source_fna_path`, which contains 
    sequences for multiple taxa.
  - It searches for headers (`>`) that contain the specified taxid
    (`|${taxid}|`).
  - Once a matching sequence is found, it reads the sequence data 
    until the next header or end of the file.

- Output Generation:
  - If a matching sequence is found, it is written to a FASTA file named
    ${taxid}.fa.
  - The script prints the number of sequences found for the given taxid
    for logging and verification purposes.

- BWA Indexing:
  - If a sequence is found and written to the output file, BWA is used to
    index the FASTA file (`bwa index ${taxid}.fa`), preparing it for 
    alignment and further analysis.
*/
