#!/usr/bin/env python3
# Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.
import json
import csv
import re
import argparse

############
def collect_bam_header_stats(stats:dict, bam_header_file:str):
    """
    Parses the LN attribute from the SQ line of the BAM header
    in order to obtain a value for the reference_length property
    """
    seq_len = 0    
    with open(bam_header_file, 'r') as f:
        for line in f:
            match = re.match(r'^@SQ\s+.+\s+LN:(\d+)', line)
            if match:
                seq_len = int(match.group(1))

    stats['reference_length'] = seq_len

############
def collect_mutation_stats(stats:dict, ivar_variants_file:str):
    """
    Obtains a collection of properties (counts of SNPs, indels etc)
    from the variants file produced by 'ivar variants' file. 
    Format of ivar variants file documented here: https://andersen-lab.github.io/ivar/html/manualpage.html
    """

    # Define transition pairs
    transitions = {
        ('A', 'G'), ('G', 'A'), # Purine <-> Purine
        ('C', 'T'), ('T', 'C')  # Pyrimidine <-> Pyrimidine
    }
    transversion = {
        ('A', 'T'), ('T', 'A'), # Purine <-> Pyrimidine
        ('A', 'C'), ('C', 'A'),
        ('G', 'T'), ('T', 'G'),
        ('G', 'C'), ('C', 'G')
    }

    with open(ivar_variants_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        
        tis = 0
        tvs = 0

        # count mutations
        for row in reader:
            # Skip if not PASSing
            if row['PASS'] != 'TRUE':
                continue
            
            stats['mutations'] += 1
            ref = row['REF'].upper()
            alt = row['ALT'].upper()
            
            # Check for insertion/deletion/SNP
            # assuming insertion being represented as:
            #   REF ALT
            #   A   AT
            # and deletions as:
            #   REF ALT
            #   AT  A
            # 
            if len(ref) < len(alt):
                stats['insertions'] += 1
            elif len(ref) > len(alt):
                stats['deletions'] += 1
            else:
                stats['snps'] += 1
                # Check for transition/transversion
                if (ref, alt) in transitions:
                    tis += 1
                elif (ref, alt) in transversion:
                    tvs += 1
    
    stats['ti_tv_ratio'] = f"{tis}/{tvs}"
 

############ 
def collect_read_stats( stats: dict, flagstat_file: str):
    """
    Parses read and base counts from the samtools flagstat file. 
    Example of samtools flagstat output:
    
    6326 + 0 in total (QC-passed reads + QC-failed reads)
    5772 + 0 primary
    0 + 0 secondary
    554 + 0 supplementary
    0 + 0 duplicates
    0 + 0 primary duplicates
    6326 + 0 mapped (100.00% : N/A)
    5772 + 0 primary mapped (100.00% : N/A)
    5772 + 0 paired in sequencing
    2886 + 0 read1
    2886 + 0 read2
    5460 + 0 properly paired (94.59% : N/A)
    5772 + 0 with itself and mate mapped
    0 + 0 singletons (0.00% : N/A)
    0 + 0 with mate mapped to a different chr
    0 + 0 with mate mapped to a different chr (mapQ>=5)
    """
    with open(flagstat_file) as flagstat_reader:
        flagstat_values = [line.split(" ")[0] for line in flagstat_reader.readlines()]

    # "total" is the count of alignments; "primary" is the count of reads
    total_reads = int(flagstat_values[1])

    stats['reads_mapped'] = int(flagstat_values[7])
    stats['reads_mapped_in_proper_pairs'] = int(flagstat_values[11])
    stats['reads_unmapped'] = total_reads - int(stats['reads_mapped'])

############
def collect_depth_stats( stats: dict, depths_file: str ):
    """
    Parses and calculates a collection of properties using
    data from the output of 'samtools depth'. 
    Format of aamtools depth output: <sequence name>  <position>  <coverage>
    """
    total_aligned_bases = 0
    positions_exceeding_min_depth = 0
    for i in range(0, 101, 5):
        stats['positions_exceeding_depth'][str(i)] = 0

    with open(depths_file, "r") as depth_file_reader:
        for line in depth_file_reader:
            parts = line.rstrip().split()
            depth = int(parts[2])
            total_aligned_bases += depth
            for i in range(0, 101, 5):
                if (int(depth) >= i):
                    stats['positions_exceeding_depth'][str(i)] += 1

    stats['percent_positions_exceeding_depth_10'] = round(100 * int(stats['positions_exceeding_depth']["10"]) / int(stats['reference_length']), 2)
    stats['bases_mapped'] = total_aligned_bases
    stats['mean_depth_per_position'] = int(total_aligned_bases / int(stats['reference_length']))

############
def collect_consensus_sequence_stats( stats: dict, fasta_file : str):
    """
    Obtains metrics from the consensus sequence itself
    """
    pct_non_N_bases = 0
    longest_non_N_stretch = 0

    # Read in FASTA and if not empty then perform QC calculations
    seq = parse_fasta_to_string( fasta_file )
    if len(seq) != 0:
        pct_non_N_bases = get_pct_non_N_bases(seq)
        longest_non_N_stretch = get_longest_non_n_stretch(seq)

    stats['consensus_length'] = len(seq)
    stats['percent_non_n_bases'] = round(pct_non_N_bases, 2)
    stats['longest_non_n_subsequence'] = longest_non_N_stretch

############
def get_pct_non_N_bases(seq : str) -> float:
    """
    Determines the % of bases of a FASTA that are 'n'.
    """
    n_count = 0
    non_n_count = 0

    for char in seq:
        if char != 'n':
            non_n_count += 1
        else:
            n_count += 1
       
    pct_non_n_bases = (non_n_count / len(seq) * 100)

    return pct_non_n_bases

############
def get_longest_non_n_stretch(seq : str) -> int:
    """
    Returns the length of the longest stretch of non-'n' sequence
    """
    max_len = 0
    current_len = 0

    for char in seq:
        if char == 'n':
            current_len = 0  # reset if the letter is found
        else:
            current_len += 1
            max_len = max(max_len, current_len)

    return max_len


############
def parse_fasta_to_string(filepath):
    """
    Parses a single-entry fasta file into a string
    """
    sequence = []

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('>'):
                continue  # Skip header and blank lines
            sequence.append(line.lower())

    return ''.join(sequence)

############
def generate_qc_file(args: argparse.ArgumentParser.parse_args):
    stats = {
        # ivar mutation stats
        'mutations' : 0,
        'insertions': 0,
        'deletions': 0,
        'snps': 0,
        "ti_tv_ratio": "0/0",
        # read alignment stats
        "reads_mapped" : 0,
        "bases_mapped" : 0,
        "reads_mapped_in_proper_pairs" : 0,
        # coverage stats
        "reference_length" : 0,
        "positions_exceeding_depth" : {},
        "percent_positions_exceeding_depth_10" : 0.00,
        # consensus stats
        "percent_non_n_bases" : 0.00,
        "longest_non_n_subsequence" : 0
    }

    collect_bam_header_stats( stats, args.samtools_bam_header_file )
    collect_mutation_stats( stats, args.ivar_variants_file )
    collect_read_stats( stats, args.samtools_flagstat_file)
    collect_depth_stats ( stats, args.samtools_depth_file )
    collect_consensus_sequence_stats( stats, args.fasta_file)

    # Write output header & QC columns to a JSON file
    with open(args.outfile, 'w') as jsonfile:
        json.dump( stats, jsonfile, indent=4)
       
############
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--outfile', required=True, type=str,
        help='''The path of the output QC summary file''')
    parser.add_argument('--fasta_file', required=True, type=str,
        help='''The path of a consensus fasta file produced by ivar.''')
    parser.add_argument('--samtools_bam_header_file', required=True, type=str,
        help='''Output of samtools view -H on the BAM file source of the consensus sequence''')
    parser.add_argument('--samtools_depth_file', required=True, type=str,
        help='''Output of samtools depth on source BAM file for the consensus.''')
    parser.add_argument('--samtools_flagstat_file', required=True, type=str,
        help='''Output of samtools flagstat on source BAM file for the consensus''')
    parser.add_argument('--ivar_variants_file', required=True, type=str,
        help='''Output of ivar variants command, to allow summary mutation information to be reported''')
  
    args = parser.parse_args()
    generate_qc_file(args)

if __name__ == "__main__":
    main()
