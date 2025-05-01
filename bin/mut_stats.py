#!/usr/bin/env python3
# Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.
import csv
import argparse
from collections import defaultdict

def calculate_mutation_stats(tsv_file):
    # Initialize counters
    stats = {
        'total_mutations': 0,
        'insertions': 0,
        'deletions': 0,
        'snps': 0,
        'transitions': 0,
        'transversions': 0,
        "ti_tv_ratio": 0
    }
    
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

    with open(tsv_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        
        # for each row, skip any which PASS equal False
        # count mutations
        for row in reader:
            # Skip if not PASSing
            if row['PASS'] != 'TRUE':
                continue
            
            stats['total_mutations'] += 1
            ref = row['REF'].upper()
            alt = row['ALT'].upper()
            
            # Check for insertion/deletion/SNP
            # assuming insertion being represented as:
            #   REF ALT
            #   A   AT
            #
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
                    stats['transitions'] += 1
                elif (ref, alt) in transversion:
                    stats['transversions'] += 1
    
    # Calculate Ti/Tv ratio
    try:
        ti_tv_ratio = stats['transitions'] / stats['transversions']
    except ZeroDivisionError:
        ti_tv_ratio = float('inf')

    # Print results
    print("Mutation Statistics (PASS=True only):")
    print(f"Total mutations: {stats['total_mutations']}")
    print(f"Insertions: {stats['insertions']}")
    print(f"Deletions: {stats['deletions']}")
    print(f"SNPs: {stats['snps']}")
    print(f"Transitions (Ti): {stats['transitions']}")
    print(f"Transversions (Tv): {stats['transversions']}")
    print(f"Ti/Tv ratio: {ti_tv_ratio:.2f}")
    
    stats["ti_tv_ratio"] = round(ti_tv_ratio,2)
    
    return stats
def main():
    parser = argparse.ArgumentParser(description='Calculate mutation statistics from TSV file.')
    parser.add_argument('input_file', help='Path to the input TSV file')
    args = parser.parse_args()
    
    calculate_mutation_stats(args.input_file)

if __name__ == "__main__":
    main()