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
        'transversions': 0
    }
    
    # Define transition pairs
    transitions = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}

    with open(tsv_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        
        # for each row, skip any which PASS equal False
        # count mutations
        for row in reader:
            # Skip if not PASSing
            if row['PASS'] != 'TRUE':
                continue
            
            stats['total_mutations'] += 1
            ref = row['REF']
            alt = row['ALT']
            
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
                else:
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

def main():
    parser = argparse.ArgumentParser(description='Calculate mutation statistics from TSV file.')
    parser.add_argument('input_file', help='Path to the input TSV file')
    args = parser.parse_args()
    
    calculate_mutation_stats(args.input_file)

if __name__ == "__main__":
    main()