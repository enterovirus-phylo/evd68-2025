#!/usr/bin/env python3
import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser(
        description='Add source_2025 column for sequences from 2025 ENPEN data'
    )
    parser.add_argument('--metadata', required=True, help='Input metadata file')
    parser.add_argument('--matching-file', required=True, help='Pre-made matching file (strain -> site)')
    parser.add_argument('--output', required=True, help='Output metadata file')

    args = parser.parse_args()

    # Read files
    meta = pd.read_csv(args.metadata, sep='\t')

    # Read the matching file (tab-separated: strain -> site)
    matching = pd.read_csv(args.matching_file, sep='\t', header=None, names=['strain', 'site'])

    # Create a mapping of strain names to source_2025 values
    source_2025_map = {}
    for _, row in matching.iterrows():
        strain = row['strain']
        site = row['site']
        source_2025_map[strain] = f"{site} 2025"

    # Add source_2025 column to metadata
    if 'source_2025' not in meta.columns:
        meta['source_2025'] = ''

    # Fill in the source_2025 values for matching strains
    for strain, source_value in source_2025_map.items():
        meta.loc[meta['strain'] == strain, 'source_2025'] = source_value

    # Write output
    meta.to_csv(args.output, sep='\t', index=False)
    print(f"Added source_2025 column to {len(source_2025_map)} sequences")
    print(f"Output written to {args.output}")

if __name__ == '__main__':
    main()
