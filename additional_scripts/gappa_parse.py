#! /usr/bin/env python
import argparse
import pandas as pd

def process_taxonomic_data(input_file, output_file):
    # Load input file
    df = pd.read_csv(input_file, sep='\t', header=None)
    
    # Split last column into separate taxa levels
    taxon_levels = df.iloc[:, -1].str.split(';', expand=True)
    
    # Setup out df
    output_df = pd.DataFrame(columns=['taxon_level', 'taxon_label', 'count_of_occurrences', 'percent_of_total_calls'])
    
    # Process each taxa level
    for i, column in enumerate(taxon_levels.columns):
        counts = taxon_levels[column].fillna('Missing').value_counts(dropna=False)
        percentages = (counts / counts.sum()) * 100
        level_data = pd.DataFrame({
            'taxon_level': f'Level_{i+1}',
            'taxon_label': counts.index,
            'count_of_occurrences': counts.values,
            'percent_of_total_calls': percentages.values
        })
        output_df = pd.concat([output_df, level_data], ignore_index=True)
    
    # df to CSV
    output_df.to_csv(output_file, index=False)

if __name__ == "__main__":
    # arg parsing
    parser = argparse.ArgumentParser(description="Process single gene gappa assign best hit outputs and generate a summary CSV")
    parser.add_argument('-i', '--input', required=True, help="Input file path (TSV format)")
    parser.add_argument('-o', '--output', required=True, help="Output file path (CSV format)")
    
    args = parser.parse_args()
    
    # Do the thing
    process_taxonomic_data(args.input, args.output)

