import pandas as pd
import argparse
import os
from functools import reduce

# Parse command-line arguments
parser = argparse.ArgumentParser(
    description="Combine enrichment score datasets across replicates and calculate the average and standard deviation."
)

parser.add_argument(
    "-i", "--inputs",
    required=True,
    nargs='+',
    help="Input Excel files with enrichment scores for each replicate (at least 2 required)"
)

parser.add_argument(
    "-o", "--output",
    required=True,
    help="Output Excel file path for the merged enrichment scores"
)

args = parser.parse_args()

# Allow --output to be either a directory (a default filename is appended)
# or a full file path ending in .xlsx
output_path = args.output
if os.path.splitext(output_path)[1] == '':
    os.makedirs(output_path, exist_ok=True)
    output_path = os.path.join(output_path, 'merged_enrichment_scores.xlsx')

if len(args.inputs) < 2:
    parser.error("At least 2 input files are required to combine across replicates.")

# Read each replicate file, renaming enrichment_score to keep replicate identity
dfs = []
for i, input_file in enumerate(args.inputs, start=1):
    df = pd.read_excel(input_file, usecols=['variant', 'enrichment_score'])
    df = df.rename(columns={'enrichment_score': f'enrichment_score_rep_{i}'})
    dfs.append(df)
    print(f"Replicate {i} file: {input_file}, shape: {df.shape}")

# Merge all replicate DataFrames on 'variant'
merged_df = reduce(lambda left, right: pd.merge(left, right, on='variant'), dfs)
print(merged_df)

# Calculate the average and standard deviation of the enrichment scores across replicates
enrichment_score_columns = [f'enrichment_score_rep_{i}' for i in range(1, len(args.inputs) + 1)]
merged_df['average_enrichment_score'] = merged_df[enrichment_score_columns].mean(axis=1)
merged_df['std_dev_enrichment_score'] = merged_df[enrichment_score_columns].std(axis=1)

# Save the merged DataFrame to an Excel file
merged_df.to_excel(output_path, index=False)
print(f"File saved as {output_path}")
