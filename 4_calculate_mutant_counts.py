#Take fastq file (containing AAV variants as lines) as inputfile and load as pandas df
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
import glob
import os
from tqdm import tqdm
import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser(
    description="Calculate mutant counts from NGS reads."
)

parser.add_argument(
    "-i", "--input",
    required=True,
    help="Input directory containing sample folders"
)

parser.add_argument(
    "-o", "--output",
    required=True,
    help="Output directory for mutant counts"
)

args = parser.parse_args()
input_folder = args.input
output_folder = args.output

# loop through filtered and trimmed files
def fastq_to_dataframe(fastq_file):

    records = []
    for record in tqdm(SeqIO.parse(fastq_file, "fastq"), total=total, desc="Parsing FASTQ"):
        records.append({
            "ID": record.id,
            "variant": str(record.seq),
            "quality": record.letter_annotations["phred_quality"]
        })

    df = pd.DataFrame(records)
    return df

for file in glob.glob(os.path.join(input_folder, "*filtered_trimmed.fastq")):
    base = os.path.basename(file).replace('_filtered_trimmed.fastq', '')
    print(base)

    # import fastq into pandas dataframe
    df = fastq_to_dataframe(file)

    # calculate the length of each variant
    df["variant_length"] = df["variant"].apply(len)

    # filter by variants of the right length
    df_cleaned = df[df['variant_length'].between(21, 63)]

    # table variants and return variant count
    df_cleaned['variant_count'] = df_cleaned.groupby('variant')['variant'].transform('count')
    df_cleaned = df_cleaned.sort_values('variant_count', ascending=False)
    print("Shape of cleaned dataframe: ", df_cleaned.shape)
    print(df_cleaned.head())
    df_value_counts_no_duplicates = df_cleaned.drop_duplicates(subset='variant')
    df_value_counts_no_duplicates = df_value_counts_no_duplicates.sort_values('variant_count', ascending=False)
    print("Shape of df_value_counts_no_duplicates: ", df_value_counts_no_duplicates.shape)
    print(df_value_counts_no_duplicates.head())

    try:
        # Try to save the DataFrame as a whole
        df_value_counts_no_duplicates.to_excel(f'{output_folder}/{base}_variant_count.xlsx', index=False)
        
    except Exception as e:
        print(f"Error saving DataFrame: {e}")
        print("Splitting DataFrame into parts...")

        # Calculate the index to split on
        split_index = len(df_value_counts_no_duplicates) // 2

        # Split the DataFrame into two
        df1 = df_value_counts_no_duplicates.iloc[:split_index]
        df2 = df_value_counts_no_duplicates.iloc[split_index:]

        # Save the DataFrames to Excel files
        df1.to_excel(f'{output_folder}/{base}_variant_count_part1.xlsx', index=False)
        df2.to_excel(f'{output_folder}/{base}_variant_count_part2.xlsx', index=False)

    # print the maximum variant count
    max_count = df_value_counts_no_duplicates['variant_count'].max()
    print("Maximum variant count: ", max_count)
