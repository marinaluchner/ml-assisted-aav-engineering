import pandas as pd
import argparse
import glob
import os
from Bio.Seq import Seq

# Parse command-line arguments
parser = argparse.ArgumentParser(
    description="Translate the 'variant' nucleotide sequence column of Excel files into amino acid sequences."
)

parser.add_argument(
    "-i", "--input_dir",
    required=True,
    help="Directory containing input Excel files with a 'variant' column of nucleotide sequences"
)

parser.add_argument(
    "-o", "--output_dir",
    required=True,
    help="Directory to save the output Excel files with the added 'amino_acid_sequence' column"
)

args = parser.parse_args()

os.makedirs(args.output_dir, exist_ok=True)

# Function to translate DNA sequence to amino acid sequence
def translate_dna_to_protein(dna_seq):
    dna_seq_obj = Seq(dna_seq)
    protein_seq = str(dna_seq_obj.translate())
    return protein_seq

# Find all merged_enrichment_scores.xlsx files in the input directory (including subdirectories)
excel_files = glob.glob(os.path.join(args.input_dir, '**', 'merged_enrichment_scores.xlsx'), recursive=True)

output_file = os.path.join(args.output_dir, "merged_enrichment_scores_with_amino_acid_sequences.xlsx")

for file in excel_files:
    # Read the Excel file into a DataFrame
    df = pd.read_excel(file)

    # Check if 'variant' column exists
    if 'variant' in df.columns:
        # Translate the 'variant' column and save in a new column 'amino_acid_sequence'
        df['amino_acid_sequence'] = df['variant'].apply(translate_dna_to_protein)

        # Save the modified DataFrame to the output directory
        print(f"Saving updated file to: {output_file}")
        df.to_excel(output_file, index=False)
