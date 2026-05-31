import pandas as pd
import glob
from Bio.Seq import Seq

# Function to translate DNA sequence to amino acid sequence
def translate_dna_to_protein(dna_seq):
    dna_seq_obj = Seq(dna_seq)
    protein_seq = str(dna_seq_obj.translate())
    return protein_seq

# Find all Excel files matching the pattern
excel_files = glob.glob('data/Excel/ml_assessment_library/EVAAV/merged_enrichment_scores.xlsx')

for file in excel_files:
    # Read the Excel file into a DataFrame
    df = pd.read_excel(file)
    
    # Check if 'variant' column exists
    if 'variant' in df.columns:
        # Translate the 'variant' column and save in a new column 'amino_acid_sequence'
        df['amino_acid_sequence'] = df['variant'].apply(translate_dna_to_protein)
        
        # Save the modified DataFrame back to an Excel file
        # You can choose to overwrite the original file or save it as a new file
        new_file_name = file.replace('.xlsx', '_with_aa_column.xlsx')

        print(f"Saving updated file to: {new_file_name}")
        df.to_excel(new_file_name, index=False)