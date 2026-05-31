import pandas as pd
import itertools
import os
import matplotlib.pyplot as plt

# Define the folder path
folder_path = "data/Excel/ml_assessment_library/EVAAV"
outputfile_path = "data/Excel/ml_assessment_library/EVAAV"

# Set the global font to be Arial, size 10 (or any other size you prefer)
plt.rcParams['font.size'] = 18
plt.rcParams['font.family'] = 'Arial'

# define numerator and denominator files to calculate enrichment
numerator_files = [["post_encapsulation_rep_1_variant_count.xlsx"], ["post_encapsulation_rep_2_variant_count.xlsx"], ["post_encapsulation_rep_3_variant_count.xlsx"]]
denominator_files = [["pre_encapsulation_rep_1_variant_count.xlsx"], ["pre_encapsulation_rep_2_variant_count.xlsx"], ["pre_encapsulation_rep_3_variant_count.xlsx"]]

# Function to read and combine parts of Excel files
def read_and_combine_excel_parts(base_file_name, file_path):
    combined_df = pd.DataFrame()
    parts_found = False
    
    for part in ['_part1.xlsx', '_part2.xlsx']:
        file_name = base_file_name + part
        full_path = os.path.join(file_path, file_name)
        if os.path.exists(full_path):
            parts_found = True
            print("Found split excel file: ", file_name)
            df_part = pd.read_excel(full_path)
            combined_df = pd.concat([combined_df, df_part], ignore_index=True)
    
    # If no parts were found, try reading the file with just the base_file_name
    if not parts_found:
        print("No split excel files found for: ", base_file_name)
        full_path = os.path.join(file_path, base_file_name + '.xlsx')
        if os.path.exists(full_path):
            combined_df = pd.read_excel(full_path)
    
    return combined_df

for numerator_list, denominator_list in zip(numerator_files, denominator_files):

    for numerator, denominator in itertools.product(numerator_list, denominator_list):
        print(f"Numerator file: {numerator}, Denominator file: {denominator}")

        # get file name without extension
        denominator_filename_without_ext = os.path.splitext(denominator)[0]
        numerator_filename_without_ext = os.path.splitext(numerator)[0]
        print(f"Denominator filename without extension: {denominator_filename_without_ext}")
        print(f"Numerator filename without extension: {numerator_filename_without_ext}")

        # read the excel files into pandas dataframes
        df_numerator = read_and_combine_excel_parts(numerator_filename_without_ext, folder_path)
        df_denominator = read_and_combine_excel_parts(denominator_filename_without_ext, folder_path)
        print("Numerator dataframe shape: ", df_numerator.shape)
        print("Denominator dataframe shape: ", df_denominator.shape)

        # Identify Common Variants
        print(df_numerator.head())
        print(df_denominator.head())
        common_variants_df = pd.merge(df_numerator, df_denominator, on='variant', how='inner')
        num_common_variants = common_variants_df.shape[0]
        print(f"Number of variants in both dataframes: {num_common_variants}")

        # Calculate the Percentage of Variants that are in both dataframes compared to df_numerator
        # e.g., all variants in the EV fraction should be in the OL fraction, therefore if the percentage is not 100%,
        # this means that there is a loss through subsampling
        total_variants_numerator = df_numerator.shape[0]
        percentage_common = (num_common_variants / total_variants_numerator) * 100
        print(f"Percentage of variants in both dataframes compared to df_numerator: {percentage_common:.2f}%")

        # Calculate the sum of all integers in the column "variant_count"
        total_variant_count_denominator = df_denominator['variant_count'].sum()
        total_variant_count_numerator = df_numerator['variant_count'].sum()

        # Normalize the "variant_count" column
        df_denominator['variant_count_normalized'] = df_denominator['variant_count'] / total_variant_count_denominator
        df_numerator['variant_count_normalized'] = df_numerator['variant_count'] / total_variant_count_numerator

        # Sort variant count columns in descending order
        df_denominator = df_denominator.sort_values('variant_count_normalized', ascending=False)
        df_numerator = df_numerator.sort_values('variant_count_normalized', ascending=False)

        # sanity check
        print(df_denominator.head())
        print(df_numerator.head())

        # merging both dataframes
        merged_df = pd.merge(df_denominator, df_numerator, on='variant', suffixes=('_df_denominator', '_df_numerator'))

        # Create a new column for the division result
        merged_df['enrichment_score'] = merged_df['variant_count_normalized_df_numerator'] / merged_df['variant_count_normalized_df_denominator']
        print(merged_df.sort_values('enrichment_score', ascending=False).head())

        # Construct the output file name
        output_file_name = f"enrichment_score_{numerator_filename_without_ext}_{denominator_filename_without_ext}.xlsx"
        output_file = os.path.join(outputfile_path, output_file_name)

        # Save the merged_df as an Excel file
        merged_df.to_excel(output_file, index=False)

        print(f"File saved as {output_file}")