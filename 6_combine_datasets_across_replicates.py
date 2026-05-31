import pandas as pd
import glob
import os
import sys

# Directory paths
rep1 = 'data/Excel/ml_assessment_library/EVAAV/enrichment_score_pre_encapsulation_rep_1_variant_count_post_encapsulation_rep_1_variant_count.xlsx'
rep2 = 'data/Excel/ml_assessment_library/EVAAV/enrichment_score_pre_encapsulation_rep_2_variant_count_post_encapsulation_rep_2_variant_count.xlsx'
rep3 = 'data/Excel/ml_assessment_library/EVAAV/enrichment_score_pre_encapsulation_rep_3_variant_count_post_encapsulation_rep_3_variant_count.xlsx'

df1 = pd.read_excel(rep1, usecols=['variant', 'enrichment_score'])
df2 = pd.read_excel(rep2, usecols=['variant', 'enrichment_score'])
df3 = pd.read_excel(rep3, usecols=['variant', 'enrichment_score'])

# Merge the first two DataFrames
merged_df_1_2 = pd.merge(df1, df2, on='variant', suffixes=('_rep_1', ''))
print(merged_df_1_2)


# Merge the result with the third DataFrame
merged_df = pd.merge(merged_df_1_2, df3, on='variant', suffixes=('_rep_2', '_rep_3'))
print(merged_df)

# Calculate the average enrichment score
merged_df['average_enrichment_score'] = merged_df[['enrichment_score_rep_1', 'enrichment_score_rep_2', 'enrichment_score_rep_3']].mean(axis=1)

# Calculate the standard deviation of the enrichment scores
merged_df['std_dev_enrichment_score'] = merged_df[['enrichment_score_rep_1', 'enrichment_score_rep_2', 'enrichment_score_rep_3']].std(axis=1)

# Save the merged DataFrame to an Excel file
merged_file_name = 'data/Excel/merged_enrichment_scores.xlsx'
merged_df.to_excel(merged_file_name, index=False)