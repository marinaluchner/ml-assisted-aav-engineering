import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import ShuffleSplit
import seaborn as sns


def LR():
    from sklearn.linear_model import LinearRegression
    return LinearRegression()

def AA_hotencoding(variant):
    
    """
    
    AA_hotencoding takes an amino acid sequence 'variant' of an arbitrary length, 
    and returns a 20xlength one-hot encoding matrix 'onehot_encoded'.   
    
    """
       
    AAs = 'ARNDCQEGHILKMFPSTWYV'
    encoding_length = len(AAs)
    variant_length = len(variant)

    # Define a mapping of chars to integers
    char_to_int = dict((c, i) for i, c in enumerate(AAs))
    int_to_char = dict((i, c) for i, c in enumerate(AAs))

    # Encode input data 
    integer_encoded = [char_to_int[char] for i, char in enumerate(variant) if i <variant_length]
    
    # Start one-hot-encoding
    onehot_encoded = list()
    
    for value in integer_encoded:
        letter = [0 for _ in range(encoding_length)]
        letter[value] = 1
        onehot_encoded.append(letter)
                
    return onehot_encoded
# This takes 6 min to run with 2 folds and 10 test set percentages

# Paul Tol's Notes color scheme
colors = ['#4477AA','#66CCEE','#228833','#CCBB44','#EE6677','#AA3377','#BBBBBB']

# Figure specifications
# Set font type to Arial and default size to 16
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 16


# Load the Excel file into a DataFrame
first_round = 'data/Excel/first_round/enrichment_score_with_amino_acid_sequences.xlsx'
third_round_cutoff = 'data/Excel/third_round/enrichment_score_with_amino_acid_sequence_threshold.xlsx'
third_round = 'data/Excel/third_round/enrichment_score_with_amino_acid_sequence.xlsx'
pseudo_count = 'data/Excel/third_round/enrichment_score_with_amino_acid_sequence_pseudo_count.xlsx'
weighing = 'data/Excel/third_round/enrichment_score_with_amino_acid_sequence_weighing.xlsx'

df_first_round_ = pd.read_excel(first_round)
df_third_round_cutoff_ = pd.read_excel(third_round_cutoff)
df_pseudo_count_ = pd.read_excel(pseudo_count)
df_third_round_ = pd.read_excel(third_round)
df_weighing_ = pd.read_excel(weighing)

# log2 transform the 'averaged_enrichment_score' column
df_first_round_['log2_averaged_enrichment_score'] = df_first_round_['averaged_enrichment_score'].apply(lambda x: np.log2(x))
df_third_round_cutoff_['log2_averaged_enrichment_score'] = df_third_round_cutoff_['averaged_enrichment_score'].apply(lambda x: np.log2(x))
df_pseudo_count_['log2_averaged_enrichment_score'] = df_pseudo_count_['averaged_enrichment_score'].apply(lambda x: np.log2(x))
df_third_round_['log2_averaged_enrichment_score'] = df_third_round_['averaged_enrichment_score'].apply(lambda x: np.log2(x))
df_weighing_['log2_averaged_enrichment_score'] = df_weighing_['averaged_enrichment_score'].apply(lambda x: np.log2(x))

# Define the models to be used
models = {
    'Round_1_LR': [LR, df_first_round_],
    'Round_3_LR_Threshold': [LR, df_third_round_cutoff_],
    'Round_3_LR_Pseudo_count': [LR, df_pseudo_count_],
    'Round_3_LR': [LR, df_third_round_],
    'Round_3_LR_Weight': [LR, df_weighing_]
}

# Define number of rows to leave out for testing
shuffle_split = ShuffleSplit(n_splits=5, test_size=0.1, random_state=0)

# Define percentages
percentages = [100,90,80,70,60,50,40,30,20,10]

# Initialize dictionaries to store results
pearson_results = {model_name: {f'test_set_{percentage}': [] for percentage in percentages} for model_name in models.keys()}

for model_name, model_func in models.items():

    print(f'Model: {model_name}')

    model = model_func[0]
    df = model_func[1]

    for train_index, test_index in shuffle_split.split(df):

        print(f'size train, test: {len(train_index)}, {len(test_index)}')

        train = df.iloc[train_index]
        test = df.iloc[test_index]

        len(train), len(test)
        # Sort the DataFrame by 'averaged_enrichment_score' in descending order
        test_sorted = test.sort_values(by='averaged_enrichment_score', ascending=False)

        # Create a dictionary to hold the different test sets
        test_sets = {}

        # Generate different test set based on top X% of variants
        for percentage in percentages:

            num_rows = int(len(test_sorted) * (percentage / 100))
            test_sets[f'test_set_{percentage}'] = test_sorted.iloc[:num_rows]

            print(f'Test set {percentage}%: {len(test_sets[f"test_set_{percentage}"])} rows')

        for key, test_set in test_sets.items():

            print(f'Test set: {key}')

            # do one hot encoding of amino acids
            train_x = np.asarray([AA_hotencoding(variant) for variant in train["amino_acid_sequence"]])
            train_y = np.asarray([score for score in train['log2_averaged_enrichment_score']])
            test_x = np.asarray([AA_hotencoding(variant) for variant in test_set["amino_acid_sequence"]])
            test_y = np.asarray([score for score in test_set['log2_averaged_enrichment_score']])

            # Reshape the 3D arrays into 2D arrays
            X_train_reshaped = train_x.reshape(train_x.shape[0], -1)
            X_test_reshaped = test_x.reshape(test_x.shape[0], -1)

            # Initialize the Linear Regression model
            lin_reg = LinearRegression()

            # Fit the model on the training data
            lin_reg.fit(X_train_reshaped, train_y)

            # Apply the model to the one-hot encoded test set
            y_pred = lin_reg.predict(X_test_reshaped) 
            pearson_corr = np.corrcoef(test_y, y_pred)[0, 1]

            pearson_results[model_name][key].append(pearson_corr)

# Now pearson_results contains the Pearson correlation coefficients for each model, test set percentage, and fold
print(pearson_results)

# Calculate the average Pearson correlation coefficients for each test set percentage across all folds
average_corrs_per_test_set = {model_name: [] for model_name in models.keys()}
test_set_percentages = sorted([int(key.split('_')[-1]) for key in pearson_results[model_name].keys()])

for model_name in models.keys():
    for percentage in test_set_percentages:
        key = f'test_set_{percentage}'
        average_corrs_per_test_set[model_name].append(np.mean(pearson_results[model_name][key]))

# Prepare data for seaborn
plot_data = []
for model_name in models.keys():
    for percentage in percentages:
        print(percentage)
        key = f'test_set_{percentage}'
        for value in pearson_results[model_name][key]:
            plot_data.append({'Model': model_name, 'Test Set Percentage': percentage, 'Pearson Correlation': value})

plot_df = pd.DataFrame(plot_data)

# Plot the results using seaborn
plt.figure(figsize=(8.25, 6))
sns.lineplot(data=plot_df, x='Test Set Percentage', y='Pearson Correlation', hue='Model', palette=colors, marker='None', errorbar=('ci', 95))

plt.xlabel('Top X% of Test Set', fontsize=16)
plt.ylabel("5-Fold CV test Pearson's correlation", fontsize=16)
# invert the x-axis so it plots from larger to smaller
plt.gca().invert_xaxis()

# legend above plot centered and adjusted for smaller figure
# reduce font and columns, move legend slightly closer to plot for 8x6 figure
plt.legend(fontsize=16, loc='upper center', bbox_to_anchor=(0.5, 1.30), ncol=2, frameon=False)
plt.gca().tick_params(labelsize=16)
plt.ylim(0, 1)
plt.tight_layout()

# save figure; change filename here to preserve the original or create a new plot
output_file = "plots/figure_2c.png"
plt.savefig(output_file)
