import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import ShuffleSplit
import seaborn as sns
import argparse

def LongShortTermMemoryModel(L1=160, L2=20):
    
    """
    
    parent_model builds an LSTM model with paramters that work accross all functional fitness models in the Fit4Function study. 
    L1 and L2 define the sizes of the model two layers. 
    
    """
    import os
    os.environ['KERAS_BACKEND'] = 'tensorflow'
    from keras.models import Sequential
    from keras.layers import LSTM, Dense

    model = Sequential()
    model.add(LSTM(L1, return_sequences=True, input_shape=(8, 20)))
    model.add(LSTM(L2, return_sequences=False))
    model.add(Dense(units=1))
    model.compile(loss='mean_squared_error', optimizer='adam', metrics=['mae'])

    return model

def Perceptron(L1=160):

    import os
    os.environ['KERAS_BACKEND'] = 'tensorflow'
    from keras.models import Sequential
    from keras.layers import Dense, Flatten

    model = Sequential()
    model.add(Flatten(input_shape=(8, 20)))
    model.add(Dense(L1))
    model.add(Dense(L1))
    model.add(Dense(units=1))
    model.compile(loss='mean_squared_error', optimizer='adam', metrics=['mae'])
    return model

def CNN(L1=160):
    import os
    os.environ['KERAS_BACKEND'] = 'tensorflow'
    from keras.models import Sequential
    from keras.layers import Conv2D, Flatten, Dense

    model = Sequential()
    # Add a Conv2D layer with filters that convolve along the height (8) while preserving the width (20)
    model.add(Conv2D(L1, (8, 1), activation='relu', input_shape=(8, 20, 1)))
    # Flatten the output to feed into Dense layers
    model.add(Flatten())
    # Output layer with a single unit for continuous output
    model.add(Dense(1))
    model.compile(loss='mean_squared_error', optimizer='adam', metrics=['mae'])
    return model

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

# Paul Tol's Notes color scheme
colors = ['#4477AA','#66CCEE','#228833','#CCBB44','#EE6677','#AA3377','#BBBBBB']

# Figure specifications
# Set font type to Arial
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 16

# Parse command-line arguments
parser = argparse.ArgumentParser(
    description="Compare different model architectures for sequence-to-function prediction."
)

parser.add_argument(
    "-i", "--input",
    default="data/Excel/third_round/enrichment_score_with_amino_acid_sequence_threshold.xlsx",
    help="Input directory containing sample folders"
)
parser.add_argument(
    "-o", "--output",
    default="plots",
    help="Output directory containing ML model performance plot"
)

# Load the Excel file into a DataFrame
args = parser.parse_args()
file_path = args.input
df = pd.read_excel(file_path)

# Create output directory if needed
output_dir = args.output
os.makedirs(output_dir, exist_ok=True)

# log2 transform the 'averaged_enrichment_score' column
df['log2_averaged_enrichment_score'] = df['averaged_enrichment_score'].apply(lambda x: np.log2(x))

# Define the models to be used
models = {
    'Round_3_LR_Threshold': LR,
    'Round_3_FNN_Threshold': Perceptron,
    'Round_3_LSTM_Threshold': LongShortTermMemoryModel,
    'Round_3_CNN_Threshold': CNN
}

# Define number of rows to leave out for testing
shuffle_split = ShuffleSplit(n_splits=5, test_size=0.1, random_state=0)

# Define percentages
percentages = [100, 90, 80, 70, 60, 50, 40, 30, 20, 10]

# Initialize dictionaries to store results
pearson_results = {model_name: {f'test_set_{percentage}': [] for percentage in percentages} for model_name in models.keys()}

for model_name, model_func in models.items():

    print(f'Model: {model_name}')

    for train_index, test_index in shuffle_split.split(df):

        print((len(train_index), len(test_index)))

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

        for key, test_set in test_sets.items():

            print(f'Test set: {key}')

            # do one hot encoding of amino acids
            train_x = np.asarray([AA_hotencoding(variant) for variant in train["amino_acid_sequence"]])
            train_y = np.asarray([score for score in train['log2_averaged_enrichment_score']])
            test_x = np.asarray([AA_hotencoding(variant) for variant in test_set["amino_acid_sequence"]])
            test_y = np.asarray([score for score in test_set['log2_averaged_enrichment_score']])

            if model_name == 'Round_3_LR_Threshold':

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

            else: 

                # Create the model
                model = model_func()

                # Train the model
                model.fit(
                    train_x,
                    train_y,
                    epochs=50,
                    batch_size=500,
                    verbose=0
                )

                # Store the results for the current model
                y_pred = model.predict(test_x)
                y_pred = np.reshape(y_pred, (1,y_pred.shape[0]))[0]
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
        key = f'test_set_{percentage}'
        for value in pearson_results[model_name][key]:
            plot_data.append({'Model': model_name, 'Test Set Percentage': percentage, 'Pearson Correlation': value})

plot_df = pd.DataFrame(plot_data)

# Plot the results using seaborn
plt.figure(figsize=(8.25, 6))
sns.lineplot(data=plot_df, x='Test Set Percentage', y='Pearson Correlation', hue='Model', palette=colors, marker=None, errorbar=('ci', 95))

plt.xlabel('Top X% of Test Set', fontsize=16)
plt.ylabel("5-Fold CV test Pearson's correlation", fontsize=16)
# invert the x-axis so it plots from larger to smaller
plt.gca().invert_xaxis()

# legend above plot centered
plt.legend(fontsize=16, loc='upper center', bbox_to_anchor=(0.5, 1.30), ncol=2, frameon=False)
plt.gca().tick_params(labelsize=16)
plt.ylim(0, 1)
plt.tight_layout()

plt.savefig(f"{output_dir}/figure_2d.png")
