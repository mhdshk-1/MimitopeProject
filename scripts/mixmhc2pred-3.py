# This script takes the csv output from mixmhc2pred-2.py in the script-created folder
# and combines the %Rank columns into a single column for each allele, categorising
# the peptides into Strong Binders (<=2% Rank) or not (>2% Rank). The resulting
# single csv file is saved to the same folder and can then be submitted to 
# sequence-alignment.py for further analysis.

import pandas as pd
import os

# set working directory to the output folder of mixmhc2pred-2.py
working_directory = "C:\\Users\\test\\Documents\\GitHub\\MimitopeDiscovery\\Muhammad\\data\\t1dm"
os.chdir(working_directory)

# Tell the script where the mixMHCPred input csv is, for indexing of the peptides against their SeqIDs
csv_file = "C:\\Users\\test\\Documents\\GitHub\\MimitopeDiscovery\\Muhammad\\data\\mixmhc2pred-input.csv"

# Get the list of files in the directory
files = os.listdir(working_directory)
# Sort the files alphanumerically
sorted_files = sorted(files)
# Select the first file
first_file = sorted_files[0]
# Construct the full path to the first file
file_path = os.path.join(f"{working_directory}/{first_file}")
# open the output file for the allele
with open(file_path, "r") as file:
    lines = file.readlines()[19:]
# Convert the subsequent lines into a TSV format
tsv = [line.strip().split("\t") for line in lines]
# define new tsv keeping only column 0
tsv = [[line[0]] for line in tsv]
# Convert new_tsv list to DataFrame
df = pd.DataFrame(tsv[1:], columns=tsv[0])

# Iterate over the sorted files
for file in sorted_files:
    # Construct the full path to the file
    file_path = os.path.join(f"{working_directory}/{file}")
    
    # Open the file and read the lines
    with open(file_path, "r") as file:
        lines = file.readlines()[19:]
    
    # Convert the subsequent lines into a TSV format
    tsv = [line.strip().split("\t") for line in lines]
    
    # Extract column 7 and set it as a new column in the DataFrame
    column_name = tsv[0][7]
    column_data = [line[7] for line in tsv[1:]]
    df[column_name] = column_data

# Append df to mixmhc2pred-input.csv
csv_df = pd.read_csv(csv_file)
df = pd.merge(df, csv_df, left_on='Peptide', right_on='Protein')

# Remove rows with NA values in column 1
df = df[df.iloc[:, 1] != "NA"]

# Drop the "Protein" and "ContextSequence" columns
df = df.drop(columns=['Protein', 'ContextSequence'])

# Make the last two columns (SeqID, Positions) the first two columns
columns = df.columns.tolist()
df = df.reindex(columns=columns[-2:] + columns[:-2])

# Get the column names starting with '%' (i.e. all %Rank columns)
columns_to_scan = [col for col in df.columns if col.startswith('%')]

# Iterate over the selected columns 
for col in columns_to_scan:
    df[col] = df[col].astype(float)  # Convert column values to float

# Reset the index after removing rows
df = df.reset_index(drop=True)

# Create a 4th column indicating whether the peptide is a strong binder
# Strong binders defined in our code as being in the top two centiles 
# The model ranks peptides by comparing affinities to a random distribution of natural peptides 
# 10.1038/s43018-023-00548-5 chose =<2% as strong binders
# 10.1038/s41588-022-01273-y chose =<1% as strong binders
df.insert(3, 'StrongBinder', df.iloc[:, 3:].apply(lambda x: any(val <= 2 for val in x), axis=1))

# Display the resulting DataFrame
print(df)
print(df.dtypes)
df.to_csv("viral_epitopes.csv", index=False)
print("Saved to working directory.")