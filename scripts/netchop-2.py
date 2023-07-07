# This script converts all html netchop results to csv files
# in a given directory. The csv files are saved in the same directory.
# These files are then read back in and joined together into a single
# simpler csv file for further analysis. Ensure 'html_folder' is set
# to the correct directory. Backslashes are escaped with another backslash

import os
import csv
from bs4 import BeautifulSoup
import pandas as pd
import numpy as np

# Set the working directory to the folder containing HTML files
html_folder = "C:\\Users\\test\\Documents\\GitHub\\MimitopeDiscovery\\Muhammad\\data\\chunked_fastas_with_html"
# Name the resulting dataframe that will be created from combining all csv files
output_file = "viral_proteome_cleaved.csv"

# Set the working directory
os.chdir(html_folder)

# Define a function to extract text from between HTML <pre> tags
def extract_text_between_delimiters(html_file):
    with open(html_file, "r") as file:
        html_content = file.read()

    soup = BeautifulSoup(html_content, "html.parser")
    pre_tags = soup.find_all("pre")
    extracted_text = []
    for tag in pre_tags:
        text = tag.get_text()
        extracted_text.append(text)

    return extracted_text

# Get a list of HTML files in the directory
html_files = [file for file in os.listdir() if file.endswith(".html")]

# Process each HTML file and save the corresponding CSV file
for html_file in html_files:
    # Extract the text from the HTML file
    extracted_text = extract_text_between_delimiters(html_file)

    # Split the text by newline character to get a list of lines
    lines = str(extracted_text).split("\\n")

    # Create a DataFrame from the list of lines, specifying the column positions
    colspecs = 0, 4, 8, 14, 23, 34
    df = pd.DataFrame(
        [[row[i:k] for i, k in zip(colspecs[:-1], colspecs[1:])] for row in lines[1:]],
        columns=["position", "AA", "cleaved", "score", "SeqID"]
    )

    # Convert the position column to numeric type, coercing non-numeric values to NaN
    df["position"] = pd.to_numeric(df["position"], errors="coerce")

    # Use boolean indexing to remove rows with NaN values in the position column
    df_filtered = df[pd.notnull(df["position"])]

    # Remove space characters from all cells in the DataFrame
    df_filtered["AA"] = df_filtered["AA"].str.strip()
    df_filtered["cleaved"] = df_filtered["cleaved"].str.strip()

    # Specify the column data types
    df_filtered = df_filtered.astype(
        {"position": int, "AA": str, "cleaved": str, "score": float, "SeqID": str}
    )

    # Change the data type of the cleaved column to boolean
    df_filtered["cleaved"] = np.where(df_filtered["cleaved"] == "S", True, False)

    # Save the DataFrame as a CSV file with the same name as the HTML file
    csv_file = os.path.splitext(html_file)[0] + ".csv"
    df_filtered.to_csv(csv_file, index=False)

    print(f"Processed: {html_file} -> Saved as: {csv_file}")

# Now we read the csv files back in and join them together into a single dataframe containing SeqIDs, cleavage produces, and positions
# Create a dataframe incorporating all csv files in the folder
def merge_csv_files():
    # Get the current directory
    directory = os.getcwd()

    # Get a list of CSV files in the directory
    csv_files = [file for file in os.listdir(directory) if file.endswith(".csv")]

    # Create an empty list to store DataFrames
    dfs = []

    # Read each CSV file and append its DataFrame to the list
    for file in csv_files:
        file_path = os.path.join(directory, file)
        df = pd.read_csv(file_path)
        dfs.append(df)

    # Merge the DataFrames into one
    df = pd.concat(dfs, ignore_index=True)

    return df
df = merge_csv_files()

# double check that the dataframe is correct
print(df)
print(df.dtypes)

# Define a function to return, for each SeqID, tuples of the AA sequences and their positions
def separate_protein(df):
    protein = []
    positions = []
    current_seq_id = None
    
    # Iterate over the rows of the DataFrame
    for _, row in df.iterrows():
        # Check if the current SeqID has changed
        if current_seq_id != row['SeqID']:
            # If it has changed, yield the previous protein and positions (if any)
            if protein:
                yield current_seq_id, "".join(protein), positions
            
            # Reset the protein and positions lists
            protein = []
            positions = []
            current_seq_id = row['SeqID']
        
        # Add the AA to the protein list
        protein.append(row['AA'])
        positions.append(row['position'])
        
        # Check if the cleaved column is True
        if row['cleaved']:
            # Yield the protein and positions as tuples
            yield current_seq_id, "".join(protein), positions
            
            # Reset the protein and positions lists
            protein = []
            positions = []
    
    # If there are remaining amino acids after the last cleaved row, yield them as well
    if protein:
        yield current_seq_id, "".join(protein), positions

# Call the iterator function and pass the DataFrame
protein_iterator = separate_protein(df)

# Create a list to store the protein data
protein_data = []

# Iterate over the separated proteins and store them in the list
for seq_id, protein, positions in protein_iterator:
    protein_data.append([seq_id, protein, positions])

# Write the protein data to the CSV file
with open(output_file, mode="w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(["SeqID", "Protein", "Positions"])  # Write the header
    writer.writerows(protein_data)  # Write the protein data rows

print("Protein data has been written to", output_file)