# This script takes the output from netchop-2.py and creates a csv file that can be used as input for mixMHC2pred-2.py
# which then converts this to txt and submits it to mixMHC2pred via command line. Strictly speaking, these could both 
# be combined into one script, but I wanted to keep them separate for clarity. Lines 41-82 extract the three amino acids
# either side of the cleavage point of the cleaved peptides by accessing the original fasta file (i.e. the context sequence,
# see https://github.com/GfellerLab/MixMHC2pred for the readme which explains this). # However the script
# bugs out when the cleaved peptide occurs within 3 AAs of the start or end of the protein, as the script tries to 
# access a negative index/inappropriately high index because it reads 3 AAs at a time. This could be fixed by reading 
# 1 AA at a time and inserting an X when out of bounds. I have not done this as it is not necessary for the current
# method where context sequences are not relevant in the setting of cleaved peptides. I've kept the buggy code in 
# just in case we decide context sequences are relevant, in which case it can be relatively easily modified. 

from Bio import SeqIO
import pandas as pd
import os

# set working directory to folder containing all csv outputs from netchop-2.py
os.chdir("C:\\Users\\test\\Documents\\GitHub\\MimitopeDiscovery\\Muhammad\\data")

# read in netchop-2.py output
df = pd.read_csv("chunked_fastas_with_html\\viral_proteome_cleaved.csv")

# tell the script where the original fasta file is
fasta_file = "viral_proteome.fasta"

# delete the first three characters of the SeqID column from the netchop-2.py output csv
df.iloc[:, 0] = df.iloc[:, 0].str.slice(3)

# remove any underscores from the SeqID column (these two steps give us a usable accession code)
df.iloc[:, 0] = df.iloc[:, 0].str.replace("_", "")

# Remove any rows where protein <9 peptides (this is the binding core length for MHCII, see 10.1093/bioinformatics/btl479)
# N.b. mixMHC2pred's input contstraints: 12 < peptide =< 21 AAs
df = df[df['Protein'].str.len() >= 9]

#create a column to store the context sequence, filling it with the first three and last three amino acids of the peptide
df['ContextSequence'] = df['Protein'].str[:3] + df['Protein'].str[-3:]

# check that the dataframe is correct
print(df)
print(df.dtypes)

# Iterate over the rows of the CSV
for index, row in df.iterrows():
    # Get the SeqID from the 'SeqID' column
    seq_id = row['SeqID']
    
    # Get the positions from the 'Positions' column
    positions = eval(row['Positions'])
    
    # Extract the three amino acids before the first position by parsing the fasta file
    start_position = positions[0] - 4
    end_position = positions[0] - 1
    extracted_sequence = ""
    for record in SeqIO.parse(fasta_file, "fasta"):
        if seq_id in record.id:
            extracted_sequence = str(record.seq[start_position:end_position])
            break
    # append them to the start of the context sequence
    df.at[index, 'ContextSequence'] = extracted_sequence + row['ContextSequence']

# Print the updated dataframe
print(df)
print(df.dtypes)

# repeat the process for the last three amino acids
# Iterate over the rows of the CSV
for index, row in df.iterrows():
    # Get the SeqID from the 'SeqID' column
    seq_id = row['SeqID']
    
    # Get the positions from the 'Positions' column
    positions = eval(row['Positions'])

    # Extract the three amino acids after the last position by parsing the fasta file
    start_position = positions[-1]
    end_position = positions[-1] + 3
    extracted_sequence = ""
    for record in SeqIO.parse(fasta_file, "fasta"):
        if seq_id in record.id:
            extracted_sequence = str(record.seq[start_position:end_position])
            break
    # append them to the end of the context sequence
    df.at[index, 'ContextSequence'] = row['ContextSequence'] + extracted_sequence

# Print the updated dataframe
print(df)
print(df.dtypes)

# save as a csv
df.to_csv("mixmhc2pred-input.csv", index=False)
print("Done!")