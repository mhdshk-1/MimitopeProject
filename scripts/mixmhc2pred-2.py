# This script takes the appropriately formatted csv from mixmhc2pred-2.py,
# converts it to a txt file whic is then saved in the same directory,
# and runs mixmhc2pred.exe on it via command line. Set your alleles of interest
# in the alleles list using the defintions from Alleles_list_Human.txt in the
# PMWdef folder in the exe directory. Having mixMHCpred-2.0.exe installed is a 
# prerequisite for this script to work. Follow the instructions on 
# https://github.com/GfellerLab/MixMHC2pred and ensure the programme directory
# is in PATH. The output files are saved to a subfolder called mixmhc2pred_output.

import pandas as pd
import subprocess
import os

# Set working directory to folder containing all the csv input file formatted by mixmhc2pred-2.py
# N.b. there are two working directory variables in this script, differently formatted for command line and Python
# Set input_file to the appropriate csv file
# Double backslashes as single backslash is an escape character
wd = "C:\\Users\\test\\Documents\\GitHub\\MimitopeDiscovery\\Muhammad\\data" # This is for Python
working_directory = "C:/Users/test/Documents/GitHub/MimitopeDiscovery/Muhammad/data" # This is for command line
input_file = "mixmhc2pred-input.csv"
output_folder_name = "nt1" # All of the output files will be saved in this subfolder in the working directory
# set executable path to folder containing mixmhc2pred.exe
executable_path = "C:/Users/test/MixMHC2pred-2.0" # This is for command line
# Set alleles of interest, this is at line 55 in this script

# Set working directory
os.chdir(wd)

# read in the csv file
df = pd.read_csv(input_file)
# Print the dataframe
print(df)
print(df.dtypes)
# Drop the first, third, and fourth columns (context sequence is not required for pre-cleaved peptides)
df = df.drop(df.columns[[0, 2, 3]], axis=1)
# save df as a tsv but with txt extension so mixmhc2pred will accept it
df.to_csv('mixmhc2pred-input.txt', sep='\t', index=False, header=False)

# Run mixmhc2pred via command line
# The code below builds the allele list from all >6000 allele definitions in Alleles_list_Human.txt
# This is not necessary if you only want to run a few alleles, in which case you can just set the alleles list manually
# On my computer running all the alleles took >12hrs

"""
# now we want to grab the allele definitions from mixmhc2pred's Alleles_list_Human.txt file
with open(f"{executable_path}/PWMdef/Alleles_list_Human.txt", "r") as file:
    lines = file.readlines()[4:]
# Convert the subsequent lines into a TSV format
tsv_lines = [line.strip().split("\t") for line in lines]
# Store the values of the first column as a list
alleles = [line[0] for line in tsv_lines]
"""
# Set the alleles list manually; refer to Alleles_list_Human.txt in the PMWdef folder in the exe directory
# Each as a text string

"""
alleles = ["DQA1_05_01__DQB1_02_01", # for T1DM, top 4 lines risk, bottom 4 lines protective, see worksheet.docx in docs for citations
           "DQA1_03_01__DQB1_03_02",
           "DRB1_04_01",
           "DRB1_04_09",

           "DRB1_13_03",
           "DRB1_15_01",
           "DQA1_01_01__DQB1_05_03",
           "DQA1_02_01__DQB1_03_03"] 
           """
"""
alleles = ["DRB1_15_01", # for MS, top 6 lines risk, bottom 6 lines protective
           "DRB1_03_01",
           "DRB1_13_03",
           "DRB1_04_01",
           "DRB1_04_09",
           "DRB1_14_01",

           "DRB1_01_01",
           "DRB1_11_01",
           "DRB1_11_02",
           "DRB1_11_03",
           "DRB1_11_04",
           "DRB1_11_05"] 
           """

alleles = ["DQA1_01_01__DQB1_06_02", # for NT1, top line risk, bottom line protective
           
           "DQA1_01_01__DQB1_05_01"]

# Create a subfolder to save the output files
output_folder = os.path.join(working_directory, output_folder_name)
os.makedirs(output_folder, exist_ok=True)

# Iterate over each allele
for allele in alleles:
    # Generate the output file path for the current allele
    output_file = os.path.join(output_folder, f"output_{allele}.txt")
        # Construct the command to run the executable
    command = [f"{executable_path}/MixMHC2pred.exe", "-i", f"{working_directory}/mixmhc2pred-input.txt",
               "-o", output_file, "-a", allele, "--no_context"]
        # Run the command
    subprocess.run(command, shell=True)
