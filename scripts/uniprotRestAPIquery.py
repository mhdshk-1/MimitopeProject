import requests
import csv
import os

# set working directory to same as blast.py
working_directory = "C:/Users/test/Documents/GitHub/MimitopeDiscovery/Muhammad"
os.chdir(working_directory)

# When you obtain a proteome from the Human Protein Atlas, it is in the form of a TSV file.
# This code snippet uses the TSV SeqIDs to download the FASTA file for the proteome from UniProt.
tsv_file = "humanproteinatlas_pancreas_proteome.tsv" # Filename for the input TSV file
fasta_filename = 'pancreas_proteome.fasta' # Filename for the output FASTA file
uniprot_ids = [] # List to store the UniProt IDs

# Extract the UniProt IDs from 5th column in the TSV file and add them to the uniprot_ids list
with open(tsv_file, "r") as file:
    reader = csv.reader(file, delimiter="\t")
    next(reader)  # Skip the header row
    for row in reader:
        uniprot_ids.append(row[4])

print(uniprot_ids)

# For each ID in the list, download the FASTA file and append it to the output file
with open(fasta_filename, 'w') as fasta_file:
    for uniprot_id in uniprot_ids:
        url = f'https://www.uniprot.org/uniprot/{uniprot_id}.fasta'
        response = requests.get(url)
        fasta_content = response.text

        # Write the FASTA content to the output file
        fasta_file.write(fasta_content)

print("FASTA file with all sequences saved.")

"""
# This code can be used to build the url for the query to view the fasta file in your browser and save it
# N.b. there's an upper limit on the number of Seqs you can query at once
uniprot_ids_str = '%2C'.join(uniprot_ids)
url = f"https://rest.uniprot.org/uniprotkb/accessions?accessions={uniprot_ids_str}&compressed=false&format=fasta"
print(url)
"""