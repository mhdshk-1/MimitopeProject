import re
import os

os.chdir("C:\\Users\\test\\Documents\\GitHub\\MimitopeDiscovery\\Muhammad\\tests")

# Read the input text file
input_file = "pancreas_proteome.fasta"
with open(input_file, "r") as file:
    content = file.read()

# Define the regex pattern to match the tags and everything in between
pattern = r"<!doctype html>.*?</html>"

# Remove the matched pattern from the content using regex substitution
clean_content = re.sub(pattern, "", content, flags=re.DOTALL)

# Write the cleaned content to an output file
output_file = "pancreas_proteome_cleaned.fasta"
with open(output_file, "w") as file:
    file.write(clean_content)

print("done")
