import pandas as pd
from Bio import SeqIO
from Bio import Align

# Read the dataframe and extract the query sequences
queries_path = "C:\\Users\\test\\Documents\\GitHub\\MimitopeDiscovery\\Muhammad\\data\\t1dm\\viral_epitopes.csv"
queries = pd.read_csv(queries_path)
queries = queries[queries.iloc[:, 3]]  # Keep rows where the fourth column is True
query_sequences = queries.iloc[:, 2].tolist()  # Builds peptide list from 3rd column
print("query_sequences: ", query_sequences)

# Read the FASTA file to be searched against
fasta_file = "C:\\Users\\test\\Documents\\GitHub\\MimitopeDiscovery\\Muhammad\\data\\t1dm\\control_proteome_NOT_pancreas.fasta"
fasta_sequences = list(SeqIO.parse(fasta_file, 'fasta')) # The fasta_sequences iterator is converted to a list, 
    # ensuring that the sequences can be iterated over multiple times in the subsequent loops
    # otherwise the iterator will be exhausted after the first loop

# Create a pairwise aligner object
aligner = Align.PairwiseAligner()

"""
This object stores the match and mismatch scores, as well as the gap scores.
Typically, match scores are positive, while mismatch scores and gap scores 
are negative or zero. By default, the match score is 1 (I have used 2 to 
weight it more heavily), and the mismatch and gap scores are zero. Based on 
the values of the gap scores, a PairwiseAligner object automatically chooses
the appropriate alignment algorithm (the Needleman-Wunsch, Smith-Waterman, 
Gotoh, or Waterman-Smith-Beyer global or local alignment algorithm). Read 
the documentation at: https://biopython.org/docs/1.76/api/Bio.Align.html
"""

# Set the parameters for local alignment
aligner.mode = 'local'
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -1
aligner.extend_gap_score = -1
window = 5   # The length of the AA match region

# Create an empty list to store the match information
matches = []

# Iterate over each query sequence
for query in query_sequences:
    # Generate all possible {window} AA subsequences from the query
    subsequences = []
    for i in range(len(query) - (window - 1)):
        subsequence = query[i:i+window]
        subsequences.append(subsequence)
        print("all possible", window, "AA subsequences for", query, "are", subsequences)

    # Create an empty list to store the match information for the current query
    current_matches = []

    # Iterate over each sequence in the fasta file
    for fasta_seq in fasta_sequences:
        sequence = str(fasta_seq.seq)
        seq_len = len(sequence)

        # Iterate over each possible {window} AA region in the sequence
        for i in range(seq_len - (window - 1)):
            subsequence = sequence[i:i+window]

            # Perform local alignment for each {window} AA subsequence of the query
            for subseq_query in subsequences:
                alignments = aligner.align(subsequence, subseq_query)

                # Iterate over each alignment
                for alignment in alignments:
                    if alignment.score == 2 * len(subseq_query):
                        start_pos = i + 1
                        end_pos = i + window
                        match_info = {
                            'SeqID': fasta_seq.id,
                            'Query': query,
                            'Subsequence': subseq_query,
                            'Hit': alignment.target, # This doesn't really need to be printed but I just did it for validation
                            'Start': start_pos,
                            'End': end_pos,
                        }
                        current_matches.append(match_info)

    # Add the matches for the current query to the overall matches list
    matches.extend(current_matches)
    print("current matches:", current_matches)

# Create a DataFrame from the matches list
df = pd.DataFrame(matches)

# Print the DataFrame
print(df)

# Save the DataFrame to a CSV file
df.to_csv("C:\\Users\\test\\Documents\\GitHub\\MimitopeDiscovery\\Muhammad\\data\\t1dm\\control_alignment.csv", index=False)
print("Done")