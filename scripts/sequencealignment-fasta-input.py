import pandas as pd
from Bio import SeqIO
from Bio import Align

# Read the query sequences from the fasta file
query_fasta_file = "C:\\Users\\test\\Documents\\GitHub\\MimitopeDiscovery\\Muhammad\\data\\viral_proteome.fasta"
query_sequences = []
query_seq_ids = []
for record in SeqIO.parse(query_fasta_file, "fasta"):
    query_sequences.append(str(record.seq))
    query_seq_ids.append(record.id)

print("query_sequences:", query_sequences)
print("query_seq_ids:", query_seq_ids)

# Read the FASTA file to be searched against
fasta_file = "C:\\Users\\test\\Documents\\GitHub\\MimitopeDiscovery\\Muhammad\\data\\t1dm\\control_proteome_NOT_pancreas.fasta"
fasta_sequences = list(SeqIO.parse(fasta_file, 'fasta'))  # Convert fasta_sequences iterator to a list

# Create a pairwise aligner object
aligner = Align.PairwiseAligner()

# Set the parameters for local alignment
aligner.mode = 'local'
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -1
aligner.extend_gap_score = -1
window = 6  # The length of the AA match region

# Create an empty list to store the match information
matches = []

# Iterate over each query sequence
for query, query_seq_id in zip(query_sequences, query_seq_ids):
    # Generate all possible {window} AA subsequences from the query
    subsequences = [query[i:i+window] for i in range(len(query) - (window - 1))]
    print("all possible", window, "AA subsequences for", query, "are", subsequences)

    # Create an empty list to store the match information for the current query
    current_matches = []

    # Iterate over each sequence in the fasta file
    for fasta_seq in fasta_sequences:
        sequence = str(fasta_seq.seq)
        seq_len = len(sequence)

        # Iterate over each possible {window}-AA region in the sequence
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
                            'Query_SeqID': query_seq_id,
                            'Query': query,
                            'Subsequence': subseq_query,
                            'Hit': alignment.target,
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
df.to_csv("C:\\Users\\test\\Documents\\GitHub\\MimitopeDiscovery\\Muhammad\\data\\t1dm\\control_direct_alignment.csv", index=False)
print("Done")
