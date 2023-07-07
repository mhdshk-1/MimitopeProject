# This script reads your fasta viral proteome and chunks it for
# submission to the NetChop server. Each submission should be
# at most 100 sequences and 100,000 AAs per submission.
# Netchop does not run locally on WSL, and using urllib to submit
# a HTML POST request to both the Netchop DTU and IEDB servers
# does not work due to a two-step submission process, so this
# is the only workaround I could effectuate.

# Import the necessary modules
from Bio import SeqIO
import os

# Set your working directory path
working_dir = "C:\\Users\\test\\Documents\\GitHub\\MimitopeDiscovery\\Muhammad\\data"
viral_proteome = "viral_proteome.fasta"
output_dir = "chunked_fastas_with_html"

# Set working directory to folder containing script
os.chdir(working_dir)

# Define a function to chunk our fasta file
def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.
    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.Align.parse(...), or simply
    lines from a file handle.
    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    batch = []
    for entry in iterator:
        batch.append(entry)
        if len(batch) == batch_size:
            yield batch
            batch = []
    if batch:
        yield batch

# Open the fasta file and chunk it into smaller files using
# the batch_iterator function. ensure fasta is in script directory
record_iter = SeqIO.parse(open(viral_proteome), "fasta")
for i, batch in enumerate(batch_iterator(record_iter, 100)):
    filename = os.path.join(output_dir, "group_%i.fasta" % (i + 1))
    with open(filename, "w") as handle:
        count = SeqIO.write(batch, handle, "fasta")
    print("Wrote %i records to %s" % (count, filename))

# The resulting files should now be small enough to submit to the webserver.
# Use C-term 3.0 as the prediction method and 0.5 as the threshold.
# Save the results to html with the same filename (eg group_1.fasta becomes group_1.html)
# in the same directory. Then, run netchop-2.py to convert all the html files to csv files.