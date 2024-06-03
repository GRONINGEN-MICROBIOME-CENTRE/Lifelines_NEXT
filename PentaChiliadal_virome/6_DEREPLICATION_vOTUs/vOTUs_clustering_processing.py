import argparse

# Define command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('Clusters', help='input file with clusters of viral sequences')
parser.add_argument('Negcontrol_sequences', help='input file with predicted viral sequences in the negative control sample')
args = parser.parse_args()

#Read the files
with open(args.Clusters, 'r') as clusters, open(args.Negcontrol_sequences, 'r') as negcontrol, open('Final_viral_sequences_after_dereplication.txt', 'w') as output:
    # Read the sequence names to exclude from negative control file into a set
    exclude = set(line.strip() for line in negcontrol)

# Loop over the lines in clusters file:
    for line in clusters:
        # Split the line by tabs and discard the first column (representative sequence - it is already included in the 2nd column)
        cols = line.strip().split('\t')[1:]
        # Split each sequence by comma
        seqs = [seq.strip() for col in cols for seq in col.split(',')]
        # Check if any of the sequences in the line are in the exclude set
        if any(seq in exclude for seq in seqs):
            # If so, skip the line
            continue
        # If no sequence is present in the exclude set, write the sequences to output file
        for seq in seqs:
            output.write(seq + '\n')
