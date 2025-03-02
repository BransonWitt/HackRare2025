from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import random

# Read in the input FASTA file
input_file = "C:\\Users\\brans\\Downloads\\Harvard_Hackathon\\GCF_000001405.13_GRCh37_genomic.fna"
output_file = "C:\\Users\\brans\\Downloads\\Harvard_Hackathon\\modified_proteins.fna"

#function to artificially mutate specific genes
def artificialMutation(sequence, mutation_rate=0.04):
    print(sequence)
    
    seq_list = list(sequence)
    
    i = 0
    while i < len(seq_list):
        if seq_list[i] in ['A', 'T', 'C', 'G']:
            rand_val = random.random()
            
            if rand_val < mutation_rate:  # Point mutation
                current_base = seq_list[i]
                possible_bases = ['A', 'T', 'C', 'G']
                possible_bases.remove(current_base)
                new_base = random.choice(possible_bases)
                seq_list[i] = new_base

        i += 1
    print(''.join(seq_list), '\n')
    
    return ''.join(seq_list)

# List of sequence IDs to translate (modify this as needed)
selected_ids = ['NC_000016.9', 'NC_000009.11', 'NC_000008.10']  # Replace with the sequence IDs you want to translate

# Function to replace a segment of the DNA sequence
def replace_dna_segment(dna_sequence, start, end):
    """
    Replaces a segment from start to end in the DNA sequence with a replacement sequence.
    """
    # Replace the segment in the sequence
    return dna_sequence[:start] + artificialMutation(dna_sequence[start:end]) + dna_sequence[end:]

# Specify the segments to replace for each sequence
segments_to_replace = pd.read_csv("C:\\Users\\brans\\Downloads\\ncbi_dataset\\ncbi_dataset\\data\\GCF_000001405.13\\epilepsy.csv", index_col="Index")

# Open the output file
with open(output_file, 'w') as output_fasta:
    # Parse the FASTA file
    for record in SeqIO.parse(input_file, 'fasta'):
        # Check if the current sequence needs a replacement
        if record.id in selected_ids:
            specified_df = segments_to_replace[segments_to_replace['Gene'] == str(record.id)]
            
            for index, row in specified_df.iterrows():
                # Replace the segment in the DNA sequence
                record.seq = replace_dna_segment(record.seq, row["Start"], row["End"])
        
        # Write the modified or original record to the output file
        SeqIO.write(record, output_fasta, 'fasta')

print(f"Modified sequences saved to {output_file}")