#!/usr/bin/env python

# coding: utf-8 

# Written by ChatGPT 3.5

def read_fasta(file_path):
    sequences = {}
    current_seq = ""
    with open(file_path, "r") as fasta_file:
        for line in fasta_file:
            line = line.strip()
            if line.startswith(">"):
                current_seq = line[1:]
                sequences[current_seq] = ""
            else:
                sequences[current_seq] += line
    return sequences

def remove_columns_with_gaps(sequences):
    alignment_length = len(list(sequences.values())[0])
    columns_to_remove = set()

    for i in range(alignment_length):
        column = [seq[i] for seq in sequences.values()]
        if "-" in column:
            columns_to_remove.add(i)

    filtered_sequences = {
        name: "".join(seq[i] for i in range(alignment_length) if i not in columns_to_remove)
        for name, seq in sequences.items()
    }

    return filtered_sequences

def write_fasta(output_path, sequences):
    with open(output_path, "w") as output_file:
        for name, seq in sequences.items():
            output_file.write(f">{name}\n{seq}\n")

input_fasta = "Combined_CDS_selected_ungapped.fasta"
output_fasta = "Combined_CDS_selected_ungapped.fasta"

sequences = read_fasta(input_fasta)
filtered_sequences = remove_columns_with_gaps(sequences)
write_fasta(output_fasta, filtered_sequences)

print(f"Columns with gaps removed and saved to {output_fasta}.")