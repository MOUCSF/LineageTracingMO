import sys
from Bio import pairwise2
import pandas as pd

def create_state_vector(alignment, original_seq):
    aligned_seq = alignment[1]
    state_vector = []
    for o, a in zip(original_seq, aligned_seq):
        if o == a:
            state_vector.append(0)
        else:
            state_vector.append(1)
    return state_vector

def main(input_file, output_file):
    # original sequence
    original_seq = "TTTGCTACCTATTACTAGGACAAGTGGAGGCGTTGAGAATGTCTCGGAATGTGCCGGAAATTTGCATCGTATTACTAGGACAATTCGAGTCCTTGAGGAAGTCTCTCGATGTGCCGGAAA"
    
    # input file reading
    data = []

    with open(input_file, 'r') as file:
        header = file.readline().strip().split()
        for line in file:
            data.append(line.strip().split())

    # Process each mutate sequence and generate the state matrix
    output_data = []

    for row in data:
        tenx = row[0]
        static = row[1]
        mutate_seq = row[2]

        # Align mutate_seq with original_seq using Needleman-Wunsch algorithm
        alignments = pairwise2.align.globalms(original_seq, mutate_seq, 2, -1, -0.5, -0.1)

        # Choose the best alignment (first one in the list)
        best_alignment = alignments[0]

        # Create the state vector based on the best alignment
        state_vector = create_state_vector(best_alignment, original_seq)

        # Add row to output data
        output_data.append([tenx, static] + state_vector)

    # Create DataFrame for output
    columns = ['tenx', 'static'] + [f'V{i+1}' for i in range(len(state_vector))]
    output_df = pd.DataFrame(output_data, columns=columns)

    # Write output to file
    output_df.to_csv(output_file, index=False, sep=',')

    print(f"Output written to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python daisyalignment.py input.txt output.txt")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    main(input_file, output_file)