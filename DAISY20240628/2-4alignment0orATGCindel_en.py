from Bio.Align import PairwiseAligner
import pandas as pd
import sys

def create_state_vector(alignment, original_seq_length):
    aligned_original_seq = alignment[0]
    aligned_mutate_seq = alignment[1]
    state_vector = []
    for o, a in zip(aligned_original_seq, aligned_mutate_seq):
        if o == a:
            state_vector.append('0')
        elif o == '-' or a == '-':
            state_vector.append('-')
        else:
            state_vector.append(a)
    
    # Ensure the state vector is the same length as the original sequence
    # Fill with '0' if the state vector is shorter
    if len(state_vector) < original_seq_length:
        state_vector += ['0'] * (original_seq_length - len(state_vector))
    # Trim if it's longer (this case shouldn't normally happen in global alignment)
    elif len(state_vector) > original_seq_length:
        state_vector = state_vector[:original_seq_length]
    
    return state_vector

def main(input_file, output_file):
    # original sequence
    original_seq = "TTTGCTACCTATTACTAGGACAAGTGGAGGCGTTGAGAATGTCTCGGAATGTGCCGGAAATTTGCATCGTATTACTAGGACAATTCGAGTCCTTGAGGAAGTCTCTCGATGTGCCGGAAA"
    original_seq_length = len(original_seq)
    
    # input file reading
    data = []

    with open(input_file, 'r') as file:
        header = file.readline().strip().split()
        for line in file:
            data.append(line.strip().split())

    # Process each mutate sequence and generate the state matrix
    output_data = []

    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -0.5
    aligner.extend_gap_score = -0.1

    for row in data:
        tenx = row[0]
        static = row[1]
        mutate_seq = row[2]

        # Align mutate_seq with original_seq using PairwiseAligner
        alignments = aligner.align(original_seq, mutate_seq)

        # Choose the best alignment (first one in the list)
        best_alignment = alignments[0]

        # Create the state vector based on the best alignment
        state_vector = create_state_vector(best_alignment, original_seq_length)

        # Add row to output data
        output_data.append([tenx, static] + state_vector)

    # Create DataFrame for output
    columns = ['tenx', 'static'] + [f'V{i+1}' for i in range(original_seq_length)]
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