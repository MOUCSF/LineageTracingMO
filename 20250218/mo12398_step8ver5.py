import subprocess
import pandas as pd
import tempfile
import os
import re
import concurrent.futures
from concurrent.futures import ThreadPoolExecutor

# サンプル名のリストとデータパスの設定
samples = [
    "non1live_sub_S275_L006",
    #"non1sort_sub_S276_L006",
    #"non2live_sub_S277_L006",
    #"non2sort_sub_S278_L006",
    #"non3live_sub_S279_L006",
    #"non3sort_sub_S280_L006"
]
base_path="/Volumes/TRANSCEND/Gupta/MO12398"

# select needle(global) or water(local), parameter gapopen 13,gapextend0.5
method = 'needle'
gapopen = '13'
gapextend = '0.5'

def run_water(seq1, seq2):
    with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.fasta') as f1, \
         tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.fasta') as f2:
        f1.write(f">seq1\n{seq1}")
        f2.write(f">seq2\n{seq2}")
        seq1_path = f1.name
        seq2_path = f2.name
    
    with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.txt') as f3:
        water_out = f3.name

    cmd = ["water", "-asequence", seq1_path, "-bsequence", seq2_path, "-gapopen", gapopen, "-gapextend", gapextend, "-outfile", water_out, "-aformat", "fasta"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        raise RuntimeError(f"water failed: {result.stderr}")

    with open(water_out, "r") as f:
        lines = f.readlines()

    os.remove(seq1_path)
    os.remove(seq2_path)
    os.remove(water_out)

    seq1_aligned = ""
    seq2_aligned = ""
    collecting_seq1 = False
    collecting_seq2 = False

    for line in lines:
        line = line.strip()
        if line.startswith(">seq1"):
            collecting_seq1 = True
            collecting_seq2 = False
        elif line.startswith(">seq2"):
            collecting_seq1 = False
            collecting_seq2 = True
        elif collecting_seq1:
            seq1_aligned += line
        elif collecting_seq2:
            seq2_aligned += line
    
    return seq1_aligned, seq2_aligned

def run_needle(seq1, seq2):
    with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.fasta') as f1, \
         tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.fasta') as f2:
        f1.write(f">seq1\n{seq1}")
        f2.write(f">seq2\n{seq2}")
        seq1_path = f1.name
        seq2_path = f2.name
    
    with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.txt') as f3:
        needle_out = f3.name

    cmd = ["needle", "-asequence", seq1_path, "-bsequence", seq2_path, "-gapopen", gapopen, "-gapextend", gapextend, "-outfile", needle_out, "-aformat", "fasta"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        raise RuntimeError(f"needle failed: {result.stderr}")

    with open(needle_out, "r") as f:
        lines = f.readlines()

    os.remove(seq1_path)
    os.remove(seq2_path)
    os.remove(needle_out)

    seq1_aligned = ""
    seq2_aligned = ""
    collecting_seq1 = False
    collecting_seq2 = False

    for line in lines:
        line = line.strip()
        if line.startswith(">seq1"):
            collecting_seq1 = True
            collecting_seq2 = False
        elif line.startswith(">seq2"):
            collecting_seq1 = False
            collecting_seq2 = True
        elif collecting_seq1:
            seq1_aligned += line
        elif collecting_seq2:
            seq2_aligned += line
    
    return seq1_aligned, seq2_aligned

def calculate_cigar(seq1, seq2):
    cigar = []
    count = 0
    current_op = ''
    
    for s1, s2 in zip(seq1, seq2):
        if s1 == s2:
            op = 'M'
        elif s1 == '-':
            op = 'D'
        elif s2 == '-':
            op = 'I'
        else:
            op = 'M'
        
        if op == current_op:
            count += 1
        else:
            if current_op:
                cigar.append(f"{count}{current_op}")
            current_op = op
            count = 1
    
    if current_op:
        cigar.append(f"{count}{current_op}")
    
    return ''.join(cigar)

def split_cigar_by_region(cigar, seq2_aligned, region_size=30, num_regions=4):
    regions_cigar = [''] * num_regions  # 各領域のCIGAR文字列を初期化
    region_limits = [region_size * (i + 1) for i in range(num_regions)]  # 各領域の参照位置の境界を設定

    ref_pos = 0  # 参照シーケンスの現在位置
    current_region = 0  # 現在の領域インデックス

    blocks = re.findall(r'(\d+)([MID])', cigar)

    for length_str, op in blocks:
        length = int(length_str)  # 長さは常に整数型として扱う

        while length > 0 and current_region < num_regions:
            region_end = region_limits[current_region]
            remaining_in_region = region_end - ref_pos

            if op in ('M', 'D'):
                assign_length = min(length, remaining_in_region)
                regions_cigar[current_region] += f"{assign_length}{op}"
                ref_pos += assign_length
                length -= assign_length

                if ref_pos == region_end:
                    current_region += 1

            elif op == 'I':
                regions_cigar[current_region] += f"{length}{op}"
                length = 0

    return regions_cigar


def main(input_file, output_file, method=method):
    input_df = pd.read_csv(input_file, sep='\t')
    reference_sequence = 'TTTGCTACCTATTACTAGGACAAGTGGAGGCGTTGAGAATGTCTCGGAATGTGCCGGAAATTTGCATCGTATTACTAGGACAATTCGAGTCCTTGAGGAAGTCTCTCGATGTGCCGGAAA'

    def process_sequence(mutate_seq):
        if method == 'needle':
            seq1_aligned, seq2_aligned = run_needle(mutate_seq, reference_sequence)
        elif method == 'water':
            seq1_aligned, seq2_aligned = run_water(mutate_seq, reference_sequence)
        else:
            raise ValueError("Invalid method. Choose 'needle' or 'water'.")

        cigar = calculate_cigar(seq1_aligned, seq2_aligned)
        regions_cigar = split_cigar_by_region(cigar, seq2_aligned)
        return cigar, regions_cigar

    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = list(executor.map(process_sequence, input_df['mutableBC']))

    cigars = [result[0] for result in results]
    region1_cigars = [result[1][0] for result in results]
    region2_cigars = [result[1][1] for result in results]
    region3_cigars = [result[1][2] for result in results]
    region4_cigars = [result[1][3] for result in results]

    input_df['whole_cigar'] = cigars
    input_df['region1_cigar'] = region1_cigars
    input_df['region2_cigar'] = region2_cigars
    input_df['region3_cigar'] = region3_cigars
    input_df['region4_cigar'] = region4_cigars

    output_df = input_df[['cellBC', 'staticBC', 'mutableBC', 'whole_cigar', 'region1_cigar', 'region2_cigar', 'region3_cigar', 'region4_cigar', 'lineageGroup', 'selected_staticBCs', 'cellNumber', 'UMInormalizedCount', 'ReadCount']]
    output_df.to_csv(output_file, sep='\t', index=False)


if __name__ == "__main__":
    # 各サンプルに対して処理を実行
    for sample in samples:
        input_file = f"{base_path}/{sample}_processed_merge_staticBCcorrected_UMIsort_normalized_clone.txt"
        output_file = f"{base_path}/{sample}_processed_merge_staticBCcorrected_UMIsort_normalized_clone_alleletable.txt"
        main(input_file, output_file, method)

        # 実行完了のメッセージを表示
        alignment = "needle" if method == "needle" else "water"
        print(f"{sample}, {alignment} alignment to CIGAR.")

