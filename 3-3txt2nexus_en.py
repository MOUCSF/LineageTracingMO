from Bio import SeqIO
import os

def convert_to_nexus(input_file):
    # 出力ファイル名を生成
    output_file = os.path.splitext(input_file)[0] + ".nexus"

    # DNAの配列を含む行列を読み込む
    data = {}
    with open(input_file, 'r') as f:
        header = f.readline().strip().split(',')
        for line in f:
            parts = line.strip().split(',')
            taxon = parts[0]
            sequence = ''.join(parts[2:-1])  # 最初の2つと最後の1つを除外
            data[taxon] = sequence

    # Nexus形式の出力
    with open(output_file, 'w') as f:
        f.write("#NEXUS\n")
        f.write("BEGIN DATA;\n")
        f.write("DIMENSIONS NTAX={} NCHAR={};\n".format(len(data), len(next(iter(data.values())))))
        f.write("FORMAT DATATYPE=DNA MISSING=? GAP=-;\n")
        f.write("MATRIX\n")
        for taxon, sequence in data.items():
            f.write("{} {}\n".format(taxon, sequence))
        f.write(";\n")
        f.write("END;\n")

    print("Conversion completed. Output file:", output_file)

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python convert_to_nexus.py input.txt")
        sys.exit(1)
    input_file = sys.argv[1]
    convert_to_nexus(input_file)