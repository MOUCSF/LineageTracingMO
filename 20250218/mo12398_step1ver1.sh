#!/bin/bash

# merge fastq R1 and R2 by fastp 
# サンプル名のリストを定義、bashでは,とかspace不要
samples=(
    "non1live_sub_S275_L006"
    "non1sort_sub_S276_L006"
    "non2live_sub_S277_L006"
    "non2sort_sub_S278_L006"
    "non3live_sub_S279_L006"
    "non3sort_sub_S280_L006"
)

# データパスの基本パス
base_path="/Volumes/TRANSCEND/Gupta/MO12398"


# 検索する配列を直接変数に指定
sequence="TAGTTACGCCAAGCTTGAATTC"



# 各サンプルに対してfastpを実行
for sample in "${samples[@]}"; do
  echo "Processing $sample"

  # 入力ファイルの指定
  R1="${base_path}/${sample}_common_R1.fastq.gz"
  R2="${base_path}/${sample}_common_R2.fastq.gz"

  # 出力ファイルの指定
  unmerged_R1="${base_path}/${sample}_unmerged_R1.fastq.gz"
  unmerged_R2="${base_path}/${sample}_unmerged_R2.fastq.gz"
  merged="${base_path}/${sample}_merged.fastq.gz"

  # fastpの実行
  fastp -i "$R1" -I "$R2" -o "$unmerged_R1" -O "$unmerged_R2" --merge --merged_out "$merged" -q 20 -u 30 -n 5 -e 20

  # seqkit grepを実行して、特定の配列を含むリードを抽出
  grep_output="${base_path}/${sample}_merged_homology.fastq.gz"
  echo "Running: seqkit grep -P -p ${sequence} -m 2 -j 12 ${merged} -o ${grep_output}"
  seqkit grep -P -p "${sequence}" -m 2 -j 12 "${merged}" -o "${grep_output}"

  # 抽出されたFASTQファイルをFASTA形式に変換
  fasta_output="${base_path}/${sample}_merged_homology.fasta"
  echo "Running: seqkit fq2fa ${grep_output} -o ${fasta_output}"
  seqkit fq2fa "${grep_output}" -o "${fasta_output}"

  echo "Finished processing ${sample}."

done
