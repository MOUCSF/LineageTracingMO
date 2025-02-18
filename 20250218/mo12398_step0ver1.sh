#!/bin/bash

#fastq preprocessing
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

# 検索する配列
sequenceR1="TAGTTACGCCAAGCTTGAATTC"
sequenceR2="GAATTCAAGCTTGGCGTAACTA"

# ミスマッチ許容数
mismatch=2

# 各サンプルに対して処理
for sample in "${samples[@]}"; do
    r1_file="${base_path}/${sample}_R1_001.fastq.gz"
    r2_file="${base_path}/${sample}_R2_001.fastq.gz"
    
    # 出力ファイル名
    output_r1="${base_path}/${sample}_${sequenceR1}_R1.fastq.gz"
    output_r2="${base_path}/${sample}_${sequenceR2}_R2.fastq.gz"
    
    # R1から特定配列を含むリードを抽出
    echo "Running: seqkit grep -p ${sequence} -m ${mismatch} -o ${output_r1} ${r1_file}"
    seqkit grep -p "${sequenceR1}" -m "${mismatch}" -o "${output_r1}" "${r1_file}"

    # R1のリードIDを一時ファイルに保存
    r1_id_file="${base_path}/${sample}_r1_ids.txt"
    seqkit seq -n -i "${output_r1}" > "${r1_id_file}"
    #seqkit seq -n -i "${output_r1}" | sed 's/ \[[12]:N:.*//g' > "${r1_id_file}"  # 最初の部分を抽出、スペースの後に1or2:N:で始まる文字を削除
    #(スペース): 最初の部分に空白を探します。
    #\[ : 開き角括弧をエスケープします。
    #[12] : 1または2のいずれかを探します。
    #:N: : 文字列:N:を探します。
    #.* : このパターンに続く任意の文字をすべて（0回以上）マッチさせます。

    # R2から特定配列の相補鎖を含むリードを抽出
    echo "Running: seqkit grep -p ${sequenceR2} -m ${mismatch} -o ${output_r2} ${r2_file}"
    seqkit grep -p "${sequenceR2}" -m "${mismatch}" -o "${output_r2}" "${r2_file}"

    # R2のリードIDを一時ファイルに保存
    r2_id_file="${base_path}/${sample}_r2_ids.txt"
    seqkit seq -n -i "${output_r2}" > "${r2_id_file}"

    # 出力ファイル名
    common_output_r1="${base_path}/${sample}_common_R1.fastq.gz"
    common_output_r2="${base_path}/${sample}_common_R2.fastq.gz"

    # R1とR2の共通IDを抽出
    common_id_file="${base_path}/${sample}_common_ids.txt"
    comm -12 <(sort "${r1_id_file}") <(sort "${r2_id_file}") > "${common_id_file}"

    # 共通するリードIDを使ってR1とR2のリードを抽出し、新しいファイルに保存
    echo "Running: seqkit grep -f ${common_id_file} -o ${common_output_r1} ${output_r1}"
    seqkit grep -f "${common_id_file}" -o "${common_output_r1}" "${output_r1}"

    echo "Running: seqkit grep -f ${common_id_file} -o ${common_output_r2} ${output_r2}"
    seqkit grep -f "${common_id_file}" -o "${common_output_r2}" "${output_r2}"

    # 一時ファイルの削除
    rm "${r1_id_file}" "${r2_id_file}" "${common_id_file}"

    echo "Common R1 and R2 reads for ${sample} have been saved to ${common_output_r1} and ${common_output_r2}."
done


#以下で実行
#chmod +x cassiopeia_daisypilot_step0.sh
#./cassiopeia_daisypilot_step0.sh