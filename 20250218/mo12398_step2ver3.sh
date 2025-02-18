#!/bin/bash

# サンプル名のリストを定義
samples=(
    "non1live_sub_S275_L006"
    #"non1sort_sub_S276_L006"
    "non2live_sub_S277_L006"
    #"non2sort_sub_S278_L006"
    "non3live_sub_S279_L006"
    #"non3sort_sub_S280_L006"
)

# データパスの基本パス
base_path="/Volumes/TRANSCEND/Gupta/MO12398"

# 各サンプルに対して処理を行う
for sample in "${samples[@]}"; do

    # ファイルパスの定義
    input_fasta_file="${base_path}/${sample}_merged_homology.fasta"
    input_txt_file="${base_path}/${sample}_merged_homology.txt"
    # ここでは '_' で区切られた最初のフィールドを取得　# cellBC_listファイル名は、接頭辞 + _main.whitelist.singlet.txt という形式
    sample_prefix="${sample%%_*}"
    cellBC_list="${base_path}/${sample_prefix}_main.whitelist.singlet.txt"

    seqkit fx2tab "$input_fasta_file" > "$input_txt_file" 

    # プロセスIDを格納する配列waitに使う
    pids=()

    ## 10xのbarcodeの抽出
    #awk '
    #NR==FNR { barcodes[$0]; next }
    #{
    #    found = 0
    #    for (barcode in barcodes) {
    #        if (index($0, barcode) > 0) {
    #            print barcode
    #            found = 1
    #            break
    #        }
    #    }
    #    if (!found) print "NA"
    #}' "$cellBC_list" "$input_txt_file" > "${base_path}/${sample}_processed_cellBC.txt" &
    #pids+=($!)

    # 10xのbarcodeと仮想UMIを同時抽出
    awk '
    NR==FNR { barcodes[$0]; next }
    {
        found_barcode = "NA"
        found_umi = "NA"
        for (barcode in barcodes) {
            pos = index($0, barcode)
            if (pos > 0) {
                found_barcode = barcode  # 一致したバーコードを記録
                found_umi = substr($0, pos + 16, 12)  # 16文字後から12文字をUMIとして抽出
                break  # 一致するバーコードが見つかればループを抜ける
            }
        }
    
        # バーコードとUMIを出力
        print (found_barcode == "NA" ? "NA" : found_barcode) "\t" (found_umi == "NA" ? "NA" : found_umi)
    }' "$cellBC_list" "$input_txt_file" > "${base_path}/${sample}_processed_cellBC_UMI.txt" &
    pids+=($!)

    # AGCTTGAATTCとTTTGCTに挟まれた配列抽出
    awk -F"\t" '
    {
        split($0, fields, "\t")
        if (length(fields) < 2) {
            print "NA"
            next
        }
        if ($2 ~ /AGCTTGAATTC.*TTTGCT/) {
            start_idx = index($2, "AGCTTGAATTC") + 11
            end_idx = index($2, "TTTGCT") - 1
            if (start_idx <= end_idx) {
                extracted_seq = substr($2, start_idx, end_idx - start_idx + 1)
                print extracted_seq
            } else {
                print "NA"
            }
        } else {
            print "NA"
        }
    }' "$input_txt_file" > "${base_path}/${sample}_processed_staticBC.txt" &
    pids+=($!)

    # TTTGCTより後ろのmutableBC配列抽出
    awk -F"\t" '
    {
        split($0, fields, "\t")
        if (length(fields) < 2) {
            print "NA"
            next
        }
        if ($2 ~ /TTTGCT(.*)/) {
            match($2, /TTTGCT.*/)
            print substr($2, RSTART, RLENGTH)
        } else {
            print "NA"
        }
    }' "$input_txt_file" > "${base_path}/${sample}_processed_mutableBC.txt" &
    pids+=($!)

    # すべてのバックグラウンド処理が終了するのを待つ
    for pid in "${pids[@]}"; do
        wait "$pid" || { echo "Process $pid failed"; exit 1; }
    done

    # ヘッダーを作成
    echo -e "cellBC\tUMI\tstaticBC\tmutableBC" > "${base_path}/${sample}_processed_merge.txt"

    # ファイルをマージして、`NA`を含む行を削除して追加
    paste "${base_path}/${sample}_processed_cellBC_UMI.txt" "${base_path}/${sample}_processed_staticBC.txt" "${base_path}/${sample}_processed_mutableBC.txt" | awk '!/NA/' >> "${base_path}/${sample}_processed_merge.txt"

done