import pandas as pd
import numpy as np
from rapidfuzz import process, fuzz, utils #process.cdistの方が汎用性高いがHammingの方が楽
from rapidfuzz.distance import Hamming

# UMIcount
# Hamming Distanceを指定、距離を指定
hamming_threshold = 1  # 変更可能 (例: 1文字、2文字など)

# サンプル名のリストとデータパスの設定
samples = [
    "non1live_sub_S275_L006",
    #"non1sort_sub_S276_L006",
    "non2live_sub_S277_L006",
    #"non2sort_sub_S278_L006",
    "non3live_sub_S279_L006",
    #"non3sort_sub_S280_L006"
]
base_path="/Volumes/TRANSCEND/Gupta/MO12398"

# UMI同士のHamming距離を計算しグループ化する関数
def find_similar_groups(umi_list, threshold):
    groups = []
    used = set()
    
    for i, umi in enumerate(umi_list):
        if i in used:
            continue
        group = [umi]
        used.add(i)
        
        for j, candidate in enumerate(umi_list):
            if j not in used:
                # Hamming距離を計算して閾値内であればグループ化
                if Hamming.distance(umi, candidate) <= threshold:
                    group.append(candidate)
                    used.add(j)
        
        groups.append(group)
    
    return groups

# サンプルごとにデータを処理
for sample in samples:
    input_file = f"{base_path}/{sample}_processed_merge_staticBCcorrected.txt"
    output_file = f"{base_path}/{sample}_processed_merge_staticBCcorrected_UMIsort.txt"
    
    # データの読み込み
    df = pd.read_csv(input_file, sep='\t')
    
    results = []
    
    # 'cellBC'ごとにグループ化
    for cellBC, group in df.groupby('cellBC'):
        umi_list = group['UMI'].tolist()
        staticBC_list = group['staticBC'].tolist()
        mutableBC_list = group['mutableBC'].tolist()
        
        # 似たUMIでグループ化
        similar_groups = find_similar_groups(umi_list, hamming_threshold)
        
        # 各グループごとにまとめる
        for similar_group in similar_groups:
            umi_count = len(similar_group)
            # 出現回数が最も多いUMIを選択
            representative_umi = max(similar_group, key=similar_group.count) #similar_group[0]だと行の前の方
            indices = [umi_list.index(umi) for umi in similar_group]

            # representative_umiに対応するstaticBCとmutableBCを取得
            rep_index = umi_list.index(representative_umi)
            
            # 結果を保存
            results.append({
                'cellBC': cellBC,
                'UMI': representative_umi,
                'staticBC': staticBC_list[rep_index], #[indices[0]],
                'mutableBC': mutableBC_list[rep_index], #[indices[0]],
                'UMIcount': umi_count
            })
    
    # 結果を新しいデータフレームに変換
    result_df = pd.DataFrame(results)
    
    # 結果を新しいファイルに出力
    result_df.to_csv(output_file, sep='\t', index=False)
    
    print(f"Processed file saved to: {output_file}")