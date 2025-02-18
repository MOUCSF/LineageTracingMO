import pandas as pd
import numpy as np
from rapidfuzz import process, fuzz, utils #process.cdistの方が汎用性高いがHammingの方が楽
from rapidfuzz.distance import Hamming

# UMInormalized-staticBCcount
# Hamming Distanceを指定、距離を指定
hamming_threshold = 1  # 変更可能 (例: 1文字、2文字など)

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

# staticBC同士のHamming距離を計算しグループ化する関数
def find_similar_groups(staticBC_list, threshold):
    groups = []
    used = set()
    
    for i, staticBC in enumerate(staticBC_list):
        if i in used:
            continue
        group = [staticBC]
        used.add(i)
        
        for j, candidate in enumerate(staticBC_list):
            if j not in used:
                # Hamming距離を計算して閾値内であればグループ化
                if Hamming.distance(staticBC, candidate) <= threshold:
                    group.append(candidate)
                    used.add(j)
        
        groups.append(group)
    
    return groups

# サンプルごとにデータを処理
for sample in samples:
    input_file = f"{base_path}/{sample}_processed_merge_staticBCcorrected_UMIsort.txt"
    output_file = f"{base_path}/{sample}_processed_merge_staticBCcorrected_UMIsort_normalized.txt"
    
    # データの読み込み
    df = pd.read_csv(input_file, sep='\t')

    # UMI数が2以上のものを抽出
    df_filtered = df[df['UMIcount'] >= 1]
    
    results = []
    
    # 'cellBC'ごとにグループ化
    for cellBC, group in df_filtered.groupby('cellBC'):
        staticBC_list = group['staticBC'].tolist()
        mutableBC_list = group['mutableBC'].tolist()
        umi_count_list = group['UMIcount'].tolist()
        
        # staticBCで似たものをグループ化
        similar_groups = find_similar_groups(staticBC_list, hamming_threshold)
        
        # 各グループごとにまとめる
        for similar_group in similar_groups:
            staticBC_count = len(similar_group)
            # 出現回数が最も多いstaticBCを選択
            representative_staticBC = max(similar_group, key=similar_group.count)
            indices = [staticBC_list.index(staticBC) for staticBC in similar_group]

            # representative_staticBCに対応するmutableBCとUMIcountを取得
            rep_index = staticBC_list.index(representative_staticBC)

            # UMIの数は行数
            umi_normalized_count = len(similar_group)            

            # UMIの出現回数を合計
            read_count = sum(umi_count_list[i] for i in indices)

            # 結果を保存
            results.append({
                'cellBC': cellBC,
                'staticBC': representative_staticBC,
                'mutableBC': mutableBC_list[rep_index],
                'UMInormalizedCount': umi_normalized_count,
                'ReadCount': read_count
            })
    
    # 結果を新しいデータフレームに変換
    result_df = pd.DataFrame(results)
    
    # 結果を新しいファイルに出力
    result_df.to_csv(output_file, sep='\t', index=False)
    
    print(f"Processed file saved to: {output_file}")