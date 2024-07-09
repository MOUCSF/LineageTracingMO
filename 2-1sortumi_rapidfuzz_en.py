import pandas as pd
from rapidfuzz import process, fuzz
import sys

# コマンドライン引数からファイル名を取得
input_file = sys.argv[1]
output_file = sys.argv[2]

# ファイルの読み込み
df = pd.read_csv(input_file, sep='\t')

# 重複を1文字ミスマッチを許容してまとめる関数
def find_similar_groups(umi_list):
    groups = []
    for umi in umi_list:
        found = False
        for group in groups:
            if fuzz.ratio(umi, group[0]) >= 80:  # 0.8 * 100
                group.append(umi)
                found = True
                break
        if not found:
            groups.append([umi])
    return groups

# 結果を保存するリスト
results = []

# 'tenx'ごとにグループ化
for tenx, group in df.groupby('tenx'):
    umi_list = group['umi'].tolist()
    static = group['static'].tolist()
    mutate = group['mutate'].tolist()
    
    # 似たUMIでグループ化
    similar_groups = find_similar_groups(umi_list)
    
    # 各グループごとにまとめる
    for similar_group in similar_groups:
        count = len(similar_group)
        representative_umi = similar_group[0]
        indices = [umi_list.index(umi) for umi in similar_group]
        
        # 結果を保存
        results.append({
            'tenx': tenx,
            'umi': representative_umi,
            'static': static[indices[0]],
            'mutate': mutate[indices[0]],
            'overlapumi': count
        })

# 結果をデータフレームに変換
result_df = pd.DataFrame(results)

# ファイルに出力
result_df.to_csv(output_file, sep='\t', index=False)