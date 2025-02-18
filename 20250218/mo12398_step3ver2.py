import pandas as pd
import numpy as np

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

# 修正処理をベクトル化
def correct_staticBC(df):
    # 条件に応じた修正を行う
    cond1 = df['staticBC'].str.len() <= 8
    cond2 = df['staticBC'].str.len() >= 12
    
    df.loc[cond1, 'new_staticBC'] = df['staticBC'] + df['mutableBC'].str[:6]
    df.loc[cond1, 'new_mutableBC'] = df['mutableBC'].str[6:]
    
    df.loc[cond2, 'new_staticBC'] = df['staticBC'].str[:10]
    df.loc[cond2, 'new_mutableBC'] = df['staticBC'].str[10:] + df['mutableBC']
    
    df.loc[~cond1 & ~cond2, 'new_staticBC'] = df['staticBC']
    df.loc[~cond1 & ~cond2, 'new_mutableBC'] = df['mutableBC']
    
    df['staticBC'] = df['new_staticBC']
    df['mutableBC'] = df['new_mutableBC']
    return df.drop(columns=['new_staticBC', 'new_mutableBC'])

# サンプルごとにデータを処理
for sample in samples:
    input_file = f"{base_path}/{sample}_processed_merge.txt"
    output_file = f"{base_path}/{sample}_processed_merge_staticBCcorrected.txt"
    
    # データの読み込み
    df = pd.read_csv(input_file, sep='\t')
    
    # ベクトル化された修正処理を適用
    df = correct_staticBC(df)
    
    # 結果を新しいファイルに出力（元の列名のまま）
    df.to_csv(output_file, sep='\t', index=False)
    
    print(f"Processed file saved to: {output_file}")