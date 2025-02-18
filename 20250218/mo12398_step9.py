import pandas as pd

# alleletable to cassiopeia input 
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

# 各サンプルに対して処理を行う
for sample in samples:
    input_file = f"{base_path}/{sample}_processed_merge_staticBCcorrected_UMIsort_normalized_clone_alleletable.txt"
    output_file = f"{base_path}/{sample}_processed_merge_staticBCcorrected_UMIsort_normalized_clone_alleletable_organized.txt"
    
    # データを読み込む
    df = pd.read_csv(input_file, sep="\t")
    
    # 列名を変更
    df = df.rename(columns={
        'staticBC': 'intBC',
        'region1_cigar': 'r1',
        'region2_cigar': 'r2',
        'region3_cigar': 'r3',
        'region4_cigar': 'r4',
        'mutableBC': 'allele',
        'lineageGroup': 'LineageGroup',
        'ReadCount': 'readCount',
        'UMInormalizedCount': 'UMI'
    })
    
    # sampleID 列を追加
    df['sampleID'] = sample
    
    # 列を最終的な順序に並び替え
    df = df[['cellBC', 'intBC', 'r1', 'r2', 'r3', 'r4', 'allele', 'LineageGroup', 'sampleID', 'readCount', 'UMI']]
    
    # 新しいCSVファイルに保存
    df.to_csv(output_file, index=False, sep="\t")
    
    print(f"Processed {sample} and saved to {output_file}")