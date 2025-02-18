import pandas as pd

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

# staticBCの割合（例えば50%の場合は0.5）
staticBC_threshold = 0.5

# 各サンプルに対して処理を実行
for sample in samples:
    input_file = f"{base_path}/{sample}_processed_merge_staticBCcorrected_UMIsort_normalized.txt"
    output_file = f"{base_path}/{sample}_processed_merge_staticBCcorrected_UMIsort_normalized_clone.txt"
    
    # データを読み込む
    df = pd.read_csv(input_file, sep="\t")
    
    # cellBCごとにUMInormalizedCountでソート
    sorted_df = df.sort_values(by=['cellBC', 'UMInormalizedCount'], ascending=[True, False])
    
    # lineageGroupを追加するための空のリスト
    lineage_groups = []
    lineage_group_map = {}  # lineageGroup番号を管理する辞書
    lineage_group_counter = 1  # lineageGroup番号のカウンター
    selected_staticBCs_map = {}  # cellBCごとのselected_staticBCsを保存

    # cellBCごとに処理を行う
    for cell_bc, group in sorted_df.groupby('cellBC'):
        total_count = group['UMInormalizedCount'].sum()
        cumulative_count = 0
        selected_staticBCs = []
        
        # UMInormalizedCountの累積を計算
        for _, row in group.iterrows():
            cumulative_count += row['UMInormalizedCount']
            selected_staticBCs.append(row['staticBC'])
            if cumulative_count / total_count >= staticBC_threshold:
                break
        
        # staticBCsのリストを保存
        selected_staticBCs_map[cell_bc] = (
            "_".join(sorted(selected_staticBCs)) if len(selected_staticBCs) > 1 else selected_staticBCs[0]
        )
        
        # selected_staticBCsをタプルとして他のcellBCとの一致を確認
        staticBC_tuple = tuple(sorted(selected_staticBCs))
        
        if staticBC_tuple not in lineage_group_map:
            # 新しいlineageGroupを作成
            lineage_group_map[staticBC_tuple] = lineage_group_counter
            lineage_group_counter += 1
        
        # 各staticBCに対してlineageGroupを付与
        lineage_group_number = lineage_group_map[staticBC_tuple]
        for staticBC in selected_staticBCs:
            lineage_groups.append((cell_bc, staticBC, int(lineage_group_number)))

    # lineageGroupをデータフレームに変換
    lineage_df = pd.DataFrame(lineage_groups, columns=['cellBC', 'staticBC', 'lineageGroup'])
    
    # 元のデータフレームにlineageGroupをマージ
    final_result = pd.merge(df, lineage_df, on=['cellBC', 'staticBC'], how='left')

    # lineageGroupが欠損している行に対して適切な番号を付与
    for cell_bc in final_result['cellBC'].unique():
        lineage_group = final_result.loc[final_result['cellBC'] == cell_bc, 'lineageGroup'].dropna().unique()
        if len(lineage_group) > 0:
            final_result.loc[final_result['cellBC'] == cell_bc, 'lineageGroup'] = lineage_group[0]

    # lineageGroupを多い順に並べ替え
    lineage_counts = final_result.groupby('lineageGroup')['cellBC'].nunique().reset_index()  # nunique()で重複を排除して数える
    lineage_counts.columns = ['lineageGroup', 'cellNumber']
    lineage_counts = lineage_counts.sort_values(by='cellNumber', ascending=False)
    
    # lineageGroupに番号を付け直す
    lineage_counts['newLineageGroup'] = range(1, len(lineage_counts) + 1)
    
    # 新しい番号とマージ
    final_result['lineageGroup'] = final_result['lineageGroup'].astype(int)
    final_result = final_result.merge(lineage_counts[['lineageGroup', 'newLineageGroup', 'cellNumber']], on='lineageGroup', how='left')
    
    # newLineageGroupを適切なlineageGroupに名前変更
    final_result['lineageGroup'] = final_result['newLineageGroup']
    final_result = final_result.drop(columns=['newLineageGroup'])

    # cellBCごとのselected_staticBCsを新しい列として追加
    final_result['selected_staticBCs'] = final_result['cellBC'].map(selected_staticBCs_map)

    # 結果を保存
    final_result.to_csv(output_file, sep="\t", index=False)

print("lineageGrouping completed.")
