from IPython.display import Image

import numpy as np
import pandas as pd

import cassiopeia as cas

import os

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

# read in an allele table
sample = samples[0]  # サンプルを指定0で一つ目、1で二つ目

allele_table_processed_path = os.path.join(base_path, f"{sample}_processed_merge_staticBCcorrected_UMIsort_normalized_clone_alleletable_organized.txt")

allele_table = pd.read_csv(allele_table_processed_path, sep='\t',
                           usecols = ['cellBC', 'intBC', 'r1', 'r2', 'r3', 'r4', 'allele', 'LineageGroup', 'sampleID', 'readCount', 'UMI'])

# データ型を変換（float（数字）とstr（文字）が混在しているカラムを str に変換）
allele_table['r1'] = allele_table['r1'].astype(str)
allele_table['r2'] = allele_table['r2'].astype(str)
allele_table['r3'] = allele_table['r3'].astype(str)
allele_table['r4'] = allele_table['r4'].astype(str)
allele_table['allele'] = allele_table['allele'].astype(str)

# 保存フォルダを作成
output_dir = os.path.join(base_path, f"{sample}_Cassiopeia")
os.makedirs(output_dir, exist_ok=True)

#Estimating indel priors
indel_priors = cas.pp.compute_empirical_indel_priors(allele_table, grouping_variables=['intBC', 'LineageGroup'])

# priortable を保存
priortable_output_path = os.path.join(output_dir, "indel_priors.txt")
indel_priors.to_csv(priortable_output_path, sep='\t', index=True)#Trueにしないと列名indel出ない
print(f"Indel priors saved to {priortable_output_path}")

# すべてのクローン（LineageGroup）の番号を取得
unique_clones = allele_table['LineageGroup'].unique()

# 各クローンごとにツリーを作成
for CLONE in unique_clones:
    clone_allele_table = allele_table[allele_table['LineageGroup'] == CLONE]

    n_cells = clone_allele_table['cellBC'].nunique()
    n_intbc = clone_allele_table['intBC'].nunique()
    print(f"Clonal population #{CLONE} has {n_cells} cells and {n_intbc} intBCs ({n_intbc * 4}) characters.")

    # Character matrixを作成
    character_matrix, priors, state_2_indel = cas.pp.convert_alleletable_to_character_matrix(
        clone_allele_table,
        allele_rep_thresh=0.9,  # 90%以上で同じものは無意味として削除
        mutation_priors=indel_priors
    )

    # character matrix が空でないかを確認
    if character_matrix.empty:
        print(f"Skipping clone {CLONE} due to empty character matrix.")
        continue
    if n_cells < 5:  # セル数が少なすぎる場合スキップ
        print(f"Skipping clone {CLONE} due to insufficient cellBCs or intBCs.")
        continue

    # Character matrix を保存
    character_matrix_output_path = os.path.join(output_dir, f"clone_{CLONE}_character_matrix.txt")
    character_matrix.to_csv(character_matrix_output_path, sep='\t', index=True)#Trueにしないと列名cellBC出ない
    print(f"Character matrix for clone {CLONE} saved to {character_matrix_output_path}")

    # ツリーを構築
    cas_tree = cas.data.CassiopeiaTree(character_matrix=character_matrix, priors=priors)

    # Vanilla Greedy Solverでツリーを構築
    vanilla_greedy = cas.solver.VanillaGreedySolver()
    vanilla_greedy.solve(cas_tree, collapse_mutationless_edges=True)

    # Newick フォーマットでツリーを保存
    tree_output_path = os.path.join(output_dir, f"vanilla_greedy_clone_{CLONE}.newick")
    newick_tree = cas_tree.get_newick()

    with open(tree_output_path, 'w') as f:
        f.write(newick_tree)

    print(f"Tree for clone {CLONE} saved to {tree_output_path}")