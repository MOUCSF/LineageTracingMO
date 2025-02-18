import os
import pandas as pd
import cassiopeia as cas

import matplotlib.pyplot as plt


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

# クローン番号を指定
clone_number = 2  # ここでクローン番号を設定, select appropriate clone_number


# Plotting and visualizing trees with Cassiopeia
allele_table_processed_path = os.path.join(base_path, f"{sample}_processed_merge_staticBCcorrected_UMIsort_normalized_clone_alleletable_organized.txt")
# Newick形式のツリーファイルを読み込む
newick_file_path = os.path.join(base_path, f"{sample}_Cassiopeia/vanilla_greedy_clone_{clone_number}.newick")


allele_table = pd.read_csv(allele_table_processed_path, sep='\t')
with open(newick_file_path, 'r') as f:
    tree = f.read()

# Initialize tree
tree = cas.data.CassiopeiaTree(tree=tree)

# Basic visualization
cas.pl.plot_matplotlib(tree, add_root=True)
plt.savefig(f"/Users/masahirookada/{sample}_clone_{clone_number}_vanilla_greedy_tree_circos.png")  # pngファイルに保存、パス変えてください。


cas.pl.plot_matplotlib(tree, orient='right', add_root=True)
plt.savefig(f"/Users/masahirookada/{sample}_clone_{clone_number}_vanilla_greedy_tree_right.png")

print(tree.nodes)

#clade_colors = {
#    '2|2|2|0|2|0|0|2|0|2|0|0|0|0|2|2|0|2|2|0|0|0|0|0|0|0|0|0|0': 'red',
#    '2|4|0|4|4|3|5|4|4|4|4|6|0|3|4|5|0|2|4|0|3|3|3|3|3|0|3|3|0': 'red'
#}

# 指定したclone_numberと同じLineageGroupを持つ行をフィルタリング
allele_table_clone_number = allele_table[allele_table['LineageGroup'] == clone_number]

cas.pl.plot_matplotlib(tree, orient='right', allele_table=allele_table_clone_number)#, clade_colors=clade_colors) # specific clone_number only
plt.savefig(f"/Users/masahirookada/{sample}_clone_{clone_number}_vanilla_greedy_tree_right_allele.png")