import os
import shutil
from collections import Counter

def main(input_file, output_folder):
    # 出力先フォルダが存在しない場合は作成
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # データを読み込む
    with open(input_file, 'r') as f:
        lines = f.readlines()

    # ヘッダーとデータを分割
    header = lines[0].strip().split(',')
    data = [line.strip().split(',') for line in lines[1:]]

    # staticの列のインデックスを取得
    static_index = header.index('static')

    # staticごとのカウントを取得
    static_count = Counter([line[static_index] for line in data])

    # staticを多い順にソート
    sorted_statics = [static for static, _ in static_count.most_common()]

    # staticBC列を追加
    header.append('BCnumber')
    for i in range(len(data)):
        # staticごとに新しい名前を付けてBCnumber列を追加
        static = data[i][static_index]
        data[i].append(f"BC{sorted_statics.index(static) + 1}")

    # BC1のデータを抽出して別ファイルに保存
    bc1_data = [line for line in data if line[-1] == 'BC1']
    bc1_file = os.path.join(output_folder, "BC1.txt")
    with open(bc1_file, 'w') as f:
        # ヘッダー書き込み
        f.write(','.join(header) + '\n')
        # データ書き込み
        for line in bc1_data:
            f.write(','.join(line) + '\n')

    # BCall.txtを作成
    bcall_file = os.path.join(output_folder, "BCall.txt")
    with open(bcall_file, 'w') as f:
        # ヘッダー書き込み
        f.write(','.join(header) + '\n')
        # データ書き込み
        for line in data:
            f.write(','.join(line) + '\n')

    # staticごとにデータを分割して保存
    for i, static in enumerate(sorted_statics, start=1):
        # staticごとのデータを抽出
        static_data = [line for line in data if line[static_index] == static]
        # 出力ファイル名
        output_file = os.path.join(output_folder, f"BC{i}.txt")
        with open(output_file, 'w') as f:
            # ヘッダー書き込み
            f.write(','.join(header) + '\n')
            # データ書き込み
            for line in static_data:
                f.write(','.join(line) + '\n')

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python BCsort.py input_file output_folder")
        sys.exit(1)

    input_file = sys.argv[1]
    output_folder = sys.argv[2]
    main(input_file, output_folder)