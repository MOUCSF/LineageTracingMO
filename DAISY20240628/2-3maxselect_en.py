import sys
import pandas as pd

def main(input_file, output_file):
    # 入力ファイルの読み込み
    df = pd.read_csv(input_file, sep='\t')

    # 各 tenx のグループ内で overlapmutate が最大の行を選択
    result = df.loc[df.groupby('tenx')['overlapmutate'].idxmax()]

    # 結果を出力ファイルに書き込み
    result.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    # コマンドライン引数のチェック
    if len(sys.argv) != 3:
        print("Usage: python maxselect.py input.txt output.txt")
        sys.exit(1)

    # コマンドライン引数の取得
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # メイン関数の実行
    main(input_file, output_file)
