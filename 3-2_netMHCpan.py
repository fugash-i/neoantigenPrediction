import pandas as pd
from Bio.Seq import Seq
from pyfaidx import Fasta
import gffutils
import subprocess
import os
from tqdm import tqdm
import sys
import re
import yaml
import io

'''
0. 設定ファイルの読み込み
1. rMATSの実行
    1.1. コマンドの組み立て
    1.2. rMATSを実行
2. rMATSからnetMHCpanへのデータの受け渡し(変換)
    2.1. 必要ファイルの読み込み
        2.1.0. ファイルパスの指定
        2.1.1. rMATSの結果ファイル
        2.1.2. ゲノムFASTAファイル
        2.1.3. GTFアノテーションファイル
    2.2. rMATSデータのフィルタリング
        2.2.1. FDRとIncLevelDifferenceでフィルタリング, フィルタリング結果を書き出し
    2.3. rMATSデータの変換
        2.3.1. rMATSの各イベントの座標をmRNA配列に変換, 保存
        2.3.2. mRNA配列を翻訳, ペプチド配列としてFASTA形式で保存
3. netMHCpanの実行 # --- ココ ---
4. netMHCpanの結果をTPM(transcripts pre million)、RANK_EL(結合能)でフィルタリング
    4.1. TPMを計算
        4.1.1. 入出力ファイルの指定
        4.1.2. リードカウント
        4.1.3. サンプルごとにTPMを計算
    4.2. フィルタリング
        4.2.1. 各ファイルのパスと、フィルタリングの閾値を指定
        4.2.2. データの読み込み
        4.2.3. NetMHCpanの結果を結合能でフィルタリング
        4.2.4. NetMHCpanの結果をTPMでフィルタリング
        4.2.5. 最終結果をファイルに保存


'''

# 0. 設定ファイルの読み込み
try:
    with open('config.yaml', 'r') as f:
        config = yaml.safe_load(f)
    print("✅設定ファイルが正常に読み込まれました. ")
except FileNotFoundError:
    print("❌エラー: 設定ファイル 'config.yaml' が見つかりません。")
    sys.exit(1)

paths = config['paths']
params = config['params']
thresholds = config['thresholds']
b1_file = paths['b1_file'] 
b2_file = paths['b2_file'] 
gtf_file = paths['gtf_file'] 
base_dir = paths['base_dir']
output_dir = os.path.join(base_dir, "result")
tmp_dir = os.path.join(base_dir, "temp")
all_bam_files = paths['all_bam_files_dir']
sample_groups = config['sample_groups']
fasta_output_file_name = "abnormal_proteins_biopython.fasta"
fasta_output_file = os.path.join(output_dir, fasta_output_file_name)

# 3.netMHCpanの実行
# コマンドの組立
netmhcpan_executable = paths['netmhcpan_executable']

fasta_file = os.path.abspath(fasta_output_file) # 絶対パスで指定
netmhcpan_params = config['params']['netmhcpan']
hla_alleles = netmhcpan_params['hla_alleles']
peptide_lengths = netmhcpan_params['peptide_lengths']

netmhcpan_output_file_name = "netmhcpan_predictions.txt"
netmhcpan_output_file = os.path.join(output_dir, netmhcpan_output_file_name)

command = [
    netmhcpan_executable,
    '-f', fasta_file,
    '-a', hla_alleles,
    '-l', peptide_lengths,
    '-BA', # Binding Affinityの予測
    '-s'   # 結果をスコアでソート
]

# 実行前に、ファイルが存在するかチェック
if not os.path.exists(fasta_file):
    print(f"❌エラー: 入力ファイルが見つかりません: {fasta_file}")
    exit()

print("--- NetMHCpan実行 ---")
print("実行するコマンド:\n", " ".join(command))


# subprocessでコマンドを実行
try:
    print(f"\n予測を実行中...（ペプチドの数によっては時間がかかる場合があります）")
    
    # コマンドを実行し、標準出力と標準エラーをキャプチャする
    # check=True は、コマンドがエラーで終了した場合に例外を発生させるための設定
    # netMHCpanが自分の作業ディレクトリにあるdataフォルダ等を読み込めるように指定してやる(cwd = netmhcpan_dir)

    netmhcpan_dir = os.path.dirname(netmhcpan_executable)

    result = subprocess.run(
        command, 
        capture_output=True, 
        text=True, 
        check=True,
        cwd=netmhcpan_dir
    )
    # エラー出力を表示
    if result.stderr:
        print("--- NetMHCpan 標準エラー出力 (stderr) ---")
        print(result.stderr)
        print("------------------------------------------")
    if result.stderr:
        print("--- NetMHCpan 標準エラー出力 (stderr) [正常終了時] ---")
        print(result.stderr)
        print("-------------------------------------------------")

    # 実行結果をファイルに保存
    with open(netmhcpan_output_file, 'w') as f:
        f.write(result.stdout)
        
    print(f"\n✅ 実行が完了し、結果を {netmhcpan_output_file} に保存しました。")

except FileNotFoundError:
    print(f"❌エラー: '{netmhcpan_executable}' が見つかりません。パスが正しいか確認してください。")
    exit()
except subprocess.CalledProcessError as e:
    # NetMHCpanがエラーを返した場合の処理
    print(f"❌エラー: NetMHCpanの実行に失敗しました。")
    print("--- NetMHCpanからのエラーメッセージ ---")
    print(e.stderr) # NetMHCpanが出力したエラー内容を表示
    exit()
