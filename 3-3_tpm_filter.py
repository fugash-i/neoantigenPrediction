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
3. netMHCpanの実行
4. netMHCpanの結果をTPM(transcripts pre million)、RANK_EL(結合能)でフィルタリング # --- ココカラ ---
    4.1. TPMを計算
        4.1.1. 入出力ファイルの指定
        4.1.2. リードカウント
        4.1.3. サンプルごとにTPMを計算
    4.2. フィルタリング
        4.2.1. 各ファイルのパスと、フィルタリングの閾値を指定
        4.2.2. データの読み込み
        4.2.3. NetMHCpanの結果を結合能でフィルタリング
        4.2.4. NetMHCpanの結果をTPMでフィルタリング
        4.2.5. 最終結果をファイルに保存 # --- ココマデ ---


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
paths = config['paths']
params = config['params']
thresholds = config['thresholds']
base_dir = paths['base_dir']
output_dir = os.path.join(base_dir, "result")
tmp_dir = os.path.join(base_dir, "temp")
rmats_output = os.path.join(output_dir, "rmats")
gtf_file = paths['gtf_file'] 
rmats_file = os.path.join(rmats_output, "SE.MATS.JCEC.txt")


# 4. netMHCpanの結果をTPM(transcripts pre million)、RANK_EL(結合能)でフィルタリング

# 4.1. TPMを計算
# 4.1.1. 入出力ファイルの指定
# 入力する全てのBAMファイルのリスト
all_bam_files = paths['all_bam_files_dir']

# 平均を計算したいサンプルグループを定義(キー: グループ名) 
# 値: グループに属するサンプルのベース名 (パスや拡張子を除く) のリスト
sample_groups = config['sample_groups']

# 中間ファイル: featureCountsの出力
counts_output_file_name = 'gene_read_counts_all.txt'
counts_output_file = os.path.join(tmp_dir, counts_output_file_name)

# 最終的な出力ファイル: 遺伝子IDとグループごとの平均TPM値
final_tpm_output_file_name = 'gene_tpm_average.tsv'
final_tpm_output_file =os.path.join(tmp_dir, final_tpm_output_file_name)

# 4.1.2. リードカウント
# featureCountsを実行してリードカウントを取得
featurecounts_executable = 'featureCounts'
command = [
    featurecounts_executable,
    '-p',
    '-a', gtf_file,
    '-o', counts_output_file,
    '-T', '4',
    '-g', 'gene_id',
] + all_bam_files

print("\n--- featureCountsによるリードカウントを実行 ---")
try:
    subprocess.run(command, check=True, capture_output=True, text=True)
    print(f"✅ リードカウントが完了し、{counts_output_file} に保存されました。")
except Exception as e:
    print("❌エラー: featureCountsの実行に失敗しました。")
    print(e.stderr if hasattr(e, 'stderr') else e)
    exit()

# 4.1.3. サンプルごとにTPMを計算
print("\n--- サンプルごとのTPMを計算 ---")
df_counts = pd.read_csv(counts_output_file, sep='\t', comment='#')

df_tpm_all_samples = pd.DataFrame()
df_tpm_all_samples['GeneID'] = df_counts['Geneid']

# BAMファイルごとに（列ごとに）TPMを計算するループ
for bam_path in all_bam_files:
    # featureCountsの出力では、列名はフルパスになる
    counts = df_counts[bam_path]
    lengths = df_counts['Length']

    # TPMを計算
    rpk = counts / (lengths / 1000)
    per_million_scale_factor = rpk.sum() / 1_000_000
    if per_million_scale_factor == 0: # ゼロ除算を避ける
        tpm = 0
    else:
        tpm = rpk / per_million_scale_factor
    
    # サンプル名を綺麗にして（拡張子などを除いて）列名とする
    sample_name = os.path.basename(bam_path).replace('.bam', '')
    df_tpm_all_samples[f'{sample_name}_TPM'] = tpm

print("✅ 全サンプルのTPM計算が完了しました。")

# 4.1.4. グループごとにTPMの平均値を計算
print("\n--- グループごとの平均TPMを計算 ---")
df_tpm_final = pd.DataFrame()
df_tpm_final['GeneID'] = df_tpm_all_samples['GeneID']
print("利用可能なTPMカラム:", df_tpm_all_samples.columns)

for group_name, sample_list in sample_groups.items():
    # 平均計算の対象となる列名のリストを作成 (例: ['tumor_sample1_TPM', ...])
    columns_to_average = [f'{sample_name}_TPM' for sample_name in sample_list]
    
    # 平均値を計算し、新しい列として追加
    df_tpm_final[f'{group_name}_TPM_mean'] = df_tpm_all_samples[columns_to_average].mean(axis=1)

print("✅ 平均TPMの計算が完了しました。")

# tpmの計算結果をファイルに保存
df_tpm_final.to_csv(final_tpm_output_file, sep='\t', index=False)

print(f"\n最終的な平均TPMデータを {final_tpm_output_file} に保存しました。")
print("\n---出力データの一部---")
print(df_tpm_final.head())


# 4.2. フィルタリング
# 4.2.1. 各ファイルのパスと、フィルタリングの閾値を指定
netmhcpan_file_name = 'netmhcpan_predictions.txt'
netmhcpan_file = os.path.join(output_dir, netmhcpan_file_name)

avg_tpm_file_name = 'gene_tpm_average.tsv'
avg_tpm_file = os.path.join(tmp_dir, avg_tpm_file_name)

# フィルタリングの閾値
thresholds = config['thresholds']
rank_el_threshold = thresholds['rank_el']
tpm_threshold = thresholds['tpm']

# 4.2.2. データの読み込み
# NetMHCpanの出力はヘッダーが複雑なため、データ部分だけを読み込む処理
print("フィルタリングに必要な各ファイルを読み込みます")
print(f"\nNetMHCpanファイル {netmhcpan_file} を読み込んでいます...")
try:
    with open(netmhcpan_file, 'r') as f:
        lines = f.readlines()

    # フィルタリング: コメント行、空行、区切り線 (---)、繰り返しヘッダー (Protein) を除外
    cleaned_lines = []
    header_found = False # 最初のヘッダーを見つけたかどうかのフラグ
    
    for line in lines:
        stripped_line = line.strip()
        
        # 不要な行（コメント、区切り線、集計行、空行）をスキップ
        if not stripped_line or stripped_line.startswith(('#', '---', 'Protein')):
            continue
        
        # ヘッダー行 (Pos...) の処理
        if stripped_line.startswith('Pos'):
            if not header_found:
                # 1. 最初のヘッダー行だけをリストに追加
                cleaned_lines.append(line)
                header_found = True
            
            # 2. 2回目以降のヘッダー行はスキップする
            continue 
        
        # ★★★ ここが新しい修正点 ★★★
        # " <= SB" や " <= WB" のように分離している列を、
        # " <=SB" のように1つのフィールドに結合する
        # これで、18列ではなく17列として正しく認識される
        if " <= SB" in line:
            line = line.replace(" <= SB", " <=SB")
        elif " <= WB" in line:
            line = line.replace(" <= WB", " <=WB")
        
        # データ行を追加
        cleaned_lines.append(line)

    # フィルタリングした行を一つの文字列に結合し、pandasで読み込む
    df_mhc = pd.read_csv(
        io.StringIO("".join(cleaned_lines)), 
        # delim_whitespace=True は非推奨であり、
        # このケース(<= SB)では列の分割を解決できなかった
        sep=r'\s+', # 警告(FutureWarning)の推奨に従い、正規表現の空白区切りを使用
        engine='python' # 正規表現セパレータにはpythonエンジンが推奨される
    )
    print("✅ NetMHCpanファイルの読み込みと整形が完了しました。")

except Exception as e:
    print(f"❌ NetMHCpanファイルの読み込み中にエラーが発生しました: {e}")
    sys.exit(1)


# 平均TPMファイルとrMATSファイルを読み込む
print("TPMファイル/rMATSファイルを読み込んでいます...")
df_tpm_avg = pd.read_csv(avg_tpm_file, sep='\t')
df_rmats = pd.read_csv(rmats_file, sep='\t')
print("✅TPMファイル/rMATSファイルの読み込みが完了しました")


# 4.2.3. NetMHCpanの結果を結合能でフィルタリング
df_mhc['%Rank'] = pd.to_numeric(df_mhc['%Rank'], errors='coerce')
df_mhc.dropna(subset=['%Rank'], inplace=True)
df_mhc_filtered = df_mhc[df_mhc['%Rank'] < rank_el_threshold].copy()
print(f"NetMHCpanの結果を%Rank < {rank_el_threshold} でフィルタリングしました. フィルタリング後のイベント: {len(df_mhc_filtered)} 件")

# 4.2.4. NetMHCpanの結果をTPMでフィルタリング
# データをマージするための準備
# NetMHCpanのID列 (例: event_1234|...) からイベントID (1234) を抽出
df_mhc_filtered['event_ID'] = df_mhc_filtered['Identity'].apply(lambda x: int(re.search(r'event_(\d+)', x).group(1)))

# rMATSのデータフレームから、イベントID (行番号) とGeneIDのマッピングを作成
df_id_map = df_rmats[['ID', 'GeneID']].reset_index().rename(columns={'index': 'event_ID'})


# 全データをマージし、TPMで最終フィルタリング
# (1) NetMHCpanの結果に、イベントIDをキーとしてGeneIDを結合
df_merged = pd.merge(df_mhc_filtered, df_id_map[['event_ID', 'GeneID']], on='event_ID', how='inner')

# (2) GeneIDをキーとして平均TPMデータを結合
df_final = pd.merge(df_merged, df_tpm_avg, on='GeneID', how='inner')

# (3) treatmentサンプルの平均TPM値でフィルタリング
final_candidates = df_final[df_final['treatment_TPM_mean'] >= tpm_threshold]
print(f"TPM >= {tpm_threshold} でフィルタリング: {len(final_candidates)} 件")


# 4.2.5. 最終結果をファイルに保存
output_file_name = 'final_high_confidence_neoantigens.csv'
output_file = os.path.join(output_dir, output_file_name)
# 表示したい列を整理して保存
columns_to_save = [
    'Peptide', 'MHC', '%Rank', 'Aff(nM)', 
    'GeneID', 'treatment_TPM_mean', 'Identity'
]
final_candidates[columns_to_save].sort_values(by='%Rank').to_csv(output_file, index=False)

print(f"\n最終的なネオ抗原候補を {output_file} に保存しました。")
print("\n--- 最終候補の一部 ---")
print(final_candidates[columns_to_save].sort_values(by='%Rank').head())