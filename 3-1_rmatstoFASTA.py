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
1. rMATSの実行 # --- ココカラ ---
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
        2.3.2. mRNA配列を翻訳, ペプチド配列としてFASTA形式で保存 # --- ココマデ ---
3. netMHCpanの実行
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

# 1. rMATSの実行

# 1.1. コマンドの組み立て

# 入力
b1_file = paths['b1_file'] 
b2_file = paths['b2_file'] 
gtf_file = paths['gtf_file'] 

# 出力
base_dir = paths['base_dir']
output_dir = os.path.join(base_dir, "result")
tmp_dir = os.path.join(base_dir, "temp")

# rMATSのパラメータ
rmats_params = config['params']['rmats']
read_type   = rmats_params['read_type'] 
read_length = rmats_params['read_length']
n_threads   = rmats_params['n_threads']
lib_type    = rmats_params['lib_type']
rmats_output = os.path.join(output_dir, "rmats")

# rMATSコマンドをリスト形式で組み立てる
command = [
    'rmats.py',
    '--b1', b1_file,
    '--b2', b2_file,
    '--gtf', gtf_file,
    '-t', read_type,
    '--readLength', str(read_length),  # 数値は文字列に変換する必要がある
    '--nthread', str(n_threads),        # 同上
    '--libType', lib_type,
    '--od', rmats_output,
    '--tmp', tmp_dir,
    '--variable-read-length'
]

# 1.2. rMATSを実行(時間がかかる)
print("--- 以下のコマンドでrMATSを実行(時間がかかります) ---")
print(" ".join(command))
print("-----------------------------------------------")


try:
    result = subprocess.run(
        command,
        capture_output=True,  # 標準出力と標準エラーをキャプチャする
        text=True,            # 出力をテキストとしてデコードする
        check=True            # コマンドがエラーで終了した場合、例外を発生させる
    )

    print("--- rMATSが正常に完了しました ---")
    print("--- 標準出力 (stdout) ---")
    print(result.stdout)
    if result.stderr:
        print("--- 標準エラー出力 (stderr) ---")
        print(result.stderr)

except subprocess.CalledProcessError as e:
    # check=Trueによって、rMATSがエラー終了した際にこのブロックが実行される
    print("--- rMATSの実行中にエラーが発生しました ---")
    print(f"終了コード: {e.returncode}")
    print("--- 標準出力 (stdout) ---")
    print(e.stdout)
    print("--- 標準エラー出力 (stderr) ---")
    print(e.stderr)
    sys.exit(1)

except FileNotFoundError:
    # 'rmats.py'というコマンド自体が見つからなかった場合
    print("--- エラー ---")
    print("コマンド 'rmats.py' が見つかりませんでした。")
    print("rMATSは適切にインストールされていますか？")
    print("このスクリプトを実行する前に、ターミナルでrMATS用のConda環境を有効にしましたか？")
    sys.exit(1)

# 2. rMATSからnetMHCpanへのデータの受け渡し(変換)
print("rMATSデータを加工してnetMHCpanに受け渡すプロセスを開始します")

# 2.1. 必要ファイルの読み込み
print("--- ファイルの読み込み ---")

# 2.1.0. ファイルパスの指定
rmats_file = os.path.join(rmats_output, "SE.MATS.JCEC.txt")
genome_fasta = paths['genome_fasta']
annotation_gtf = paths['gtf_file']
db_path = os.path.join(tmp_dir, 'annotation.db') # GTFを変換したデータベースファイル

# 2.1.1. rMATSの結果ファイル
try:
    df_rmats = pd.read_csv(rmats_file, sep='\t')
    print(f"\n✅ rMATSファイルを読み込みました。 {len(df_rmats)} イベント")
    print("rMATSデータの一部:")
    print(df_rmats.head())
except FileNotFoundError:
    print(f"❌エラー: rMATSファイルが見つかりません: {rmats_file}")
    exit()

# 2.1.2. ゲノムFASTAファイル
try:
    genome = Fasta(genome_fasta)
    print(f"\n✅ ゲノムFASTAファイルを読み込みました。 {len(genome.keys())} 個の染色体")
except FileNotFoundError:
    print(f"❌エラー: ゲノムFASTAファイルが見つかりません: {genome_fasta}")
    exit()

# 2.1.3. GTFアノテーションファイル
# 巨大なGTFファイルをデータベースにして読み込み時間を短縮
try:
    if not os.path.exists(db_path):
        print(f"\nGTFデータベースが存在しないため、新規作成します...（時間がかかることがあります）")
        db = gffutils.create_db(
            annotation_gtf, 
            db_path, 
            )
    else:
        print(f"\n既存のGTFデータベースを読み込みます...")
        db = gffutils.FeatureDB(db_path)
    
    gene_count = db.execute("SELECT COUNT(*) FROM features WHERE featuretype = 'gene'").fetchone()[0]
    print(f"✅ GTFデータベースを読み込みました。 {gene_count} 個の遺伝子")

except FileNotFoundError:
     print(f"❌エラー: GTFファイルが見つかりません: {annotation_gtf}")
     exit()

print("\n--- 全てのファイルの読み込みが完了しました. ---")

# 2.2. データのフィルタリング

# 2.2.1. FDRとIncLevelDifferenceでフィルタリング, フィルタリング結果を書き出し

print("--- FDRとPSIでフィルタリングを行います. ---")
print(f"\nフィルタリング前のイベント数: {len(df_rmats)}")

# 閾値の設定
fdr_threshold = thresholds['fdr']
inc_level_diff_threshold = thresholds['inc_level_diff']

# フィルタリングの実行
significant_events = df_rmats[
    (df_rmats['FDR'] < fdr_threshold) &
    (df_rmats['IncLevelDifference'] > inc_level_diff_threshold)
].copy() 

significant_events_csv_file = os.path.join(output_dir, "significant_splicing_events.csv")
significant_events.to_csv(significant_events_csv_file, index=False)
print(f"✅フィルタリングされた {len(significant_events)} 件のイベントを {significant_events_csv_file} に保存しました。")
print("\n--- フィルタリングが完了しました. ---")

# 2.3. データの変換

# 2.3.1. rMATSの各イベントの座標をmRNA配列に変換, 保存

tqdm.pandas(desc="Processing Events") # プロセスバーのための関数

# rMATSの各イベントについて, inclusion isoformの全長CDS配列およびcassette exonの位置を返す関数を定義

def get_abnormal_sequence(event_row, db, genome):
    try:
        chrom = event_row['chr']
        strand = event_row['strand']
        gene_id = event_row['GeneID']

        # rMATSのイベント座標
        cassette_exon_range = (event_row['exonStart_0base'] + 1, event_row['exonEnd']) # 1-based

        # GTFから最適なtranscriptを見つける
        target_transcript = None
        for t in db.children(gene_id, featuretype='transcript', order_by='start'):
            # transcriptがcassette exon領域と重なるCDSを持つ(inclusion isoform)かチェック
            has_cassette_cds = False
            for cds in db.children(t.id, featuretype='CDS'):
                overlap_start = max(cds.start, cassette_exon_range[0])
                overlap_end = min(cds.end, cassette_exon_range[1])
                if overlap_start < overlap_end:
                    has_cassette_cds = True
                    break
            
            if has_cassette_cds:
                target_transcript = t
                break # 最初に見つかったInclusion transcriptを使用

        if not target_transcript:
            return None, "No suitable inclusion transcript found in GTF"

        # 全長CDS配列の構築と、座標計算: gffutilsがframe情報を考慮して正しい順序でCDSを返す
        transcript_cds_list = list(db.children(target_transcript.id, featuretype='CDS', order_by='start'))

        if not transcript_cds_list:
            return None, "No CDS found in the chosen transcript"

        full_mrna_seq_parts = []
        upstream_cds_len = 0
        cassette_cds_len = 0
        
        passed_cassette_start = False   # フラグ: cassette exonを通過したか

        for cds in transcript_cds_list:
            # 完全なmRNA配列の断片を取得
            current_cds_seq = genome[cds.chrom][cds.start - 1:cds.end].seq
            full_mrna_seq_parts.append(current_cds_seq)

            # カセットエクソンとの位置関係を計算 (ストランドを考慮)
            
            # (A) このCDSがカセットと「重なる」部分
            cas_range = (max(cds.start, cassette_exon_range[0]), min(cds.end, cassette_exon_range[1]))
            if cas_range and cas_range[0] <= cas_range[1]:
                current_cassette_len = (cas_range[1] - cas_range[0] + 1)
                cassette_cds_len += current_cassette_len
                
                # ストランドに応じて、カセット開始フラグを立てるタイミングが異なる: プラス鎖はカセットに入った時点でフラグON
                if strand == '+':
                    passed_cassette_start = True 
                
            # (B) このCDSがカセットより「上流(transcript-wise)」にある部分 = カセット開始フラグがまだOFFであるか？
            if not passed_cassette_start:
                 # このCDSのうち、カセットと重ならなかった部分の長さを足す
                 # (cas_range[0] <= cas_range[1] の場合、cas_lenは上で計算済み)
                 total_len = (cds.end - cds.start + 1)
                 if cas_range and cas_range[0] <= cas_range[1]:
                     upstream_cds_len += (total_len - current_cassette_len)
                 else:
                     upstream_cds_len += total_len

            # (C) マイナス鎖の場合、カセットの「終わり」を通過した時点でフラグON
            if strand == '-' and cas_range and cas_range[0] <= cas_range[1]:
                passed_cassette_start = True 

        # --- 3. 最終的なmRNA配列の構築 (ストランド処理) ---
        full_mrna_seq = "".join(full_mrna_seq_parts)
        
        if strand == '-':
            full_mrna_seq = str(Seq(full_mrna_seq).reverse_complement())

        # --- 4. 戻り値 ---
        return {
            "full_mrna": full_mrna_seq,         # 完全なCDS配列
            "upstream_cds_len": upstream_cds_len, # カセットまでの正確な塩基長
            "cassette_cds_len": cassette_cds_len  # カセット部分の正確な塩基長
        }, "Success: Full CDS"

    except Exception as e:
        return None, str(e)
    
# 定義した関数を用いてCDSのmRNA配列を抽出
print(f"全 {len(significant_events)} 件のイベントを処理します...（時間がかかる場合があります）")

results = significant_events.progress_apply(
    lambda row: get_abnormal_sequence(row, db, genome),
    axis=1
)

print("\n--- 全てのイベントの処理が完了しました ---")

# 結果の保存
rmats_output_file_name = "rmats_processed_sequences.pkl"
rmats_output_file = os.path.join(tmp_dir, rmats_output_file_name)
results.to_pickle(rmats_output_file)

print(f"結果を {rmats_output_file} に保存しました。結果を表示します。")

print(results)

# 2.3.2. mRNA配列を翻訳, ペプチド配列としてFASTA形式で保存

print("保存した結果を翻訳して新規ペプチド配列を書き出します。")

# 保存した結果を読み込む
pkl_file_name = "rmats_processed_sequences.pkl"
pkl_file = os.path.join(tmp_dir, pkl_file_name)
results_series = pd.read_pickle(pkl_file)

# Biopythonを使って翻訳し、FASTAを作成
fasta_output_file_name = "abnormal_proteins_biopython.fasta"
fasta_output_file = os.path.join(output_dir, fasta_output_file_name)
fasta_entries = []

CONTEXT_AMINO_ACIDS = 10 # 前後10アミノ酸をコンテキストとして含める

for index, (result_data, message) in results_series.items():
    if message == "Success: Full CDS" and result_data is not None:
        
        # 必要なデータを取得
        full_mrna = result_data['full_mrna']
        up_len = result_data['upstream_cds_len'] # カセット開始までのCDS塩基長
        cas_len = result_data['cassette_cds_len'] # カセット部分のCDS塩基長

        if not full_mrna or len(full_mrna) < 3:
            continue 

        # mRNAを翻訳
        # 翻訳前に、長さを3の倍数にトリミング
        full_mrna_trimmed = full_mrna[:len(full_mrna) - (len(full_mrna) % 3)]
        
        if len(full_mrna_trimmed) < 3:
            continue 

        # 翻訳の実行
        full_protein = Seq(full_mrna_trimmed).translate(to_stop=True)
        
        if not full_protein:
            continue 

        # アミノ酸座標の計算
        # カセットエクソンが始まるタンパク質上のインデックス (0-based)
        # up_len (塩基長) を 3 で割ることで、アミノ酸座標が正確に求まる
        protein_cassette_start_idx = up_len // 3
        
        # カセットエクソンが終わるタンパク質上のインデックス (Pythonスライス用)
        protein_cassette_end_idx_exclusive = (up_len + cas_len + 2) // 3
        
        # 関心領域（カセット ± 10aa）の切り出し
        slice_start = max(0, protein_cassette_start_idx - CONTEXT_AMINO_ACIDS)
        slice_end = min(len(full_protein), protein_cassette_end_idx_exclusive + CONTEXT_AMINO_ACIDS)

        target_peptide_seq = full_protein[slice_start : slice_end]

        # FASTAに追加
        sequence = str(target_peptide_seq)
        if len(sequence) >= 8: # NetMHCpanが受け付けられる最低長
            # ヘッダーにイベントIDと、タンパク質上の座標情報を含める
            header = f">event_{index}|inclusion_isoform|prot_coords_{slice_start}-{slice_end}"
            fasta_entries.append(f"{header}\n{sequence}")

# FASTAファイルに書き出す
with open(fasta_output_file, 'w') as f:
    f.write("\n".join(fasta_entries))

print(f"Biopythonを使い、{len(fasta_entries)}個の配列を {fasta_output_file} に保存しました。")
