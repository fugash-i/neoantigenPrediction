import os
import glob
import subprocess
import sys
import yaml


try:
    with open('/Users/kent/vscode/config.yaml', 'r') as f:
        config = yaml.safe_load(f)
except FileNotFoundError:
    print("エラー: 設定ファイル 'config.yaml' が見つかりません。")
    sys.exit(1)

paths = config['paths']
params = config['params']

INPUT_DIR = paths['fastq_dir']
OUTPUT_DIR = paths['mapping_output_dir']
INDEX_PATH = paths['hisat2_index_path']

THREADS = params['hisat2']['threads']

def run_command(command):
    """コマンドを実行し、エラーがあれば停止する関数"""
    print(f"実行中: {' '.join(command)}")
    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        print(f"エラーが発生しました: {e}")
        exit(1)

def main():
    """メインの処理"""
    print("解析を開始します...")

    # 出力用のサブディレクトリを作成
    sam_dir = os.path.join(OUTPUT_DIR, "sam")
    sorted_bam_dir = os.path.join(OUTPUT_DIR, "bam_sorted")
    os.makedirs(sam_dir, exist_ok=True)
    os.makedirs(sorted_bam_dir, exist_ok=True)

    # 入力ディレクトリ内のペアエンドリードの片方 (*_1.fastq.gz) を見つけてループ処理
    # glob.globは条件に合うファイルパスのリストを返す
    gz_files = glob.glob(os.path.join(INPUT_DIR, "*_1.fastq.gz"))
    fq_files = glob.glob(os.path.join(INPUT_DIR, "*_1.fastq"))
    fastq1_files = gz_files + fq_files # 2つのリストを結合

    if not fastq1_files:
        print(f"エラー: {INPUT_DIR} 内に *_1.fastq または *_1.fastq.gz ファイルが見つかりません。")
        return

    for fastq1_path in fastq1_files:
        # --- ファイル名とサンプル名の準備 ---
        base_name = os.path.basename(fastq1_path)
        if base_name.endswith("_1.fastq.gz"):
            sample_name = base_name.replace("_1.fastq.gz", "")
            fastq2_path = fastq1_path.replace("_1.fastq.gz", "_2.fastq.gz")
        elif base_name.endswith("_1.fastq"):
            sample_name = base_name.replace("_1.fastq", "")
            fastq2_path = fastq1_path.replace("_1.fastq", "_2.fastq")
        else:
            print(f"警告: 予期しないファイル名のためスキップします: {base_name}")
            continue # ループの次の回へ
        
        # 出力ファイルパスの準備
        sam_file = os.path.join(sam_dir, f"{sample_name}.sam")
        bam_file = os.path.join(sam_dir, f"{sample_name}.bam")
        sorted_bam_file = os.path.join(sorted_bam_dir, f"{sample_name}_sorted.bam")

        # --- 処理状況を表示 ---
        print("\n" + "="*50)
        print(f">> 処理開始: サンプル名 {sample_name}")
        print("="*50)

        # 1. HISAT2によるマッピング
        print("Step 1: HISAT2でマッピングを実行中...")
        cmd_hisat2 = [
            "hisat2",
            "-p", str(THREADS),
            "--dta",
            "-x", INDEX_PATH,
            "-1", fastq1_path,
            "-2", fastq2_path,
            "-S", sam_file
        ]
        run_command(cmd_hisat2)

        # 2. SAMをBAMに変換
        print("\nStep 2: SAMをBAMに変換中...")
        cmd_view = [
            "samtools", "view",
            "-bS",
            "-o", bam_file,
            sam_file
        ]
        run_command(cmd_view)

        # 3. BAMファイルを座標でソート
        print("\nStep 3: BAMファイルをソート中...")
        cmd_sort = [
            "samtools", "sort",
            "-o", sorted_bam_file,
            bam_file
        ]
        run_command(cmd_sort)

        # 4. ソート済みBAMファイルのインデックスを作成
        print("\nStep 4: インデックスを作成中...")
        cmd_index = [
            "samtools", "index",
            sorted_bam_file
        ]
        run_command(cmd_index)

        # メモリ節約のため、中間ファイルのSAMとソート前BAMを削除（必要ならこの行をコメントアウト）
        os.remove(sam_file)
        os.remove(bam_file)
        print(f"\n>> 中間ファイルを削除しました: {os.path.basename(sam_file)}, {os.path.basename(bam_file)}")


        print(f"\n>> 処理完了: サンプル名 {sample_name}")

    print("\n" + "="*50)
    print("🎉 全てのサンプルの処理が完了しました！")
    print("="*50)

# スクリプトが直接実行された場合にmain()関数を呼び出す
if __name__ == "__main__":
    main()