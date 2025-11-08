import subprocess
import os
import argparse
import logging

# ログ設定
logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s - %(levelname)s - %(message)s')

def process_sra(sra_id, output_dir, threads=4):
    """
    指定されたSRA IDのデータをダウンロードし、FASTQに変換する関数
    
    Args:
        sra_id (str): SRAランID (例: SRR1234567)
        output_dir (str): FASTQファイルの出力先ディレクトリ
        threads (int): fasterq-dumpが使用するスレッド数
    
    Returns:
        bool: 成功した場合はTrue、失敗した場合はFalse
    """
    logging.info(f"処理開始: {sra_id}")
    
    try:
        # 1. prefetchによるSRAファイルのダウンロード
        #    -O ./sra のように出力先を指定することも可能
        logging.info(f"[{sra_id}] SRAファイルのダウンロード中...")
        prefetch_cmd = ['prefetch', sra_id]
        subprocess.run(prefetch_cmd, check=True, capture_output=True, text=True)
        logging.info(f"[{sra_id}] ダウンロード完了。")

        # 2. fasterq-dumpによるFASTQファイルへの変換
        #    --split-files: ペアエンドデータの場合、_1.fastqと_2.fastqに分割
        #    -O: 出力ディレクトリ指定
        #    -p: 進捗表示
        #    -e: 使用するスレッド数
        logging.info(f"[{sra_id}] FASTQファイルへの変換中...")
        fasterq_dump_cmd = [
            'fasterq-dump',
            '--split-files',
            '--outdir', output_dir,
            '--progress',
            '--threads', str(threads),
            sra_id
        ]
        result = subprocess.run(fasterq_dump_cmd, check=True, capture_output=True, text=True)
        logging.info(f"[{sra_id}] 変換完了。出力先: {output_dir}")
        
        # (オプション) ダウンロードしたSRAファイルを削除してディスク容量を節約
        sra_file_path = os.path.expanduser(f"~/ncbi/public/sra/{sra_id}.sra")
        if os.path.exists(sra_file_path):
            os.remove(sra_file_path)
            logging.info(f"[{sra_id}] 中間SRAファイルを削除しました。")

        return True

    except subprocess.CalledProcessError as e:
        logging.error(f"[{sra_id}] でエラーが発生しました。")
        logging.error(f"コマンド: {' '.join(e.cmd)}")
        logging.error(f"エラー出力:\n{e.stderr}")
        return False
    except Exception as e:
        logging.error(f"[{sra_id}] で予期せぬエラーが発生しました: {e}")
        return False

def main():
    """
    メイン処理
    """
    parser = argparse.ArgumentParser(description="複数のSRA IDからFASTQファイルをダウンロード・生成するスクリプト")
    parser.add_argument('--sra_list', type=str, required=True, help="SRAランIDが1行ずつ記載されたテキストファイルのパス")
    parser.add_argument('--output_dir', type=str, required=True, help="FASTQファイルの出力先ディレクトリ")
    parser.add_argument('--threads', type=int, default=4, help="fasterq-dumpが使用するスレッド数")
    
    args = parser.parse_args()
    
    # 出力ディレクトリが存在しない場合は作成
    os.makedirs(args.output_dir, exist_ok=True)
    
    # 成功・失敗したIDを記録するリスト
    success_ids = []
    failed_ids = []
    
    # SRA IDリストファイルを読み込み
    try:
        with open(args.sra_list, 'r') as f:
            sra_ids = [line.strip() for line in f if line.strip()]
    except FileNotFoundError:
        logging.error(f"SRA IDリストファイルが見つかりません: {args.sra_list}")
        return

    logging.info(f"合計 {len(sra_ids)} 個のSRA IDの処理を開始します。")

    for sra_id in sra_ids:
        if process_sra(sra_id, args.output_dir, args.threads):
            success_ids.append(sra_id)
        else:
            failed_ids.append(sra_id)
            
    # 最終結果のサマリー
    logging.info("==================== 全ての処理が終了しました ====================")
    logging.info(f"成功: {len(success_ids)}件")
    logging.info(f"失敗: {len(failed_ids)}件")
    if failed_ids:
        logging.warning("失敗したSRA ID:")
        for failed_id in failed_ids:
            logging.warning(f"- {failed_id}")

if __name__ == '__main__':
    main()