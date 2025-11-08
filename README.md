# neoantigenPrediction
Predicting splice-neoantigens (induced by chemical treatment or mutation of splicing factor) from RNA-seq data using rMATS and NetMHCPan. 

## 概要
sraファイル(1でダウンロード)をBAMに変換, ソートし(2)たのちにrMATSを用いてASEを検出, FDR, incLevelDifference, TPMでフィルタリングした上でタンパク質配列に翻訳(3), netMHCpanで特定のMHC class Iに提示されると予測されるアミノ酸配列予測を出力する(3).  
(3)のプロセスはやや時間がかかる上複雑になるため3-1~3-3のコードを順に用いても同様のことができるようになっている. 

## config.yamlについて
入力/出力ファイルのディレクトリ, HISAT2/rMATS/netMHCpanのオプションなどを設定する. 詳しくは当該ファイルのコメントを参照. 


## 各生成ファイルの保存先("*_dir"はconfig.yamlを編集して指定): 

| ファイル | 保存先 |
| -- | -- |
| RNA-seqのマップ済データ(bamファイル) | mapping_output_dir |
| rMATSの一時ファイル | base_dir/temp |
| rMATSの結果ファイル | base_dir/result/rmats |
| GTFを変換したデータベース(annotation.db) | base_dir/temp/annotation.db |
| rMATSで有意に増加したCEI eventのmRNA配列 | base_dir/temp/rmats_processed_sequences.pkl |
| rmats_processed_sequences.pklの翻訳結果(アミノ酸配列) | base_dir/result/abnormal_proteins_biopython |
| netMHCpanの実行結果 | base_dir/result/netmhcpan_predictions.txt |
| TPM値 | base_dir/temp/gene_tpm_average.tsv |
| 最終的な予測結果 | base_dir/result/final_high_confidence_neoantigens |
