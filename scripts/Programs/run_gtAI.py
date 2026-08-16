#!/usr/bin/env python3
"""
Run gtAI analysis with command-line arguments.

Usage:
    python run_gtAI.py input.fasta tRNA_counts.txt
"""

import sys
import pandas as pd
import warnings
warnings.simplefilter(action="ignore", category=FutureWarning)
import logging
logging.getLogger("gaft").setLevel(logging.WARNING)

from gtAI import Run_gtAI, gtAI


def read_tRNA_file(path):
    GtRNA = {}
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            anticodon, count = parts[0], parts[1]
            try:
                GtRNA[anticodon] = int(count)
            except ValueError:
                # ヘッダーや文字列が来た場合はスキップ
                continue
    return GtRNA

def main():
    if len(sys.argv) < 3:
        print("Usage: python run_gtai.py <input_fasta> <tRNA_file>")
        sys.exit(1)

    main_fasta = sys.argv[1]   # 入力FASTAファイル
    tRNA_file = sys.argv[2]    # tRNAコピー数ファイル

    # パラメータ設定（必要に応じて調整可能）
    ref_fasta = ""
    genetic_code_number = 1
    bacteria = False
    size_pop = 60
    generation_number = 100

    # tRNAコピー数の読み込み
    # 例: ファイルは "anticodon<TAB>copy_number" の形式を想定
    #GtRNA = gtAI.read_tRNA_file(tRNA_file)
    GtRNA = read_tRNA_file(tRNA_file)

    # gtAI解析の実行
    df_tai, final_dict_wi, rel_values = Run_gtAI.gtai_analysis(
        main_fasta=main_fasta,
        GtRNA=GtRNA,
        ref_fasta=ref_fasta,
        genetic_code_number=genetic_code_number,
        size_pop=size_pop,
        generation_number=generation_number,
        bacteria=bacteria
    )

    # 結果を保存
    #df_tai.to_csv("gtAI_results.csv", header=True)
    #print("gtAI analysis finished. Results saved to gtAI_results.csv")
    print(df_tai.to_string(index=False))
    
if __name__ == "__main__":
    main()