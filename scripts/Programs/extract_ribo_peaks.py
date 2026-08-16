import pandas as pd
import numpy as np
from scipy.stats import poisson
import csv
import sys
import argparse

def main():
    parser = argparse.ArgumentParser(description="Extract significant Ribo-seq peaks using Poisson test")
    parser.add_argument("input_tsv", help="Path to original ribodensity TSV")
    parser.add_argument("-o", "--output", default="significant_peaks.tsv", help="Output TSV filename")
    parser.add_argument("--p_threshold", type=float, default=0.001, help="Adjusted P-value threshold (default: 0.001)")
    
    args = parser.parse_args()

    # Step 1: 遺伝子ごとの平均リード密度を計算 (1回目のスキャン)
    print("Step 1: Calculating mean density per transcript...")
    gene_means = {}
    gene_site_counts = {}
    
    with open(args.input_tsv, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            tid = row['transcript_id']
            cnt = float(row['count'])
            gene_means[tid] = gene_means.get(tid, 0) + cnt
            gene_site_counts[tid] = gene_site_counts.get(tid, 0) + 1
    
    for tid in gene_means:
        gene_means[tid] /= gene_site_counts[tid]

    # Step 2: 生P値の計算 (2回目のスキャン)
    print("Step 2: Calculating Poisson P-values...")
    p_values = []
    data_buffer = [] # メモリ節約のため必要な情報のみ一時保持
    
    with open(args.input_tsv, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            tid = row['transcript_id']
            mu = gene_means[tid]
            cnt = float(row['count'])
            
            # ポアソン検定: 平均muのときにcnt以上の値が出る確率
            # sf (survival function) は 1 - cdf
            p_val = poisson.sf(cnt - 1, mu) if cnt > 0 else 1.0
            p_values.append(p_val)
            data_buffer.append(row)

    # Step 3: 全ゲノム規模での多重比較補正 (BH法)
    print("Step 3: Correcting P-values (Benjamini-Hochberg)...")
    from statsmodels.stats.multitest import multipletests
    rejected, adj_p, _, _ = multipletests(p_values, alpha=args.p_threshold, method='fdr_bh')

    # Step 4: 有意なピークのみを出力
    print(f"Step 4: Writing significant peaks to {args.output}...")
    with open(args.output, 'w', newline='') as f:
        fieldnames = ['transcript_id', 'cds_position', 'codon_seq', 'count']
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        
        peak_count = 0
        for i, is_peak in enumerate(rejected):
            if is_peak:
                # 元のデータ形式を維持して出力
                writer.writerow({
                    'transcript_id': data_buffer[i]['transcript_id'],
                    'cds_position': data_buffer[i]['cds_position'],
                    'codon_seq': data_buffer[i]['codon_seq'],
                    'count': data_buffer[i]['count']
                })
                peak_count += 1

    print(f"\nDONE! Identified {peak_count} significant peaks.")
    print(f"Threshold: Adjusted P < {args.p_threshold}")

if __name__ == "__main__":
    main()
