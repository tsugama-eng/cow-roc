import pandas as pd
from Bio import SeqIO
from collections import Counter
import sys
import argparse
import os
import csv

def main():
    parser = argparse.ArgumentParser(description="Strict Common-Gene Ribo-seq Context Aggregator")
    parser.add_argument("psite_tsv", help="Path to the filtered P-site density TSV file")
    parser.add_argument("fasta_file", help="Path to the filtered CDS FASTA file")
    parser.add_argument("-o", "--output", default="codon_context_score.tsv", help="Output filename")
    parser.add_argument("--range", type=int, default=15, help="Analysis range")
    parser.add_argument("--no_bg", action="store_true", help="Disable background codon count normalization")

    args = parser.parse_args()

    # 1. リードデータ（事前に絞り込み済み）のIDをスキャン
    print(f"Step 1: Scanning gene IDs in filtered lead data...")
    lead_gene_ids = set()
    with open(args.psite_tsv, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            lead_gene_ids.add(row['transcript_id'])

    # 2. FASTAからリードデータに存在する遺伝子のみを抽出・背景カウント
    # フィルタリング済みのファイルが入力されるため、ここでの一致数が「解析対象の全遺伝子」となる
    print(f"Step 2: Counting background for common genes found in both files...")
    fasta_dict = {}
    genomic_codon_counts = Counter()
    
    for rec in SeqIO.parse(args.fasta_file, "fasta"):
        tid = rec.id
        if tid in lead_gene_ids:
            seq = str(rec.seq).upper()
            fasta_dict[tid] = seq
            # 共通遺伝子の全域コドンを背景としてカウント
            for i in range(0, len(seq) - 2, 3):
                codon = seq[i:i+3]
                if len(codon) == 3:
                    genomic_codon_counts[codon] += 1

    # 3. 集計処理（2段階正規化の「遺伝子内正規化」ステップ）
    print(f"Step 3: Processing P-site data and performing transcript-level normalization...")
    site_range = range(-args.range, args.range + 1)
    position_scores = {offset: Counter() for offset in site_range}
    
    current_tid = None
    current_group = []
    mismatch_total = 0

    def process_group(tid, rows):
        if not rows or tid not in fasta_dict:
            return 0
        
        seq = fasta_dict[tid]
        counts = [r['count'] for r in rows]
        avg_reads = sum(counts) / len(counts)
        if avg_reads == 0: return 0

        mismatches = 0
        for r in rows:
            pos = r['pos']
            norm_score = r['count'] / avg_reads
            p_idx = pos - 1
            
            if p_idx < 0 or p_idx + 3 > len(seq) or seq[p_idx : p_idx + 3] != r['expected']:
                mismatches += 1
                continue
            
            for offset in site_range:
                start_idx = p_idx + (offset * 3)
                end_idx = start_idx + 3
                if start_idx >= 0 and end_idx <= len(seq):
                    codon = seq[start_idx:end_idx]
                    if len(codon) == 3:
                        position_scores[offset][codon] += norm_score
        return mismatches

    with open(args.psite_tsv, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            tid = row['transcript_id']
            data = {'pos': int(row['cds_position']), 'count': float(row['count']), 'expected': row['codon_seq'].upper()}

            if tid != current_tid:
                if current_tid is not None:
                    mismatch_total += process_group(current_tid, current_group)
                current_tid = tid
                current_group = [data]
            else:
                current_group.append(data)
        mismatch_total += process_group(current_tid, current_group)

    # 4. 背景補正および結果テーブルの構築（2026年最適化版）
    if args.no_bg:
        print(f"Step 4: Accumulating scores (Background correction OFF)...")
    else:
        print(f"Step 4: Applying background correction (O/E ratio based on common genes)...")

    bases = ['T', 'C', 'A', 'G']
    all_codons = [a+b+c for a in bases for b in bases for c in bases]
    
    # 1列ずつ DataFrame に追加するのではなく、まずリストに Series を溜める
    all_columns = []

    for offset in sorted(site_range):
        label = f"pos_{offset}"
        column_data = {}
        
        for codon in all_codons:
            score = position_scores[offset].get(codon, 0)
            
            if args.no_bg:
                # 背景補正なし
                column_data[codon] = score
            else:
                # 背景補正（共通遺伝子内での出現頻度で割る）
                bg_count = genomic_codon_counts.get(codon, 0)
                if bg_count > 0:
                    column_data[codon] = score / bg_count
                else:
                    column_data[codon] = 0.0
        
        # 1列分の Series を作成してリストに追加
        all_columns.append(pd.Series(column_data, name=label))

    # 全ての列を一気に結合（axis=1）して DataFrame を作成
    # これにより Fragmented DataFrame の警告を回避し、処理が高速化される
    result_table = pd.concat(all_columns, axis=1)

    # 5. 保存
    result_table.index.name = "codon"
    result_table.to_csv(args.output, sep='\t')

    # 最終レポートの表示
    common_count = len(fasta_dict)
    missing_ids = lead_gene_ids - set(fasta_dict.keys())
    
    print(f"\n" + "="*45)
    print(f"Final Report (Year: 2026 Analysis)")
    print(f"IDs in Lead TSV:      {len(lead_gene_ids)}")
    print(f"Common Genes Used:    {common_count} (Primary variants)")
    if missing_ids:
        print(f"IDs not in FASTA:     {len(missing_ids)} (Skipped)")
    print(f"Mismatches (P-site):  {mismatch_total}")
    print(f"Background Correction: {'OFF' if args.no_bg else 'ON (Common-gene base)'}")
    print(f"Output Saved to:      {args.output}")
    print("="*45)

if __name__ == "__main__":
    main()
