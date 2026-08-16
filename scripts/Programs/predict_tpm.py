#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pickle
import sys
import os

# --- Amino Acid to Codon Table ---
AA_TABLE = {
    'A': ['GCT', 'GCC', 'GCA', 'GCG'], 'C': ['TGT', 'TGC'], 'D': ['GAT', 'GAC'],
    'E': ['GAA', 'GAG'], 'F': ['TTT', 'TTC'], 'G': ['GGT', 'GGC', 'GGA', 'GGG'],
    'H': ['CAT', 'CAC'], 'I': ['ATT', 'ATC', 'ATA'], 'K': ['AAA', 'AAG'],
    'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'M': ['ATG'],
    'N': ['AAT', 'AAC'], 'P': ['CCT', 'CCC', 'CCA', 'CCG'], 'Q': ['CAA', 'CAG'],
    'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
    'T': ['ACT', 'ACC', 'ACA', 'ACG'], 'V': ['GTT', 'GTC', 'GTA', 'GTG'],
    'W': ['TGG'], 'Y': ['TAT', 'TAC'], '*': ['TAA', 'TAG', 'TGA']
}

CODONS = sorted([a+b+c for a in "ATGC" for b in "ATGC" for c in "ATGC"])
CODON_INDEX = {c: i for i, c in enumerate(CODONS)}


# ---------- 基本ユーティリティ ----------

def compute_codon_freq(seq: str) -> dict:
    seq = seq.upper().replace("U", "T")
    usable = (len(seq) // 3) * 3
    seq = seq[:usable]
    counts = {codon: 0 for codon in CODONS}
    for i in range(0, usable, 3):
        codon = seq[i:i+3]
        if codon in counts:
            counts[codon] += 1
    total = sum(counts.values())
    if total == 0:
        total = 1e-6
    return {k: v / total for k, v in counts.items()}


# ---------- モデル1個用：TPM最大化（位置ごとに独立） ----------

# ---------- 差分CLRを使った高速最適化 ----------

def clr_from_counts_fast(counts):
    """counts: np.array shape (64,)"""
    total = counts.sum()
    if total == 0:
        total = 1e-6
    p = counts / total
    p = np.maximum(p, 1e-6)
    logp = np.log(p)
    gm = logp.mean()
    return logp - gm


def predict_from_clr_fast(clr_vec, model, log_transform):
    y = model.intercept_ + np.dot(model.coef_, clr_vec)
    if log_transform == "log1p":
        return np.expm1(y)
    elif log_transform == "log10":
        return 10 ** y
    return y

def compute_counts(seq):
    seq = seq.upper().replace("U", "T")
    usable = (len(seq) // 3) * 3
    seq = seq[:usable]
    counts = np.zeros(len(CODONS), dtype=float)
    for i in range(0, usable, 3):
        codon = seq[i:i+3]
        if codon in CODON_INDEX:
            counts[CODON_INDEX[codon]] += 1
    return counts

def optimize_single_model_fast(gene_id, initial_cds, model, log_transform, max_iters=100):
    """
    単一モデル用の高速最適化。
    1箇所変更するたびにカウントとCLR状態を更新し、配列全体で変更がなくなるまで最大max_iters回スキャンする。
    """
    cds = initial_cds.upper().replace("U", "T")
    usable = (len(cds) // 3) * 3
    cds = cds[:usable]
    aa_seq = str(Seq(cds).translate(to_stop=False))

    counts = compute_counts(cds)
    
    # 初期スコアの計算
    clr_vec = clr_from_counts_fast(counts)
    current_tpm = predict_from_clr_fast(clr_vec, model, log_transform)

    for it in range(max_iters):
        improved_any_in_this_pass = False
        
        for pos, aa in enumerate(aa_seq):
            if aa not in AA_TABLE or len(AA_TABLE[aa]) <= 1:
                continue
            
            start = pos * 3
            end = start + 3
            old_codon = cds[start:end]
            old_idx = CODON_INDEX[old_codon]

            best_codon_at_pos = old_codon
            best_tpm_at_pos = current_tpm

            # 当該位置の全コドンを試行
            for new_codon in AA_TABLE[aa]:
                if new_codon == old_codon:
                    continue

                new_idx = CODON_INDEX[new_codon]

                # 差分更新を試す
                counts[old_idx] -= 1
                counts[new_idx] += 1

                new_clr = clr_from_counts_fast(counts)
                new_tpm = predict_from_clr_fast(new_clr, model, log_transform)

                if new_tpm > best_tpm_at_pos:
                    best_tpm_at_pos = new_tpm
                    best_codon_at_pos = new_codon
                
                # 一旦戻す（最善であれば後で反映）
                counts[old_idx] += 1
                counts[new_idx] -= 1

            # 改善が見つかったら即座に配列とカウントを更新
            if best_codon_at_pos != old_codon:
                new_idx = CODON_INDEX[best_codon_at_pos]
                counts[old_idx] -= 1
                counts[new_idx] += 1
                cds = cds[:start] + best_codon_at_pos + cds[end:]
                current_tpm = best_tpm_at_pos
                improved_any_in_this_pass = True

        # 配列全体を1周スキャンして、どこも変更されなければ収束とみなして終了
        if not improved_any_in_this_pass:
            break

    return cds, current_tpm

def optimize_multi_model_fast(
    gene_id, initial_cds, models, log_transform, max_iters=100  # デフォルト回数を増やす
):
    cds = initial_cds.upper().replace("U", "T")
    usable = (len(cds) // 3) * 3
    cds = cds[:usable]
    aa_seq = str(Seq(cds).translate(to_stop=False))

    counts = compute_counts(cds)
    
    # 状態の初期化
    clr_vec = clr_from_counts_fast(counts)
    current_preds = [predict_from_clr_fast(clr_vec, m, log_transform) for m in models]
    current_min = min(current_preds)

    for it in range(max_iters):
        improved_any_in_this_pass = False
        
        for pos, aa in enumerate(aa_seq):
            if aa not in AA_TABLE or len(AA_TABLE[aa]) <= 1:
                continue
            
            start, end = pos * 3, (pos + 1) * 3
            old_codon = cds[start:end]
            old_idx = CODON_INDEX[old_codon]

            best_codon_at_pos = old_codon
            best_min_at_pos = current_min

            for new_codon in AA_TABLE[aa]:
                if new_codon == old_codon:
                    continue

                new_idx = CODON_INDEX[new_codon]
                # 差分更新を試す
                counts[old_idx] -= 1
                counts[new_idx] += 1

                new_clr = clr_from_counts_fast(counts)
                new_preds = [predict_from_clr_fast(new_clr, m, log_transform) for m in models]
                new_min = min(new_preds)

                if new_min > best_min_at_pos:
                    best_min_at_pos = new_min
                    best_codon_at_pos = new_codon
                
                # 一旦戻す
                counts[old_idx] += 1
                counts[new_idx] -= 1

            # もし、この位置で改善が見つかったら即座に配列とカウントを更新する
            if best_codon_at_pos != old_codon:
                new_idx = CODON_INDEX[best_codon_at_pos]
                counts[old_idx] -= 1
                counts[new_idx] += 1
                cds = cds[:start] + best_codon_at_pos + cds[end:]
                current_min = best_min_at_pos
                improved_any_in_this_pass = True

        # 配列全体を1周スキャンして、1箇所も変更がなければ「真の収束」とみなす
        if not improved_any_in_this_pass:
            break

    # 最終予測
    final_clr = clr_from_counts_fast(counts)
    final_preds = [predict_from_clr_fast(final_clr, m, log_transform) for m in models]
    return cds, min(final_preds), final_preds

# ---------- main ----------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta")
    parser.add_argument("model_pickles", nargs='+')
    parser.add_argument("--extra_features", default=None)
    parser.add_argument("--feature_transform", choices=["clr", "none"], default="clr")
    parser.add_argument("--log_transform", choices=["log1p", "log10", "none"], default="log1p")
    parser.add_argument("--output", default="predicted_tpm.tsv")
    parser.add_argument("--optimize", action="store_true")
    parser.add_argument("--opt_fasta", default="optimized.fasta")
    parser.add_argument("--max_iters", type=int, default=100)
    args = parser.parse_args()

    # --- モデル読み込み ---
    models = []
    model_names = []
    for path in args.model_pickles:
        name = os.path.basename(path).replace(".pkl", "").replace(".pickle", "")
        with open(path, 'rb') as f:
            models.append(pickle.load(f))
        model_names.append(name)

    optimized_records = []

    # --- カラム定義（TSV用） ---
    # 最適化がある場合は Before/After、ない場合は現在の値のみ
    if args.optimize:
        tsv_cols = ["gene_id", "status", "min_tpm_before", "min_tpm_after"]
        for name in model_names:
            tsv_cols.extend([f"{name}_before", f"{name}_after"])
    else:
        tsv_cols = ["gene_id", "status", "min_tpm"] + model_names

    # 標準出力（stdout）にヘッダー書き出し
    print("\t".join(tsv_cols), file=sys.stdout)

    # --- ヘッダー作成（stderr: 人間用） ---
    header_models_err = " | ".join([f"{name:<12}" for name in model_names])
    err_header = f"{'Gene_ID':<20} | {'Status':<10} | {'Min_TPM':<12} | {header_models_err}"
    print(err_header, file=sys.stderr)
    print("-" * len(err_header), file=sys.stderr)

    # --- メインループ ---
    for record in SeqIO.parse(args.fasta, "fasta"):
        gene_id = record.id
        seq = str(record.seq)

        # 1. 現状の予測（最適化前の基準値）
        initial_counts = compute_counts(seq)
        initial_clr = clr_from_counts_fast(initial_counts)
        before_preds = [predict_from_clr_fast(initial_clr, m, args.log_transform) for m in models]
        before_min = min(before_preds)
        
        # stderr: Original 表示
        vals_before_str = " | ".join([f"{v:>12.4f}" for v in before_preds])
        print(f"{gene_id[:20]:<20} | Original   | {before_min:>12.4f} | {vals_before_str}", file=sys.stderr)

        if args.optimize:
            # 2. 最適化の実行
            if len(models) == 1:
                # 単一モデルの場合
                opt_seq, after_min = optimize_single_model_fast(
                    gene_id, seq, models[0], args.log_transform, max_iters=args.max_iters
                )
                after_preds = [after_min]
            else:
                # 複数モデルの場合
                opt_seq, after_min, after_preds = optimize_multi_model_fast(
                    gene_id, seq, models, args.log_transform, max_iters=args.max_iters
                )
            
            # stderr: Optimized 表示
            vals_after_str = " | ".join([f"{v:>12.4f}" for v in after_preds])
            print(f"{gene_id[:20]:<20} | Optimized  | {after_min:>12.4f} | {vals_after_str}", file=sys.stderr)
            
            # FASTA保存用
            optimized_records.append(SeqRecord(Seq(opt_seq), id=gene_id, description="optimized"))
            
            # stdout: TSV行 (Before/After)
            row_vals = [gene_id, "optimized", f"{before_min:.6f}", f"{after_min:.6f}"]
            for b, a in zip(before_preds, after_preds):
                row_vals.extend([f"{b:.6f}", f"{a:.6f}"])
            print("\t".join(row_vals), file=sys.stdout)
            
        else:
            # 予測のみの場合 (stdout: TSV行)
            row_vals = [gene_id, "original", f"{before_min:.6f}"] + [f"{v:.6f}" for v in before_preds]
            print("\t".join(row_vals), file=sys.stdout)

    # 最終処理
    if args.optimize and optimized_records:
        SeqIO.write(optimized_records, args.opt_fasta, "fasta")
        print(f"\n[Done] Optimized sequences: {args.opt_fasta}", file=sys.stderr)

if __name__ == "__main__":
    main()