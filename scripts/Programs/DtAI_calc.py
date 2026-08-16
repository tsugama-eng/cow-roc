#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
stAI calculator with revised wobble rules (no inosine).
Usage:
    python stAI_calc.py cds.fasta trna_counts.txt codon_freqs.txt
"""

import sys
import math

# ---------- Utilities ----------
def dna_to_rna(seq: str) -> str:
    return seq.upper().replace("T", "U")

def expand_codons():
    bases = ["A","U","G","C"]
    return [a+b+c for a in bases for b in bases for c in bases]

def wc_pair(x,y):
    return (x=="A" and y=="U") or (x=="U" and y=="A") or \
           (x=="G" and y=="C") or (x=="C" and y=="G")

# Revised wobble map (anticodon 5' base vs codon 3' base)
REVISED_WOBBLE = {
    "A": {"U": "WC", "C": "WOBBLE"},
    "G": {"C": "WC", "U": "WOBBLE"},
    "U": {"A": "WC", "G": "WOBBLE"},
    "C": {"G": "WC", "A": "WOBBLE"}
}

def pair_kind_wobble(codon3, anticodon1):
    return REVISED_WOBBLE.get(anticodon1,{}).get(codon3,"")

def pair_kind(codon, anticodon):
    # codon: RNA 5'->3', anticodon: RNA 5'->3'
    c1,c2,c3 = codon
    a1,a2,a3 = anticodon
    if not wc_pair(c1,a3): return ""
    if not wc_pair(c2,a2): return ""
    return pair_kind_wobble(c3,a1)

# ---------- Wc and weights ----------
def compute_Wc(codons, anticodon_counts, alpha_wobble=0.4):
    W = {co:0.0 for co in codons}
    for ac,n in anticodon_counts.items():
        if n<=0 or len(ac)!=3: continue
        for co in codons:
            k = pair_kind(co,ac)
            if not k: continue
            pen = 0.0 if k=="WC" else alpha_wobble
            W[co] += n*(1.0-pen)
    return W

def normalize_weights(W):
    m = max(W.values())
    return {k:(v/m if m>0 else 0.0) for k,v in W.items()}

STOP_CODONS = {"UAA","UAG","UGA"}

def calc_DtAI(seq, weights):
    rna = seq
    L = len(rna)//3
    if L==0: return float("nan")
    eps=1e-9
    logs=[]
    # ストップコドンを除外
    last_codon = rna[-3:]
    if last_codon in STOP_CODONS:
        L -= 1
        rna = rna[:-3]
    if L == 0:
        return float("nan")
    eps = 1e-9
    logs = []
    for i in range(L):
        co=rna[3*i:3*i+3]
        w=weights.get(co,0.0)
        if w<=0: w=eps
        logs.append(math.log(w))
    return math.exp(sum(logs)/L) if L>0 else float("nan")

# ---------- File readers ----------
def read_trna_counts(path):
    pairs=[]
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split()
            if len(parts)!=2: continue
            ac, n = parts
            try:
                n = int(n)
            except ValueError:
                continue
            ac_rna = dna_to_rna(ac)
            pairs.append((ac_rna, n))
    anticodon_counts={}
    total = 0
    for ac,n in pairs:
        anticodon_counts[ac]=anticodon_counts.get(ac,0)+n
        total += n
    # 正規化して頻度に変換
    if total > 0:
        anticodon_freq = {ac: v/total for ac,v in anticodon_counts.items()}
    else:
        anticodon_freq = anticodon_counts
    return anticodon_freq

def read_fasta(path):
    seqs={}
    with open(path) as f:
        name=None
        seq=[]
        for line in f:
            line=line.strip()
            if not line: continue
            if line.startswith(">"):
                if name: seqs[name]=dna_to_rna("".join(seq))
                name=line[1:].split()[0]
                seq=[]
            else:
                seq.append(line)
        if name: seqs[name]=dna_to_rna("".join(seq))
    return seqs

def read_codon_freqs(path):
    freqs = {}
    total = 0
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) != 2:
                continue
            co, n = parts
            try:
                n = int(n)
            except ValueError:
                continue
            rna_co = dna_to_rna(co)
            freqs[rna_co] = freqs.get(rna_co, 0) + n
            total += n
    # 正規化して頻度に変換
    if total > 0:
        freqs = {co: v/total for co,v in freqs.items()}
    return freqs

# ---------- Effective supply ----------
def compute_effective_supply(anticodon_counts, codon_freqs, codons, alpha_wobble=0.4):
    """Return dict of anticodon -> effective supply number"""
    eff_supply = {}
    for ac, n in anticodon_counts.items():
        if n <= 0 or len(ac) != 3:
            continue
        denom = 0.0
        for co in codons:
            k = pair_kind(co, ac)
            if not k:
                continue
            if k == "WC":
                denom += codon_freqs.get(co, 0)
            elif k == "WOBBLE":
                denom += codon_freqs.get(co, 0) * alpha_wobble
        if denom > 0:
            eff_supply[ac] = n / denom
        else:
            eff_supply[ac] = 0.0
    return eff_supply

# ---------- Alpha estimation ----------
def estimate_alpha_wobble(codons, anticodon_counts):
    wc_total = 0
    wobble_total = 0
    for ac,n in anticodon_counts.items():
        if n<=0 or len(ac)!=3: continue
        for co in codons:
            k = pair_kind(co,ac)
            if k=="WC":
                wc_total += n
            elif k=="WOBBLE":
                wobble_total += n
    if wc_total+wobble_total == 0:
        return 0.4  # fallback
    return wobble_total/(wc_total+wobble_total)

# ---------- Main ----------
def main():
    if len(sys.argv)<4:
        print("Usage: python stai_calc.py cds.fasta trna_counts.txt codon_freqs.txt")
        sys.exit(1)
    fasta_file=sys.argv[1]
    trna_file=sys.argv[2]
    codon_file=sys.argv[3]

    cds=read_fasta(fasta_file)
    anticodon_counts=read_trna_counts(trna_file)
    codon_freqs=read_codon_freqs(codon_file)
    codons_all=expand_codons()

    # 自動推定
    alpha_wobble = estimate_alpha_wobble(codons_all, anticodon_counts)

    # 実サプライ数に修正
    eff_supply = compute_effective_supply(anticodon_counts, codon_freqs, codons_all, alpha_wobble)

    Wc=compute_Wc(codons_all, eff_supply, alpha_wobble)
    weights=normalize_weights(Wc)

    # Calculate stAI for each CDS
    for name,seq in cds.items():
        val=calc_DtAI(seq,weights)
        print(f"{name}\t{val:.6f}")

if __name__=="__main__":
    main()