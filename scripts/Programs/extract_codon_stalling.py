import pysam
import sys
import os
import argparse
import gffutils
from Bio import SeqIO

def get_args():
    parser = argparse.ArgumentParser(description="Ribo-seq P-site Codon Extractor (Variant Filtering Added)")
    parser.add_argument("config", help="Config file with offsets")
    parser.add_argument("annot_dir", help="Annotation directory")
    parser.add_argument("bam", help="Input BAM file")
    parser.add_argument("out_tsv", help="Output TSV file")
    parser.add_argument("-d", "--db", help="Path to gffutils database for representative filtering")
    parser.add_argument("--all-transcripts", action="store_true", help="Analyze all variants (Skip filtering)")
    parser.add_argument("--full", action="store_true", help="Output full transcript (including 5'UTR/3'UTR)")
    return parser.parse_args()

def get_representative_tids(db_path):
    if not db_path or not os.path.exists(db_path):
        return None
    
    db = gffutils.FeatureDB(db_path)
    selected_tids = set()
    print(f"[*] 代表バリアントを選別中 (ID正規化適用)...")
    
    for gene in db.features_of_type('gene'):
        transcripts = list(db.children(gene, featuretype=['transcript', 'mRNA']))
        if not transcripts: continue

        target_id = None
        for t in transcripts:
            attrs = str(t.attributes).lower()
            if any(tag in attrs for tag in ['mane_select', 'appris_principal_1']):
                target_id = t.id
                break
        
        if not target_id:
            max_len = -1
            for t in transcripts:
                try:
                    current_len = sum(len(c) for c in db.children(t, featuretype='CDS'))
                    if current_len > max_len:
                        max_len = current_len
                        target_id = t.id
                except: continue
        
        if not target_id: 
            target_id = transcripts[0].id
        
        # --- ID形式のズレを解消する正規化 (2026年推奨設定) ---
        # 1. 'transcript:' などの接頭辞を除去
        clean_id = target_id.split(':')[-1]
        # 2. バージョン番号 (.1, .2 など) を除去
        # これにより、DB側がバージョン付き、BAM側がバージョンなしでも一致するようになります
        base_id = clean_id.split('.')[0]
        
        selected_tids.add(base_id)
    
    return selected_tids

def main():
    args = get_args()
    
    # --- 【新規追加】代表バリアントのフィルタリング ---
    target_tids = None
    if not args.all_transcripts and args.db:
        target_tids = get_representative_tids(args.db)
        if target_tids:
            print(f"[*] {len(target_tids)} 個の代表トランスクリプトをフィルタとして登録しました。")
    # -----------------------------------------------

    # 1. configからオフセットを読み取る
    offsets = {}
    with open(args.config, 'r') as f:
        for line in f:
            if line.startswith("#") or not line.strip(): continue
            cols = line.strip().split("\t")
            if len(cols) >= 5:
                lens = cols[3].split(",")
                locs = cols[4].split(",")
                for length, loc in zip(lens, locs):
                    offsets[int(length)] = int(loc)

    # 2. CDS範囲情報を読み込む (1-based)
    cds_info = {}
    txt_path = os.path.join(args.annot_dir, "transcripts_cds.txt")
    with open(txt_path, "r") as f:
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) >= 3:
                tid = cols[0]
                # フィルタリング適用
                if target_tids is not None and tid not in target_tids:
                    continue
                start = int(cols[1]) - 1 
                end = int(cols[2]) - 1   
                cds_info[tid] = (start, end)

    # 3. 配列ロード
    fasta_path = os.path.join(args.annot_dir, "transcripts_sequence.fa")
    if not os.path.exists(fasta_path):
        fasta_path = os.path.join(args.annot_dir, "trans_sc.fasta")
    seqs = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
    
    # 4. BAMからP-site抽出とコドン集計
    codon_counts = {}

    with pysam.AlignmentFile(args.bam, "rb") as sam:
        for r in sam.fetch(until_eof=True):
            tid_full = r.reference_name
            if not tid_full:
                continue

            # BAM側のIDもドットで切り落として比較する (例: ENST0000.1 -> ENST0000)
            tid = tid_full.split('.')[0]

            # 代表バリアント制限または全バリアントの判定
            if (target_tids is None or tid in target_tids) and r.query_length in offsets:
                # フィルタリング済みのtidがcds_infoやseqsに存在するか確認
                # 注: cds_infoのキーも get_representative_tids 内で同様にドット削除済みであることが前提
                if tid in cds_info and tid in seqs:
                    cds_start, _ = cds_info[tid]
                    # Transcript上のP-site絶対座標
                    psite_abs = r.reference_start + offsets[r.query_length]
                    
                    # CDS開始点(0)を基準にした相対座標
                    psite_in_cds = psite_abs - cds_start
                    
                    # コドン枠への集約 (// 3 は負の数でも枠を維持)
                    codon_start_idx = (psite_in_cds // 3) * 3
                    
                    key = (tid, codon_start_idx)
                    codon_counts[key] = codon_counts.get(key, 0) + 1


    # 5. 書き出しフェーズ (既存のロジックを維持)
    with open(args.out_tsv, "w") as f:
        f.write("transcript_id\tcds_position\tcodon_seq\tcount\n")
        for (tid, cds_pos), c in sorted(codon_counts.items()):
            cds_start, cds_end = cds_info[tid]
            abs_pos = cds_start + cds_pos
            
            if abs_pos < 0 or abs_pos + 3 > len(seqs[tid].seq):
                continue

            if not args.full:
                if not (cds_start <= abs_pos <= cds_end - 2):
                    continue

            codon_seq = str(seqs[tid].seq[abs_pos:abs_pos+3]).upper()
            if len(codon_seq) == 3:
                f.write(f"{tid}\t{cds_pos + 1}\t{codon_seq}\t{c}\n")

if __name__ == "__main__":
    main()
