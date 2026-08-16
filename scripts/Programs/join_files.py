#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import argparse
import pandas as pd

def read_with_prefix(path, sep, use_prefix=True):
    """ファイルを読み込み、列名にファイル名プレフィックスを付けるかどうかを制御"""
    prefix = os.path.splitext(os.path.basename(path))[0]
    df = pd.read_csv(path, sep=sep, header=0)
    cols = df.columns.tolist()
    # 第一列は必ず "KEY" に統一
    cols[0] = "KEY"
    if use_prefix:
        newcols = [cols[0]] + [f"{prefix}_{c}" for c in cols[1:]]
        df.columns = newcols
    else:
        df.columns = cols
    return df

def rename_duplicates(df):
    """列名が重複している場合に自動で .1, .2 を付ける"""
    cols = df.columns.tolist()
    seen = {}
    newcols = []
    for c in cols:
        if c not in seen:
            seen[c] = 1
            newcols.append(c)
        else:
            newcols.append(f"{c}.{seen[c]}")
            seen[c] += 1
    df.columns = newcols
    return df

def main():
    parser = argparse.ArgumentParser(
        description="Join multiple files by first column."
    )
    parser.add_argument("files", nargs="+", help="Input files")
    parser.add_argument("--outer", action="store_true", help="Use outer join (default)")
    parser.add_argument("--inner", action="store_true", help="Use inner join")
    parser.add_argument("--na", default="NA", help="Missing value placeholder")
    parser.add_argument("--sep", default="\t", help="Field separator for input/output (default: tab)")
    parser.add_argument("--no-prefix", action="store_true",
                        help="Do not add input file name as prefix to column headers")
    args = parser.parse_args()

    if args.outer and args.inner:
        print("Error: specify either --outer or --inner, not both.")
        sys.exit(1)

    how = "outer" if args.outer else "inner" if args.inner else "outer"

    # 最初のファイルを読み込む
    df = read_with_prefix(args.files[0], args.sep, use_prefix=not args.no_prefix)
    keycol = "KEY"
    merged = df

    # 残りのファイルを順次結合
    for f in args.files[1:]:
        df_next = read_with_prefix(f, args.sep, use_prefix=not args.no_prefix)
        merged = pd.merge(merged, df_next, on=keycol, how=how)
        merged = rename_duplicates(merged)

    # 欠測値を指定文字で埋める
    merged = merged.fillna(args.na)

    # 出力
    merged.to_csv(sys.stdout, sep=args.sep, index=False)

if __name__ == "__main__":
    main()