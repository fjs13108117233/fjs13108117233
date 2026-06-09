#!/usr/bin/env python3
"""
Aggregate esox modcall output into a per-sample 8-oxo-dG (O8G) summary.

modcall.py writes one .txt per fast5 with columns:
    read_id   basecalls_pos   oxog_score   5mer

This script:
  - pools all .txt for a sample
  - calls a G as 8-oxo-dG when oxog_score >= threshold
  - reports total Gs evaluated, O8G calls, and the O8G:G rate
  - writes a per-5mer breakdown (useful to cross-check static/kmer_performance.txt)

Usage:
    python3 scripts/summarize.py \
        --modcall-path esox_out/OE1/modcall_out \
        --threshold 0.95 \
        --sample OE1 \
        --output esox_out/summary.tsv
"""
import os
import glob
import argparse

import pandas as pd


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--modcall-path", required=True, help="Dir with modcall .txt files for one sample")
    parser.add_argument("--threshold", type=float, default=0.95, help="oxog_score cutoff to call O8G")
    parser.add_argument("--sample", required=True, help="Sample name (for the summary row)")
    parser.add_argument("--output", required=True, help="Summary TSV (appended to)")
    args = parser.parse_args()

    txt_files = sorted(glob.glob(os.path.join(args.modcall_path, "*.txt")))
    if not txt_files:
        raise SystemExit("No modcall .txt files found in {}".format(args.modcall_path))

    frames = [pd.read_csv(f, sep="\t") for f in txt_files]
    df = pd.concat(frames, ignore_index=True)

    df["is_o8g"] = df["oxog_score"] >= args.threshold

    total_g = len(df)
    n_o8g = int(df["is_o8g"].sum())
    n_reads = df["read_id"].nunique()
    rate = (n_o8g / total_g) if total_g else 0.0

    # per-sample row appended to the shared summary
    header = not os.path.isfile(args.output)
    with open(args.output, "a") as out:
        if header:
            out.write("sample\tn_reads\ttotal_G_evaluated\tn_O8G_calls\tO8G_per_G_rate\tthreshold\n")
        out.write("{}\t{}\t{}\t{}\t{:.3e}\t{}\n".format(
            args.sample, n_reads, total_g, n_o8g, rate, args.threshold))

    # per-5mer breakdown next to the modcall output
    per5 = (df.groupby("5mer")
              .agg(n_G=("oxog_score", "size"), n_O8G=("is_o8g", "sum"))
              .reset_index())
    per5["O8G_rate"] = per5["n_O8G"] / per5["n_G"]
    per5 = per5.sort_values("n_O8G", ascending=False)
    per5_out = os.path.join(args.modcall_path, "..", "{}_per5mer.tsv".format(args.sample))
    per5.to_csv(per5_out, sep="\t", index=False)

    print("[{}] reads={}  G_evaluated={}  O8G_calls={}  rate={:.3e}  (thr={})".format(
        args.sample, n_reads, total_g, n_o8g, rate, args.threshold))
    print("Per-5mer breakdown -> {}".format(os.path.normpath(per5_out)))
