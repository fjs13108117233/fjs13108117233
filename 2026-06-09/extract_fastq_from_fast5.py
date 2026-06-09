#!/usr/bin/env python3
"""
Extract per-fast5 fastq directly from basecalls already stored inside the fast5.

Many MinKNOW/Guppy runs write the basecalls into the fast5 (Basecall_1D_000).
esox/scripts/basecall.py needs, for every <name>.fast5, a matching <name>.fastq
with those basecalls. This produces exactly that, so you do NOT need to re-run
Guppy if your fast5 already contain basecalls.

Usage:
    python3 scripts/extract_fastq_from_fast5.py \
        --fast5-path  fast5/OE1 \
        --output-path matched_fastq/OE1
"""
import os
import sys
import argparse

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from esox.fast_io import read_fast5


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fast5-path", required=True, help="fast5 file or dir of fast5 (one sample)")
    parser.add_argument("--output-path", required=True, help="dir to write per-fast5 .fastq")
    args = parser.parse_args()

    os.makedirs(args.output_path, exist_ok=True)

    if os.path.isfile(args.fast5_path):
        fast5_files = [args.fast5_path]
    else:
        fast5_files = [os.path.join(args.fast5_path, f)
                       for f in sorted(os.listdir(args.fast5_path)) if f.endswith('.fast5')]

    total, missing = 0, 0
    for f5 in fast5_files:
        out_fq = os.path.join(args.output_path, os.path.basename(f5).replace('.fast5', '.fastq'))
        reads = read_fast5(f5)
        written = 0
        with open(out_fq, 'w') as out:
            for read_id, rd in reads.items():
                if not rd.is_basecalled() or rd.basecalls is None or rd.basecalls == "":
                    missing += 1
                    continue
                out.write("@{}\n{}\n+\n{}\n".format(read_id, rd.basecalls, rd.phredq))
                written += 1
        total += written
        print("{} -> {} reads".format(os.path.basename(f5), written))

    print("Done. Wrote {} reads. {} reads had no basecalls inside the fast5.".format(total, missing))
    if total == 0:
        print("WARNING: no basecalls found in any fast5 -> your fast5 are raw signal only.")
        print("         You must run Guppy/Dorado first (see instructions).")
