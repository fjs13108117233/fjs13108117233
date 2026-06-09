#!/usr/bin/env python3
"""
Split a (possibly huge, 100GB+) Dorado/Guppy fastq into per-fast5 fastq files.

esox/scripts/basecall.py needs, for every <name>.fast5 in --fast5-path, a
matching <name>.fastq in --output-path containing exactly the reads of that
fast5. This routes one big fastq to those per-fast5 files in a SINGLE streaming
pass, so memory stays small (only a read_id -> fast5 index map is held).

Usage:
    python3 scripts/split_fastq_by_fast5.py \
        --fast5-path  /path/O8G/fast5/OE1 \
        --guppy-fastq /path/O8G/fq/OE1.fastq \
        --output-path /path/O8G/esox_out/OE1/matched_fastq \
        --max-open 512
"""
import os
import sys
import gzip
import argparse
from collections import OrderedDict

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from esox.fast_io import list_reads_ids


def _open_in(path):
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'r')


def collect_fastq_files(p):
    files = []
    if os.path.isfile(p):
        files.append(p)
    else:
        for f in sorted(os.listdir(p)):
            if f.endswith(('.fastq', '.fq', '.fastq.gz', '.fq.gz')):
                files.append(os.path.join(p, f))
    return files


class HandleCache:
    """LRU cache of output file handles so we never exceed the OS open-file limit."""
    def __init__(self, out_dir, basenames, max_open=512):
        self.out_dir = out_dir
        self.basenames = basenames                 # index -> output fastq basename
        self.max_open = max_open
        self.cache = OrderedDict()                  # index -> file handle
        self.truncated = set()                      # indices already created (write mode used once)

    def get(self, idx):
        h = self.cache.get(idx)
        if h is not None:
            self.cache.move_to_end(idx)
            return h
        mode = 'w' if idx not in self.truncated else 'a'
        path = os.path.join(self.out_dir, self.basenames[idx])
        h = open(path, mode)
        self.truncated.add(idx)
        self.cache[idx] = h
        if len(self.cache) > self.max_open:
            _, old = self.cache.popitem(last=False)
            old.close()
        return h

    def close_all(self):
        for h in self.cache.values():
            h.close()
        self.cache.clear()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fast5-path", required=True, help="Dir with .fast5 files (one sample)")
    parser.add_argument("--guppy-fastq", required=True, help="Big fastq file (or dir of fastq) for this sample")
    parser.add_argument("--output-path", required=True, help="Dir to write per-fast5 .fastq files")
    parser.add_argument("--max-open", type=int, default=512, help="Max simultaneously open output files")
    args = parser.parse_args()

    os.makedirs(args.output_path, exist_ok=True)

    fast5_files = [f for f in sorted(os.listdir(args.fast5_path)) if f.endswith('.fast5')]
    print("Found {} fast5 files".format(len(fast5_files)))

    # Phase 1: read_id -> fast5 index (small memory). Also keep output basenames.
    rid2idx = {}
    out_basenames = []
    dup = 0
    for idx, f5 in enumerate(fast5_files):
        out_basenames.append(f5.replace('.fast5', '.fastq'))
        for rid in list_reads_ids(os.path.join(args.fast5_path, f5)):
            if rid in rid2idx:
                dup += 1
            rid2idx[rid] = idx
    print("Indexed {} read ids from fast5 ({} duplicate ids)".format(len(rid2idx), dup))

    # ensure an (empty) output file exists for every fast5, even if it gets 0 reads
    for bn in out_basenames:
        open(os.path.join(args.output_path, bn), 'w').close()

    # Phase 2: stream the big fastq once, route each record to its fast5 file.
    cache = HandleCache(args.output_path, out_basenames, max_open=args.max_open)
    cache.truncated.update(range(len(out_basenames)))   # files already created above

    written = 0
    missing = 0
    seen = 0
    for fq in collect_fastq_files(args.guppy_fastq):
        print("Streaming {}".format(fq))
        with _open_in(fq) as fh:
            while True:
                header = fh.readline()
                if not header:
                    break
                seq = fh.readline()
                plus = fh.readline()
                qual = fh.readline()
                if not qual:
                    break
                seen += 1
                read_id = header.split(' ')[0][1:].strip('\n')
                idx = rid2idx.get(read_id)
                if idx is None:
                    missing += 1
                    continue
                out = cache.get(idx)
                out.write("@{}\n{}{}{}".format(read_id, seq, plus, qual))
                written += 1
                if seen % 1000000 == 0:
                    print("  ...{} reads streamed, {} routed".format(seen, written))
    cache.close_all()

    print("Done. Streamed {} fastq reads: {} routed to fast5, {} not in any fast5".format(
        seen, written, missing))
    if written == 0:
        print("WARNING: 0 reads matched. fastq read names probably don't match fast5 read_ids.")
