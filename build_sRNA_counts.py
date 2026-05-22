#!/usr/bin/env python3
"""
build_sRNA_counts.py
====================
从 all_sRNA.bed 构建 sRNA 表达 count 矩阵，用于 DESeq2 差异分析。

逻辑：
  1. 读取 all_sRNA.bed，过滤 sample_count >= MIN_SAMPLE_COUNT 的 sRNA
  2. 对每个保留的 sRNA，根据其 chrom:start-end 坐标，
     回到各样本的 ShortStack Results.txt 提取 Reads 数
  3. 输出: counts 矩阵 (sRNA x samples) + sample metadata

用法:
  python3 build_sRNA_counts.py [--min-samples 2] [--bed /path/to/all_sRNA.bed]
"""

import re
import sys
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict

TAB = '\t'


def parse_args():
    p = argparse.ArgumentParser(description='Build sRNA count matrix from all_sRNA.bed + ShortStack Results')
    p.add_argument('--bed', default=None,
                   help='Path to all_sRNA.bed (default: $WORKDIR/all_sRNA.bed)')
    p.add_argument('--workdir', default='/public/home/h14166/fang/Heyufei/smrna_seq',
                   help='Working directory')
    p.add_argument('--list', default=None,
                   help='Sample list file (default: $WORKDIR/list)')
    p.add_argument('--outdir', default=None,
                   help='Output directory (default: $WORKDIR/sRNA_diff_analysis)')
    p.add_argument('--min-samples', type=int, default=2,
                   help='Minimum number of samples a sRNA must appear in (default: 2)')
    return p.parse_args()


def classify_sample(name):
    """样本分组：WT / RDM_KO / RDM_OE"""
    if re.match(r'^WT', name, re.I):
        return 'WT'
    m = re.match(r'^RDM_-?([0-9]+)bp', name)
    if m:
        return 'RDM_' + m.group(1) + 'bp_KO'
    m = re.match(r'^RDM_(OE[0-9]+)', name)
    if m:
        return 'RDM_' + m.group(1)
    return 'Unknown'


def load_bed(bed_path, min_samples):
    """
    读取 all_sRNA.bed，过滤 sample_count >= min_samples 的 sRNA。
    返回 DataFrame: chrom, start, end, sRNA_name, sRNA_type, sample_count
    """
    records = []
    with open(bed_path) as fh:
        header_line = fh.readline()
        # 解析表头
        cols = header_line.lstrip('#').rstrip().split(TAB)
        col_idx = {c: i for i, c in enumerate(cols)}

        # 必需列
        required = ['chrom', 'start', 'end', 'sRNA_name', 'sRNA_type', 'sample_count']
        for c in required:
            if c not in col_idx:
                sys.exit(f'[ERROR] bed 文件缺少列: {c}, 实际列: {cols}')

        for line in fh:
            if not line.strip():
                continue
            fields = line.rstrip().split(TAB)
            try:
                sc = int(fields[col_idx['sample_count']])
            except (ValueError, IndexError):
                continue
            if sc >= min_samples:
                records.append({
                    'chrom': fields[col_idx['chrom']],
                    'start': int(fields[col_idx['start']]),
                    'end': int(fields[col_idx['end']]),
                    'sRNA_name': fields[col_idx['sRNA_name']],
                    'sRNA_type': fields[col_idx['sRNA_type']],
                    'sample_count': sc,
                })

    df = pd.DataFrame(records)
    return df


def load_shortstack_results(results_path):
    """
    读取 ShortStack Results.txt，返回 dict:
      key = (chrom, start_0based, end) -> reads count
    坐标转换: Results.txt 中 Locus 是 1-based，转为 0-based start
    """
    locus_reads = {}
    with open(results_path) as fh:
        header = fh.readline().rstrip().split(TAB)
        header = [h.lstrip('#') for h in header]
        try:
            i_locus = header.index('Locus')
            i_reads = header.index('Reads')
        except ValueError:
            return locus_reads

        for line in fh:
            cols = line.rstrip().split(TAB)
            if len(cols) <= max(i_locus, i_reads):
                continue
            loc = cols[i_locus]
            m = re.match(r'(\S+):(\d+)-(\d+)', loc)
            if not m:
                continue
            chrom = m.group(1)
            start = int(m.group(2)) - 1  # 转为 0-based
            end = int(m.group(3))
            try:
                reads = int(float(cols[i_reads]))
            except ValueError:
                reads = 0
            locus_reads[(chrom, start, end)] = reads

    return locus_reads


def build_count_matrix(srna_df, samples, workdir):
    """
    对于每个 sRNA (chrom, start, end)，在每个样本的 ShortStack Results.txt 中
    查找坐标重叠的 locus，取 reads 之和。

    重叠定义: sRNA 区间与 ShortStack locus 有任意碱基交集
    """
    # 预加载所有样本的 locus -> reads
    print(f'[INFO] 加载 {len(samples)} 个样本的 ShortStack Results ...')
    sample_data = {}
    for s in samples:
        rpath = Path(workdir) / s / 'shortstack_out' / 'Results.txt'
        if not rpath.exists():
            print(f'  [WARNING] 缺少: {rpath}')
            sample_data[s] = {}
            continue
        sample_data[s] = load_shortstack_results(rpath)

    # 构建 count 矩阵
    print(f'[INFO] 构建 count 矩阵: {len(srna_df)} sRNAs x {len(samples)} samples ...')

    # 为加速，按 chrom 建立索引
    # sample_index[sample][chrom] = list of (start, end, reads)
    sample_index = {}
    for s in samples:
        idx = defaultdict(list)
        for (chrom, start, end), reads in sample_data[s].items():
            idx[chrom].append((start, end, reads))
        # 排序以便后续二分
        for chrom in idx:
            idx[chrom].sort()
        sample_index[s] = idx

    mat = pd.DataFrame(0, index=srna_df['sRNA_name'].values, columns=samples, dtype=int)

    for _, row in srna_df.iterrows():
        q_chrom = row['chrom']
        q_start = row['start']
        q_end = row['end']
        srna_id = row['sRNA_name']

        for s in samples:
            loci = sample_index[s].get(q_chrom, [])
            total_reads = 0
            for (ls, le, lr) in loci:
                # 重叠判断: not (le <= q_start or ls >= q_end)
                if le <= q_start:
                    continue
                if ls >= q_end:
                    break  # 已排序，后面更不会重叠
                total_reads += lr
            mat.at[srna_id, s] = total_reads

    return mat


def main():
    args = parse_args()

    workdir = args.workdir
    bed_path = args.bed or (workdir + '/all_sRNA.bed')
    list_file = args.list or (workdir + '/list')
    outdir = Path(args.outdir or (workdir + '/sRNA_diff_analysis'))
    min_samples = args.min_samples

    # 创建输出目录
    for sub in ('counts', 'results', 'plots'):
        (outdir / sub).mkdir(exist_ok=True, parents=True)

    # 读取样本列表
    samples = [ln.strip() for ln in open(list_file) if ln.strip() and not ln.startswith('#')]
    print(f'[INFO] 样本数: {len(samples)}')

    # 1. 加载并过滤 bed
    print(f'[INFO] 读取 bed 文件: {bed_path}')
    srna_df = load_bed(bed_path, min_samples)
    print(f'[INFO] 过滤后 sRNA 数 (sample_count >= {min_samples}): {len(srna_df)}')

    if len(srna_df) == 0:
        sys.exit('[ERROR] 过滤后无 sRNA，请检查 --min-samples 参数')

    # sRNA 类型分布
    print('\n[INFO] sRNA 类型分布:')
    type_counts = srna_df['sRNA_type'].value_counts()
    for t, c in type_counts.items():
        print(f'  {t}: {c}')

    # 2. 构建 count 矩阵
    mat = build_count_matrix(srna_df, samples, workdir)

    # 过滤全零行
    nonzero_mask = mat.sum(axis=1) > 0
    n_before = len(mat)
    mat = mat[nonzero_mask]
    srna_df_filtered = srna_df[nonzero_mask.values].reset_index(drop=True)
    print(f'\n[INFO] 去除全零行: {n_before} -> {len(mat)}')

    if len(mat) == 0:
        sys.exit('[ERROR] count 矩阵全为零')

    # 3. 保存 counts
    counts_file = outdir / 'counts' / 'sRNA_counts.csv'
    mat.to_csv(counts_file)
    print(f'[INFO] Counts 矩阵已保存: {counts_file}')
    print(f'       维度: {mat.shape[0]} sRNAs x {mat.shape[1]} samples')
    print(f'       总 reads: {int(mat.values.sum())}')

    # 4. 保存 sRNA 注释信息
    anno_file = outdir / 'counts' / 'sRNA_annotation.csv'
    srna_df_filtered.to_csv(anno_file, index=False)
    print(f'[INFO] sRNA 注释已保存: {anno_file}')

    # 5. 保存 sample metadata
    meta = pd.DataFrame({
        'sample': samples,
        'group': [classify_sample(s) for s in samples]
    })
    meta.set_index('sample', inplace=True)
    meta_file = outdir / 'counts' / 'sample_metadata.csv'
    meta.to_csv(meta_file)
    print(f'\n[INFO] 样本分组:')
    print(meta.groupby('group').size().to_string())
    print(f'[INFO] Metadata 已保存: {meta_file}')

    # 6. RPM 标准化（参考用）
    lib_size = mat.sum(axis=0).replace(0, np.nan)
    rpm = mat.div(lib_size, axis=1) * 1e6
    rpm = rpm.fillna(0)
    rpm.to_csv(outdir / 'counts' / 'sRNA_rpm.csv')
    print(f'[INFO] RPM 矩阵已保存: {outdir / "counts" / "sRNA_rpm.csv"}')

    print(f'\n[INFO] 完成！下一步: Rscript sRNA_deseq2_analysis.R {outdir}')


if __name__ == '__main__':
    main()
