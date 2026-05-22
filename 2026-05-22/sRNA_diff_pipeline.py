#!/usr/bin/env python3
"""
sRNA_diff_pipeline.py
=====================
sRNA 差异表达分析一体化流水线

功能:
  1. 读取 all_sRNA.bed，过滤 sample_count >= N 的 sRNA
  2. 根据坐标回到各样本 ShortStack Results.txt，构建 count 矩阵
  3. 自动调用 Rscript 执行 DESeq2 差异分析 + 出图

用法 (直接运行或 sbatch):
  python3 sRNA_diff_pipeline.py [选项]

  # 或通过 sbatch:
  sbatch --wrap="python3 /path/to/sRNA_diff_pipeline.py" -p com300 -N1 -n1 \
         --cpus-per-task=4 --mem=16G -J sRNA_diff

选项:
  --workdir   工作目录 (含各样本子目录)
  --bed       all_sRNA.bed 路径
  --list      样本列表文件
  --outdir    输出目录
  --min-samples  sRNA 最少出现样本数 (默认 2)
  --rscript   DESeq2 R 脚本路径 (默认: 同目录下 sRNA_deseq2_analysis.R)
  --skip-r    仅构建 count 矩阵，不运行 R

作者: Kiro
"""

import re
import sys
import argparse
import subprocess
import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict

TAB = '\t'


# =============================================================================
# 参数解析
# =============================================================================
def parse_args():
    p = argparse.ArgumentParser(
        description='sRNA 差异表达分析流水线: count 矩阵构建 + DESeq2',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--workdir', default='/public/home/h14166/fang/Heyufei/smrna_seq',
                   help='工作目录 (含各样本子目录)')
    p.add_argument('--bed', default=None,
                   help='all_sRNA.bed 路径 (默认: $WORKDIR/all_sRNA.bed)')
    p.add_argument('--list', default=None,
                   help='样本列表文件 (默认: $WORKDIR/list)')
    p.add_argument('--outdir', default=None,
                   help='输出目录 (默认: $WORKDIR/sRNA_diff_analysis)')
    p.add_argument('--min-samples', type=int, default=2,
                   help='sRNA 最少出现样本数 (默认: 2)')
    p.add_argument('--rscript', default=None,
                   help='DESeq2 R 脚本路径 (默认: 同目录下 sRNA_deseq2_analysis.R)')
    p.add_argument('--skip-r', action='store_true',
                   help='仅构建 count 矩阵，跳过 DESeq2')
    return p.parse_args()


# =============================================================================
# 辅助函数
# =============================================================================
def classify_sample(name):
    """样本分组: WT / RDM_Xbp_KO / RDM_OEX"""
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
    读取 all_sRNA.bed，过滤 sample_count >= min_samples。
    返回 DataFrame: chrom, start, end, sRNA_name, sRNA_type, sample_count
    """
    records = []
    with open(bed_path) as fh:
        header_line = fh.readline()
        cols = header_line.lstrip('#').rstrip().split(TAB)
        col_idx = {c: i for i, c in enumerate(cols)}

        required = ['chrom', 'start', 'end', 'sRNA_name', 'sRNA_type', 'sample_count']
        for c in required:
            if c not in col_idx:
                sys.exit(f'[ERROR] bed 文件缺少列: {c}\n  实际列: {cols}')

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

    return pd.DataFrame(records)


def load_shortstack_results(results_path):
    """
    读取 ShortStack Results.txt -> dict{(chrom, start_0based, end): reads}
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
            start = int(m.group(2)) - 1  # 1-based -> 0-based
            end = int(m.group(3))
            try:
                reads = int(float(cols[i_reads]))
            except ValueError:
                reads = 0
            locus_reads[(chrom, start, end)] = reads

    return locus_reads


def build_count_matrix(srna_df, samples, workdir):
    """
    对每个 sRNA 区间，在各样本的 ShortStack Results 中找重叠 locus，
    将重叠 locus 的 Reads 求和作为该 sRNA 在该样本的表达量。
    """
    # 加载所有样本数据
    print(f'[INFO] 加载 {len(samples)} 个样本的 ShortStack Results ...')
    sample_data = {}
    for s in samples:
        rpath = Path(workdir) / s / 'shortstack_out' / 'Results.txt'
        if not rpath.exists():
            print(f'  [WARNING] 缺少: {rpath}')
            sample_data[s] = {}
            continue
        sample_data[s] = load_shortstack_results(rpath)
        print(f'  {s}: {len(sample_data[s])} loci')

    # 按 chrom 排序建索引，加速重叠查找
    print(f'[INFO] 构建 count 矩阵: {len(srna_df)} sRNAs x {len(samples)} samples ...')
    sample_index = {}
    for s in samples:
        idx = defaultdict(list)
        for (chrom, start, end), reads in sample_data[s].items():
            idx[chrom].append((start, end, reads))
        for chrom in idx:
            idx[chrom].sort()
        sample_index[s] = idx

    # 填充矩阵
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
                if le <= q_start:
                    continue
                if ls >= q_end:
                    break  # 已排序
                total_reads += lr
            mat.at[srna_id, s] = total_reads

    return mat


# =============================================================================
# 主流程
# =============================================================================
def main():
    args = parse_args()

    workdir = args.workdir
    bed_path = args.bed or f'{workdir}/all_sRNA.bed'
    list_file = args.list or f'{workdir}/list'
    outdir = Path(args.outdir or f'{workdir}/sRNA_diff_analysis')
    min_samples = args.min_samples
    script_dir = Path(__file__).resolve().parent
    rscript_path = args.rscript or str(script_dir / 'sRNA_deseq2_analysis.R')

    print('=' * 60)
    print('  sRNA 差异表达分析流水线')
    print('=' * 60)
    print(f'  工作目录:    {workdir}')
    print(f'  BED 文件:    {bed_path}')
    print(f'  样本列表:    {list_file}')
    print(f'  输出目录:    {outdir}')
    print(f'  最小样本数:  {min_samples}')
    print(f'  R 脚本:      {rscript_path}')
    print('=' * 60 + '\n')

    # 检查输入
    if not Path(bed_path).exists():
        sys.exit(f'[ERROR] BED 文件不存在: {bed_path}')
    if not Path(list_file).exists():
        sys.exit(f'[ERROR] 样本列表不存在: {list_file}')

    # 创建输出目录
    for sub in ('counts', 'results', 'plots'):
        (outdir / sub).mkdir(exist_ok=True, parents=True)

    # ===== Step 1: 读取样本 =====
    samples = [ln.strip() for ln in open(list_file) if ln.strip() and not ln.startswith('#')]
    print(f'[Step 1] 样本数: {len(samples)}')

    # ===== Step 2: 加载并过滤 BED =====
    print(f'\n[Step 2] 读取并过滤 BED (sample_count >= {min_samples}) ...')
    srna_df = load_bed(bed_path, min_samples)
    print(f'  保留 sRNA 数: {len(srna_df)}')

    if len(srna_df) == 0:
        sys.exit('[ERROR] 过滤后无 sRNA，请降低 --min-samples')

    # 类型分布
    print('\n  sRNA 类型分布:')
    for t, c in srna_df['sRNA_type'].value_counts().items():
        print(f'    {t}: {c}')

    # ===== Step 3: 构建 count 矩阵 =====
    print(f'\n[Step 3] 构建表达 count 矩阵 ...')
    mat = build_count_matrix(srna_df, samples, workdir)

    # 去除全零行
    nonzero_mask = mat.sum(axis=1) > 0
    n_before = len(mat)
    mat = mat[nonzero_mask]
    srna_df_filtered = srna_df[nonzero_mask.values].reset_index(drop=True)
    print(f'\n  去除全零行: {n_before} -> {len(mat)}')

    if len(mat) == 0:
        sys.exit('[ERROR] count 矩阵全为零，请检查坐标体系')

    # ===== Step 4: 保存输出 =====
    print(f'\n[Step 4] 保存结果 ...')

    # counts
    counts_file = outdir / 'counts' / 'sRNA_counts.csv'
    mat.to_csv(counts_file)
    print(f'  Counts:     {counts_file}')
    print(f'              {mat.shape[0]} sRNAs x {mat.shape[1]} samples, total={int(mat.values.sum())}')

    # annotation
    anno_file = outdir / 'counts' / 'sRNA_annotation.csv'
    srna_df_filtered.to_csv(anno_file, index=False)
    print(f'  Annotation: {anno_file}')

    # metadata
    meta = pd.DataFrame({
        'sample': samples,
        'group': [classify_sample(s) for s in samples]
    })
    meta.set_index('sample', inplace=True)
    meta_file = outdir / 'counts' / 'sample_metadata.csv'
    meta.to_csv(meta_file)
    print(f'  Metadata:   {meta_file}')
    print(f'\n  样本分组:')
    for g, cnt in meta.groupby('group').size().items():
        print(f'    {g}: {cnt}')

    # RPM
    lib_size = mat.sum(axis=0).replace(0, np.nan)
    rpm = mat.div(lib_size, axis=1) * 1e6
    rpm = rpm.fillna(0)
    rpm_file = outdir / 'counts' / 'sRNA_rpm.csv'
    rpm.to_csv(rpm_file)
    print(f'  RPM:        {rpm_file}')

    # ===== Step 5: 调用 DESeq2 =====
    if args.skip_r:
        print('\n[Step 5] 跳过 DESeq2 (--skip-r)')
        print(f'\n手动运行: Rscript {rscript_path} {outdir}')
    else:
        print(f'\n[Step 5] 运行 DESeq2 差异分析 ...')
        if not Path(rscript_path).exists():
            print(f'  [WARNING] R 脚本不存在: {rscript_path}')
            print(f'  请手动运行: Rscript sRNA_deseq2_analysis.R {outdir}')
        else:
            cmd = ['Rscript', rscript_path, str(outdir)]
            print(f'  命令: {" ".join(cmd)}')
            ret = subprocess.run(cmd)
            if ret.returncode != 0:
                print(f'  [ERROR] Rscript 返回非零: {ret.returncode}', file=sys.stderr)
                sys.exit(ret.returncode)

    # ===== 完成 =====
    print('\n' + '=' * 60)
    print('  流水线完成!')
    print('=' * 60)
    print(f'  Counts:  {outdir}/counts/')
    print(f'  Results: {outdir}/results/')
    print(f'  Plots:   {outdir}/plots/')
    print('=' * 60)


if __name__ == '__main__':
    main()
