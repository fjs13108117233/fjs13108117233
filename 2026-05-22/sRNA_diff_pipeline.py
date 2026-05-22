#!/usr/bin/env python3
"""
sRNA_diff_pipeline.py
=====================
sRNA 差异表达分析一体化流水线（支持 SLURM array 并行）

功能:
  1. 读取 all_sRNA.bed，过滤 sample_count >= N 的 sRNA
  2. 根据坐标回到各样本 ShortStack Results.txt，构建 count 矩阵
  3. 自动调用 Rscript 执行 DESeq2 差异分析 + 出图

并行模式 (--array):
  将 sRNA 按 chunk 拆分，每个 SLURM array task 处理一块，
  最后由 --merge 步骤合并所有 chunk 的 count 矩阵。

用法:
  # 单机全量运行:
  python3 sRNA_diff_pipeline.py

  # 并行模式 (被 run_sRNA_diff.sh 自动调用):
  python3 sRNA_diff_pipeline.py --array --task-id 0 --n-tasks 10
  python3 sRNA_diff_pipeline.py --merge --n-tasks 10
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

    # 并行模式参数
    p.add_argument('--array', action='store_true',
                   help='并行模式: 只处理当前 chunk')
    p.add_argument('--task-id', type=int, default=None,
                   help='当前 array task ID (0-based)')
    p.add_argument('--n-tasks', type=int, default=10,
                   help='总 array task 数 (默认: 10)')
    p.add_argument('--merge', action='store_true',
                   help='合并所有 chunk 的结果')

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
    """读取 all_sRNA.bed，过滤 sample_count >= min_samples"""
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
                chrom = fields[col_idx['chrom']]
                start = int(fields[col_idx['start']])
                end = int(fields[col_idx['end']])
                # 用坐标作为唯一名称: chrom_start_end
                srna_name = f'{chrom}_{start}_{end}'
                records.append({
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'sRNA_name': srna_name,
                    'sRNA_type': fields[col_idx['sRNA_type']],
                    'sample_count': sc,
                })

    return pd.DataFrame(records)


def load_shortstack_results(results_path):
    """读取 ShortStack Results.txt -> dict{(chrom, start_0based, end): reads}"""
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
            start = int(m.group(2)) - 1
            end = int(m.group(3))
            try:
                reads = int(float(cols[i_reads]))
            except ValueError:
                reads = 0
            locus_reads[(chrom, start, end)] = reads

    return locus_reads


def build_count_matrix(srna_df, samples, workdir):
    """对每个 sRNA 区间，在各样本中找重叠 locus，Reads 求和"""
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

    # 按 chrom 排序建索引
    print(f'[INFO] 构建 count 矩阵: {len(srna_df)} sRNAs x {len(samples)} samples ...')
    sample_index = {}
    for s in samples:
        idx = defaultdict(list)
        for (chrom, start, end), reads in sample_data[s].items():
            idx[chrom].append((start, end, reads))
        for chrom in idx:
            idx[chrom].sort()
        sample_index[s] = idx

    mat = pd.DataFrame(0, index=srna_df['sRNA_name'].values, columns=samples, dtype=int)

    total = len(srna_df)
    for i, (_, row) in enumerate(srna_df.iterrows()):
        if (i + 1) % 5000 == 0:
            print(f'  进度: {i+1}/{total} ({100*(i+1)//total}%)')
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
                    break
                total_reads += lr
            mat.at[srna_id, s] = total_reads

    return mat


# =============================================================================
# 并行模式: 处理单个 chunk
# =============================================================================
def run_array_task(args):
    """array 模式: 只处理第 task_id 个 chunk"""
    workdir = args.workdir
    bed_path = args.bed or f'{workdir}/all_sRNA.bed'
    list_file = args.list or f'{workdir}/list'
    outdir = Path(args.outdir or f'{workdir}/sRNA_diff_analysis')
    min_samples = args.min_samples
    task_id = args.task_id
    n_tasks = args.n_tasks

    # 创建 tmp 目录
    tmp_dir = outdir / 'tmp_chunks'
    tmp_dir.mkdir(exist_ok=True, parents=True)

    # 读取样本
    samples = [ln.strip() for ln in open(list_file) if ln.strip() and not ln.startswith('#')]

    # 加载全部 sRNA
    srna_df = load_bed(bed_path, min_samples)
    total_srna = len(srna_df)

    # 计算当前 chunk 的范围
    chunk_size = (total_srna + n_tasks - 1) // n_tasks
    start_idx = task_id * chunk_size
    end_idx = min(start_idx + chunk_size, total_srna)

    if start_idx >= total_srna:
        print(f'[Task {task_id}] 无需处理 (start_idx={start_idx} >= total={total_srna})')
        # 写空文件
        pd.DataFrame(columns=samples).to_csv(tmp_dir / f'chunk_{task_id}.csv')
        return

    chunk_df = srna_df.iloc[start_idx:end_idx].reset_index(drop=True)
    print(f'[Task {task_id}/{n_tasks}] 处理 sRNA {start_idx}-{end_idx} (共 {len(chunk_df)} 条)')

    # 构建该 chunk 的 count
    mat = build_count_matrix(chunk_df, samples, workdir)

    # 保存 chunk
    out_file = tmp_dir / f'chunk_{task_id}.csv'
    mat.to_csv(out_file)
    print(f'[Task {task_id}] 完成，保存: {out_file}')


# =============================================================================
# 合并模式: 合并所有 chunk
# =============================================================================
def run_merge(args):
    """合并所有 chunk 的 count 矩阵，生成最终输出"""
    workdir = args.workdir
    bed_path = args.bed or f'{workdir}/all_sRNA.bed'
    list_file = args.list or f'{workdir}/list'
    outdir = Path(args.outdir or f'{workdir}/sRNA_diff_analysis')
    min_samples = args.min_samples
    n_tasks = args.n_tasks
    script_dir = Path(__file__).resolve().parent
    rscript_path = args.rscript or str(script_dir / 'sRNA_deseq2_analysis.R')

    tmp_dir = outdir / 'tmp_chunks'

    # 创建输出目录
    for sub in ('counts', 'results', 'plots'):
        (outdir / sub).mkdir(exist_ok=True, parents=True)

    samples = [ln.strip() for ln in open(list_file) if ln.strip() and not ln.startswith('#')]

    print(f'[MERGE] 合并 {n_tasks} 个 chunk ...')

    # 合并所有 chunk
    chunks = []
    for i in range(n_tasks):
        f = tmp_dir / f'chunk_{i}.csv'
        if not f.exists():
            print(f'  [WARNING] 缺少: {f}')
            continue
        df = pd.read_csv(f, index_col=0)
        if len(df) > 0:
            chunks.append(df)
        print(f'  chunk_{i}: {len(df)} sRNAs')

    if not chunks:
        sys.exit('[ERROR] 无有效 chunk')

    mat = pd.concat(chunks)
    mat = mat.astype(int)

    # 去除重复索引（chunk 边界可能重叠或 BED 有重复 sRNA_name）
    if mat.index.duplicated().any():
        n_dup = mat.index.duplicated().sum()
        print(f'  [WARNING] 发现 {n_dup} 个重复 sRNA_name，合并重复行 (取 max)')
        mat = mat.groupby(mat.index).max()

    print(f'[MERGE] 合并后: {len(mat)} sRNAs x {len(mat.columns)} samples')

    # 去除全零行
    nonzero_mask = mat.sum(axis=1) > 0
    mat = mat[nonzero_mask]
    print(f'[MERGE] 去除全零行后: {len(mat)}')

    if len(mat) == 0:
        sys.exit('[ERROR] count 矩阵全为零')

    # 保存 counts
    counts_file = outdir / 'counts' / 'sRNA_counts.csv'
    mat.to_csv(counts_file)
    print(f'[MERGE] Counts: {counts_file}')
    print(f'        {mat.shape[0]} sRNAs x {mat.shape[1]} samples, total={int(mat.values.sum())}')

    # 保存 annotation
    srna_df = load_bed(bed_path, min_samples)
    srna_df_filtered = srna_df[srna_df['sRNA_name'].isin(mat.index)].reset_index(drop=True)
    anno_file = outdir / 'counts' / 'sRNA_annotation.csv'
    srna_df_filtered.to_csv(anno_file, index=False)
    print(f'[MERGE] Annotation: {anno_file}')

    # 保存 metadata
    meta = pd.DataFrame({
        'sample': samples,
        'group': [classify_sample(s) for s in samples]
    })
    meta.set_index('sample', inplace=True)
    meta_file = outdir / 'counts' / 'sample_metadata.csv'
    meta.to_csv(meta_file)
    print(f'[MERGE] 样本分组:')
    for g, cnt in meta.groupby('group').size().items():
        print(f'    {g}: {cnt}')

    # RPM
    lib_size = mat.sum(axis=0).replace(0, np.nan)
    rpm = mat.div(lib_size, axis=1) * 1e6
    rpm = rpm.fillna(0)
    rpm.to_csv(outdir / 'counts' / 'sRNA_rpm.csv')

    # 调用 DESeq2
    if args.skip_r:
        print(f'\n[MERGE] 跳过 DESeq2 (--skip-r)')
        print(f'手动运行: Rscript {rscript_path} {outdir}')
    else:
        print(f'\n[MERGE] 运行 DESeq2 ...')
        if not Path(rscript_path).exists():
            print(f'  [WARNING] R 脚本不存在: {rscript_path}')
        else:
            cmd = ['Rscript', rscript_path, str(outdir)]
            print(f'  命令: {" ".join(cmd)}')
            ret = subprocess.run(cmd)
            if ret.returncode != 0:
                sys.exit(ret.returncode)

    # 清理 tmp
    import shutil
    shutil.rmtree(tmp_dir, ignore_errors=True)
    print(f'[MERGE] 已清理临时文件: {tmp_dir}')

    print('\n' + '=' * 60)
    print('  流水线完成!')
    print('=' * 60)
    print(f'  Counts:  {outdir}/counts/')
    print(f'  Results: {outdir}/results/')
    print(f'  Plots:   {outdir}/plots/')


# =============================================================================
# 单机模式（原始逻辑）
# =============================================================================
def run_single(args):
    """单机全量运行"""
    workdir = args.workdir
    bed_path = args.bed or f'{workdir}/all_sRNA.bed'
    list_file = args.list or f'{workdir}/list'
    outdir = Path(args.outdir or f'{workdir}/sRNA_diff_analysis')
    min_samples = args.min_samples
    script_dir = Path(__file__).resolve().parent
    rscript_path = args.rscript or str(script_dir / 'sRNA_deseq2_analysis.R')

    print('=' * 60)
    print('  sRNA 差异表达分析流水线 (单机模式)')
    print('=' * 60)
    print(f'  工作目录:    {workdir}')
    print(f'  BED 文件:    {bed_path}')
    print(f'  输出目录:    {outdir}')
    print(f'  最小样本数:  {min_samples}')
    print('=' * 60 + '\n')

    if not Path(bed_path).exists():
        sys.exit(f'[ERROR] BED 文件不存在: {bed_path}')
    if not Path(list_file).exists():
        sys.exit(f'[ERROR] 样本列表不存在: {list_file}')

    for sub in ('counts', 'results', 'plots'):
        (outdir / sub).mkdir(exist_ok=True, parents=True)

    samples = [ln.strip() for ln in open(list_file) if ln.strip() and not ln.startswith('#')]
    print(f'[Step 1] 样本数: {len(samples)}')

    print(f'\n[Step 2] 过滤 BED (sample_count >= {min_samples}) ...')
    srna_df = load_bed(bed_path, min_samples)
    print(f'  保留 sRNA 数: {len(srna_df)}')

    if len(srna_df) == 0:
        sys.exit('[ERROR] 过滤后无 sRNA')

    print('\n  sRNA 类型分布:')
    for t, c in srna_df['sRNA_type'].value_counts().items():
        print(f'    {t}: {c}')

    print(f'\n[Step 3] 构建 count 矩阵 ...')
    mat = build_count_matrix(srna_df, samples, workdir)

    nonzero_mask = mat.sum(axis=1) > 0
    n_before = len(mat)
    mat = mat[nonzero_mask]
    srna_df_filtered = srna_df[nonzero_mask.values].reset_index(drop=True)
    print(f'\n  去除全零行: {n_before} -> {len(mat)}')

    if len(mat) == 0:
        sys.exit('[ERROR] count 矩阵全为零')

    print(f'\n[Step 4] 保存结果 ...')
    mat.to_csv(outdir / 'counts' / 'sRNA_counts.csv')
    print(f'  Counts: {mat.shape[0]} x {mat.shape[1]}, total={int(mat.values.sum())}')

    srna_df_filtered.to_csv(outdir / 'counts' / 'sRNA_annotation.csv', index=False)

    meta = pd.DataFrame({'sample': samples, 'group': [classify_sample(s) for s in samples]})
    meta.set_index('sample', inplace=True)
    meta.to_csv(outdir / 'counts' / 'sample_metadata.csv')
    print(f'  样本分组:')
    for g, cnt in meta.groupby('group').size().items():
        print(f'    {g}: {cnt}')

    lib_size = mat.sum(axis=0).replace(0, np.nan)
    rpm = mat.div(lib_size, axis=1) * 1e6
    rpm.fillna(0).to_csv(outdir / 'counts' / 'sRNA_rpm.csv')

    if args.skip_r:
        print(f'\n[Step 5] 跳过 DESeq2 (--skip-r)')
        print(f'手动运行: Rscript {rscript_path} {outdir}')
    else:
        print(f'\n[Step 5] 运行 DESeq2 ...')
        if not Path(rscript_path).exists():
            print(f'  [WARNING] R 脚本不存在: {rscript_path}')
        else:
            cmd = ['Rscript', rscript_path, str(outdir)]
            print(f'  命令: {" ".join(cmd)}')
            ret = subprocess.run(cmd)
            if ret.returncode != 0:
                sys.exit(ret.returncode)

    print('\n' + '=' * 60)
    print('  流水线完成!')
    print('=' * 60)


# =============================================================================
# 入口
# =============================================================================
def main():
    args = parse_args()

    if args.merge:
        run_merge(args)
    elif args.array:
        if args.task_id is None:
            sys.exit('[ERROR] --array 模式需要 --task-id 参数')
        run_array_task(args)
    else:
        run_single(args)


if __name__ == '__main__':
    main()
