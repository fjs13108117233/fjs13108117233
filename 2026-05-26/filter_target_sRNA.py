#!/usr/bin/env python3
"""
filter_target_sRNA.py - 筛选目标 sRNA + 邻近基因注释 (精简版)

筛选: OE1+OE2 均显著上调, KO6bp+KO8bp 均显著下调
输出: chrom, start, end, sRNA_type, {group}_RPM/pvalue/log2FC/state, 邻近基因注释
state: 0=上调, 1=下调, 3=不显著

用法: python3 filter_target_sRNA.py [--workdir ...] [--gtf ...] [--max-dist 10000]
"""

import argparse
import sys
import bisect
import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser(description='筛选目标 sRNA + 邻近基因注释')
    p.add_argument('--workdir', default='/public/home/h14166/fang/Heyufei/smrna_seq/sRNA_diff_analysis')
    p.add_argument('--gtf', default='/public/home/h14166/refer/maizev4.gtf')
    p.add_argument('--padj', type=float, default=0.05)
    p.add_argument('--lfc', type=float, default=1.0)
    p.add_argument('--max-dist', type=int, default=10000)
    p.add_argument('--output', default=None)
    return p.parse_args()


# =============================================================================
# 核心逻辑
# =============================================================================

# DE 文件名 -> (组名, RPM 列名匹配关键字)
GROUPS = {
    'OE1':   ('DE_RDM_OE1_vs_WT.csv',    'RDM_OE1'),
    'OE2':   ('DE_RDM_OE2_vs_WT.csv',    'RDM_OE2'),
    'KO6bp': ('DE_RDM_6bp_KO_vs_WT.csv', 'RDM_-6bp'),
    'KO8bp': ('DE_RDM_8bp_KO_vs_WT.csv', 'RDM_-8bp'),
}


def load_de(workdir):
    """加载 4 组 DE 结果，去重，返回 dict[group] -> DataFrame (index=sRNA_name)"""
    results_dir = Path(workdir) / 'results'
    de = {}
    for grp, (fname, _) in GROUPS.items():
        df = pd.read_csv(results_dir / fname, low_memory=False)
        df = df.drop_duplicates(subset='sRNA_name', keep='first').set_index('sRNA_name')
        de[grp] = df
        print(f'  {grp}: {len(df)} sRNAs')
    return de


def get_state(padj, lfc, padj_cut, lfc_cut):
    """0=上调, 1=下调, 3=不显著"""
    if pd.isna(padj):
        return 3
    if padj < padj_cut and lfc > lfc_cut:
        return 0
    if padj < padj_cut and lfc < -lfc_cut:
        return 1
    return 3


def filter_and_merge(de, workdir, padj_cut, lfc_cut):
    """筛选 + 合并为最终表格"""
    # 四组交集
    common = set(de['OE1'].index)
    for grp in ['OE2', 'KO6bp', 'KO8bp']:
        common &= set(de[grp].index)
    print(f'  四组共有: {len(common)}')

    # 加载 RPM
    rpm_df = pd.read_csv(Path(workdir) / 'counts' / 'sRNA_rpm.csv', index_col=0)
    # 样本列映射
    rpm_cols = {}
    for grp, (_, key) in GROUPS.items():
        rpm_cols[grp] = [c for c in rpm_df.columns if key in c]

    # 向量化筛选
    rows = []
    for srna in common:
        states = {}
        for grp in ['OE1', 'OE2', 'KO6bp', 'KO8bp']:
            r = de[grp].loc[srna]
            states[grp] = get_state(r['padj'], r['log2FoldChange'], padj_cut, lfc_cut)

        # OE 必须全上调(0), KO 必须全下调(1)
        if states['OE1'] != 0 or states['OE2'] != 0:
            continue
        if states['KO6bp'] != 1 or states['KO8bp'] != 1:
            continue

        # 从 sRNA_name 解析坐标: "chrom_start_end"
        parts = srna.rsplit('_', 2)
        chrom, start, end = parts[0], int(parts[1]), int(parts[2])

        row = {'chrom': chrom, 'start': start, 'end': end}
        # sRNA_type 从任一 DE 表中获取
        row['sRNA_type'] = de['OE1'].loc[srna].get('sRNA_type', '')

        for grp in ['OE1', 'OE2', 'KO6bp', 'KO8bp']:
            r = de[grp].loc[srna]
            # RPM: 该组样本平均
            if srna in rpm_df.index and rpm_cols[grp]:
                rpm_val = rpm_df.loc[srna, rpm_cols[grp]].mean()
            else:
                rpm_val = r['baseMean']
            row[f'{grp}_RPM'] = round(rpm_val, 4)
            row[f'{grp}_pvalue'] = r['pvalue']
            row[f'{grp}_log2FC'] = round(r['log2FoldChange'], 4)
            row[f'{grp}_state'] = states[grp]

        rows.append(row)

    df = pd.DataFrame(rows)
    if len(df) > 0:
        df = df.sort_values(['chrom', 'start']).reset_index(drop=True)
    print(f'  筛选结果: {len(df)} 个目标 sRNA')
    return df


# =============================================================================
# GTF 邻近基因注释
# =============================================================================

def parse_gtf(gtf_path):
    """解析 GTF -> dict[chrom] -> sorted list of (start, end, strand, gene_id)"""
    genes = defaultdict(list)
    n = 0
    with open(gtf_path) as f:
        for line in f:
            if line[0] == '#':
                continue
            cols = line.split('\t')
            if len(cols) < 9 or cols[2] != 'gene':
                continue
            chrom = cols[0]
            start, end, strand = int(cols[3]), int(cols[4]), cols[6]
            # 提取 gene_id
            attrs = cols[8]
            gid = ''
            for a in attrs.split(';'):
                a = a.strip()
                if a.startswith('gene_id'):
                    gid = a.split('"')[1] if '"' in a else a.split()[-1]
                    break
            genes[chrom].append((start, end, strand, gid))
            n += 1
    for ch in genes:
        genes[ch].sort()
    print(f'  GTF: {n} genes, {len(genes)} chroms')
    return genes


def nearest_gene(chrom, s_start, s_end, genes, max_dist):
    """找最近基因，返回 (gene_id, distance, direction)"""
    if chrom not in genes:
        return 'NA', 'NA', 'NA'
    gl = genes[chrom]
    starts = [g[0] for g in gl]
    idx = bisect.bisect_left(starts, s_start)

    best_d, best_g = float('inf'), None
    for i in range(max(0, idx - 10), min(len(gl), idx + 10)):
        gs, ge, gstrand, gid = gl[i]
        if s_end <= gs:
            d = gs - s_end
        elif s_start >= ge:
            d = s_start - ge
        else:
            d = 0
        if d < best_d:
            best_d, best_g = d, gl[i]

    if best_g is None or best_d > max_dist:
        return 'NA', 'NA', 'NA'

    gs, ge, gstrand, gid = best_g
    if best_d == 0:
        direction = 'overlap'
    elif s_end <= gs:
        direction = 'upstream' if gstrand == '+' else 'downstream'
    else:
        direction = 'downstream' if gstrand == '+' else 'upstream'

    return gid, int(best_d), direction


def annotate(df, gtf_path, max_dist):
    """批量注释邻近基因"""
    genes = parse_gtf(gtf_path)
    results = df.apply(
        lambda r: nearest_gene(str(r['chrom']), r['start'], r['end'], genes, max_dist),
        axis=1, result_type='expand'
    )
    results.columns = ['nearest_gene_id', 'gene_distance', 'gene_direction']
    return pd.concat([df, results], axis=1)


# =============================================================================
# 主流程
# =============================================================================

def main():
    args = parse_args()
    workdir = args.workdir
    print('=' * 60)
    print('  目标 sRNA 筛选 + 邻近基因注释')
    print(f'  padj<{args.padj}, |log2FC|>{args.lfc}, max_dist={args.max_dist}bp')
    print('=' * 60)

    print('\n[1] 加载 DE 结果 ...')
    de = load_de(workdir)

    print('\n[2] 筛选 + 合并 ...')
    df = filter_and_merge(de, workdir, args.padj, args.lfc)
    if len(df) == 0:
        print('[WARNING] 无符合条件的 sRNA，尝试放宽阈值')
        sys.exit(0)

    print('\n[3] 邻近基因注释 ...')
    gtf_path = args.gtf
    if Path(gtf_path).exists():
        df = annotate(df, gtf_path, args.max_dist)
        n_ok = (df['nearest_gene_id'] != 'NA').sum()
        print(f'  注释成功: {n_ok}/{len(df)}')
    else:
        print(f'  [SKIP] GTF 不存在: {gtf_path}')

    # 保存
    out = args.output or str(Path(workdir) / 'results' / 'target_sRNA_annotated.csv')
    df.to_csv(out, index=False)
    print(f'\n[完成] {out}')
    print(f'  共 {len(df)} 个目标 sRNA')

    # 摘要
    print('\n--- 摘要 ---')
    print(f'  sRNA类型: {dict(df["sRNA_type"].value_counts())}')
    if 'gene_direction' in df.columns:
        print(f'  基因方向: {dict(df["gene_direction"].value_counts())}')


if __name__ == '__main__':
    main()
