#!/usr/bin/env python3
"""
filter_target_sRNA.py
=====================
筛选目标 sRNA 并进行邻近基因注释

筛选逻辑:
  - 在 OE1 和 OE2 两个样本中 **均显著上调**（不能有下调）
  - 在 KO(6bp) 和 KO(8bp) 两个样本中 **均显著下调**（不能有上调）

输出格式:
  chrom, start, end, sRNA_type,
  OE1_RPM, OE1_pvalue, OE1_log2FC, OE1_state,
  OE2_RPM, OE2_pvalue, OE2_log2FC, OE2_state,
  KO6bp_RPM, KO6bp_pvalue, KO6bp_log2FC, KO6bp_state,
  KO8bp_RPM, KO8bp_pvalue, KO8bp_log2FC, KO8bp_state,
  nearest_gene_id, nearest_gene_name, gene_distance, gene_direction

state 编码: 上调=0, 下调=1, 不显著=3

用法:
  python3 filter_target_sRNA.py [--workdir /path/to/sRNA_diff_analysis] [--gtf /path/to/maizev4.gtf]
"""

import argparse
import sys
import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict


# =============================================================================
# 参数
# =============================================================================
def parse_args():
    p = argparse.ArgumentParser(description='筛选目标 sRNA + 邻近基因注释')
    p.add_argument('--workdir',
                   default='/public/home/h14166/fang/Heyufei/smrna_seq/sRNA_diff_analysis',
                   help='差异分析结果目录')
    p.add_argument('--gtf',
                   default='/public/home/h14166/refer/maizev4.gtf',
                   help='GTF 注释文件路径')
    p.add_argument('--padj', type=float, default=0.05,
                   help='padj 阈值 (默认: 0.05)')
    p.add_argument('--lfc', type=float, default=1.0,
                   help='|log2FoldChange| 阈值 (默认: 1.0)')
    p.add_argument('--max-dist', type=int, default=10000,
                   help='邻近基因最大距离 bp (默认: 10000)')
    p.add_argument('--output', default=None,
                   help='输出文件路径 (默认: $WORKDIR/results/target_sRNA_annotated.csv)')
    return p.parse_args()


# =============================================================================
# Step 1: 加载 DE 结果
# =============================================================================
def load_de_results(workdir):
    """加载各组的差异表达结果"""
    results_dir = Path(workdir) / 'results'

    # 定义需要的文件
    files = {
        'OE1': results_dir / 'DE_RDM_OE1_vs_WT.csv',
        'OE2': results_dir / 'DE_RDM_OE2_vs_WT.csv',
        'KO6bp': results_dir / 'DE_RDM_6bp_KO_vs_WT.csv',
        'KO8bp': results_dir / 'DE_RDM_8bp_KO_vs_WT.csv',
    }

    de_data = {}
    for name, fpath in files.items():
        if not fpath.exists():
            sys.exit(f'[ERROR] 文件不存在: {fpath}')
        df = pd.read_csv(fpath)
        df = df.set_index('sRNA_name')
        de_data[name] = df
        print(f'  加载 {name}: {len(df)} sRNAs')

    return de_data


# =============================================================================
# Step 2: 筛选目标 sRNA
# =============================================================================
def classify_state(row, padj_cutoff, lfc_cutoff):
    """
    判断 sRNA 状态:
      上调 = 0 (padj < cutoff & log2FC > lfc_cutoff)
      下调 = 1 (padj < cutoff & log2FC < -lfc_cutoff)
      不显著 = 3
    """
    if pd.isna(row['padj']):
        return 3
    if row['padj'] < padj_cutoff and row['log2FoldChange'] > lfc_cutoff:
        return 0  # 上调
    elif row['padj'] < padj_cutoff and row['log2FoldChange'] < -lfc_cutoff:
        return 1  # 下调
    else:
        return 3  # 不显著


def filter_target_srna(de_data, padj_cutoff, lfc_cutoff):
    """
    筛选条件:
      - OE1 和 OE2: 均为显著上调 (state=0)，不能有下调
      - KO6bp 和 KO8bp: 均为显著下调 (state=1)，不能有上调
    """
    # 获取所有 sRNA 名称的交集
    all_srnas = set(de_data['OE1'].index)
    for name in ['OE2', 'KO6bp', 'KO8bp']:
        all_srnas &= set(de_data[name].index)

    print(f'\n  四组共有 sRNA 数: {len(all_srnas)}')

    # 为每个 sRNA 判断各组状态
    target_srnas = []
    for srna in all_srnas:
        states = {}
        for name in ['OE1', 'OE2', 'KO6bp', 'KO8bp']:
            row = de_data[name].loc[srna]
            states[name] = classify_state(row, padj_cutoff, lfc_cutoff)

        # 筛选条件: OE1+OE2 均上调(0), KO6bp+KO8bp 均下调(1)
        oe_ok = (states['OE1'] == 0) and (states['OE2'] == 0)
        ko_ok = (states['KO6bp'] == 1) and (states['KO8bp'] == 1)

        if oe_ok and ko_ok:
            target_srnas.append(srna)

    print(f'  筛选后目标 sRNA 数: {len(target_srnas)}')
    return target_srnas


# =============================================================================
# Step 3: 合并结果为指定格式
# =============================================================================
def build_merged_table(target_srnas, de_data, workdir, padj_cutoff, lfc_cutoff):
    """
    构建合并表格:
    chrom, start, end, sRNA_type,
    {sample}_RPM, {sample}_pvalue, {sample}_log2FC, {sample}_state
    """
    # 加载 RPM 矩阵
    rpm_file = Path(workdir) / 'counts' / 'sRNA_rpm.csv'
    if rpm_file.exists():
        rpm_df = pd.read_csv(rpm_file, index_col=0)
        print(f'  RPM 矩阵: {rpm_df.shape[0]} x {rpm_df.shape[1]}')
    else:
        print(f'  [WARNING] RPM 文件不存在: {rpm_file}, 使用 baseMean 替代')
        rpm_df = None

    # 加载注释信息
    anno_file = Path(workdir) / 'counts' / 'sRNA_annotation.csv'
    if anno_file.exists():
        anno_df = pd.read_csv(anno_file)
        anno_df = anno_df.set_index('sRNA_name')
    else:
        anno_df = None

    # 样本名映射 (DE 比较组名 -> RPM 中的实际样本列名)
    # 需要从 RPM 列中找到各组的样本
    sample_groups = {
        'OE1': 'RDM_OE1',
        'OE2': 'RDM_OE2',
        'KO6bp': 'RDM_6bp_KO',
        'KO8bp': 'RDM_8bp_KO',
    }

    # 尝试从 RPM 列名中匹配实际样本
    rpm_sample_map = {}
    if rpm_df is not None:
        for group_key, pattern in sample_groups.items():
            matched = [c for c in rpm_df.columns if pattern in c]
            if matched:
                rpm_sample_map[group_key] = matched
            else:
                # 尝试更宽松匹配
                rpm_sample_map[group_key] = []

    rows = []
    for srna in target_srnas:
        row = {}

        # 基本坐标信息
        # 从 DE 数据或 annotation 获取坐标
        ref = de_data['OE1'].loc[srna]
        if anno_df is not None and srna in anno_df.index:
            anno_row = anno_df.loc[srna]
            if isinstance(anno_row, pd.DataFrame):
                anno_row = anno_row.iloc[0]
            row['chrom'] = anno_row.get('chrom', ref.get('chrom', ''))
            row['start'] = int(anno_row.get('start', ref.get('start', 0)))
            row['end'] = int(anno_row.get('end', ref.get('end', 0)))
            row['sRNA_type'] = anno_row.get('sRNA_type', ref.get('sRNA_type', ''))
        else:
            row['chrom'] = ref.get('chrom', '')
            row['start'] = int(ref.get('start', 0))
            row['end'] = int(ref.get('end', 0))
            row['sRNA_type'] = ref.get('sRNA_type', '')

        # 各组的统计量
        for group_key in ['OE1', 'OE2', 'KO6bp', 'KO8bp']:
            de_row = de_data[group_key].loc[srna]
            state = classify_state(de_row, padj_cutoff, lfc_cutoff)

            # RPM: 取该组所有重复样本的平均 RPM
            if rpm_df is not None and group_key in rpm_sample_map and rpm_sample_map[group_key]:
                sample_cols = rpm_sample_map[group_key]
                if srna in rpm_df.index:
                    rpm_val = rpm_df.loc[srna, sample_cols].mean()
                else:
                    rpm_val = 0.0
            else:
                # 使用 baseMean 作为替代
                rpm_val = de_row.get('baseMean', 0.0)

            row[f'{group_key}_RPM'] = round(rpm_val, 4)
            row[f'{group_key}_pvalue'] = de_row.get('pvalue', np.nan)
            row[f'{group_key}_log2FC'] = round(de_row.get('log2FoldChange', 0), 4)
            row[f'{group_key}_state'] = state

        rows.append(row)

    result_df = pd.DataFrame(rows)

    # 排序
    if len(result_df) > 0:
        result_df = result_df.sort_values(['chrom', 'start']).reset_index(drop=True)

    return result_df


# =============================================================================
# Step 4: GTF 解析 + 邻近基因注释
# =============================================================================
def parse_gtf_genes(gtf_path):
    """
    解析 GTF 文件，提取 gene 记录
    返回: dict[chrom] -> list of (start, end, strand, gene_id, gene_name)
    """
    print(f'  解析 GTF: {gtf_path}')
    genes = defaultdict(list)
    n_gene = 0

    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            fields = line.rstrip().split('\t')
            if len(fields) < 9:
                continue
            if fields[2] != 'gene':
                continue

            chrom = fields[0]
            start = int(fields[3])  # 1-based
            end = int(fields[4])
            strand = fields[6]

            # 解析属性
            attrs = fields[8]
            gene_id = ''
            gene_name = ''

            # 提取 gene_id
            for attr in attrs.split(';'):
                attr = attr.strip()
                if attr.startswith('gene_id'):
                    gene_id = attr.split('"')[1] if '"' in attr else attr.split(' ')[-1]
                elif attr.startswith('gene_name'):
                    gene_name = attr.split('"')[1] if '"' in attr else attr.split(' ')[-1]

            if not gene_name:
                gene_name = gene_id

            genes[chrom].append((start, end, strand, gene_id, gene_name))
            n_gene += 1

    # 按起始位置排序
    for chrom in genes:
        genes[chrom].sort(key=lambda x: x[0])

    print(f'  共解析 {n_gene} 个基因, 分布于 {len(genes)} 条染色体')
    return genes


def find_nearest_gene(chrom, srna_start, srna_end, genes_dict, max_dist):
    """
    查找距离 sRNA 最近的基因

    返回: (gene_id, gene_name, distance, direction)
      distance: sRNA 与基因边界的最短距离 (重叠为 0)
      direction: upstream / downstream / overlap
    """
    if chrom not in genes_dict:
        return ('NA', 'NA', 'NA', 'NA')

    gene_list = genes_dict[chrom]
    srna_mid = (srna_start + srna_end) / 2.0

    best_dist = float('inf')
    best_gene = None

    # 二分搜索找起始附近位置
    import bisect
    starts = [g[0] for g in gene_list]
    idx = bisect.bisect_left(starts, srna_start)

    # 搜索附近范围
    search_range = range(max(0, idx - 50), min(len(gene_list), idx + 50))

    for i in search_range:
        g_start, g_end, g_strand, g_id, g_name = gene_list[i]

        # 计算距离
        if srna_end <= g_start:
            # sRNA 在基因上游 (按坐标)
            dist = g_start - srna_end
        elif srna_start >= g_end:
            # sRNA 在基因下游 (按坐标)
            dist = srna_start - g_end
        else:
            # 重叠
            dist = 0

        if dist < best_dist:
            best_dist = dist
            best_gene = gene_list[i]

    if best_gene is None or best_dist > max_dist:
        return ('NA', 'NA', 'NA', 'NA')

    g_start, g_end, g_strand, g_id, g_name = best_gene

    # 判断方向 (相对于基因)
    if best_dist == 0:
        direction = 'overlap'
    elif srna_end <= g_start:
        # sRNA 在基因起始之前
        direction = 'upstream' if g_strand == '+' else 'downstream'
    else:
        # sRNA 在基因结束之后
        direction = 'downstream' if g_strand == '+' else 'upstream'

    return (g_id, g_name, int(best_dist), direction)


def annotate_nearest_genes(result_df, gtf_path, max_dist):
    """为每个 sRNA 注释最近基因"""
    genes_dict = parse_gtf_genes(gtf_path)

    gene_ids = []
    gene_names = []
    distances = []
    directions = []

    for _, row in result_df.iterrows():
        g_id, g_name, dist, dirn = find_nearest_gene(
            row['chrom'], row['start'], row['end'], genes_dict, max_dist
        )
        gene_ids.append(g_id)
        gene_names.append(g_name)
        distances.append(dist)
        directions.append(dirn)

    result_df['nearest_gene_id'] = gene_ids
    result_df['nearest_gene_name'] = gene_names
    result_df['gene_distance'] = distances
    result_df['gene_direction'] = directions

    return result_df


# =============================================================================
# 主流程
# =============================================================================
def main():
    args = parse_args()

    workdir = args.workdir
    gtf_path = args.gtf
    padj_cutoff = args.padj
    lfc_cutoff = args.lfc
    max_dist = args.max_dist

    print('=' * 70)
    print('  目标 sRNA 筛选 + 邻近基因注释')
    print('=' * 70)
    print(f'  工作目录:    {workdir}')
    print(f'  GTF 文件:    {gtf_path}')
    print(f'  padj 阈值:   {padj_cutoff}')
    print(f'  |log2FC| 阈值: {lfc_cutoff}')
    print(f'  最大注释距离: {max_dist} bp')
    print('=' * 70)

    # ---- Step 1: 加载 DE 结果 ----
    print('\n[Step 1] 加载差异表达结果 ...')
    de_data = load_de_results(workdir)

    # ---- Step 2: 筛选目标 sRNA ----
    print('\n[Step 2] 筛选目标 sRNA ...')
    print(f'  条件: OE1+OE2 均显著上调, KO6bp+KO8bp 均显著下调')
    target_srnas = filter_target_srna(de_data, padj_cutoff, lfc_cutoff)

    if len(target_srnas) == 0:
        print('\n[WARNING] 未找到符合条件的目标 sRNA!')
        print('  可尝试放宽阈值: --padj 0.1 --lfc 0.5')
        sys.exit(0)

    # ---- Step 3: 合并结果表 ----
    print('\n[Step 3] 构建合并结果表 ...')
    result_df = build_merged_table(target_srnas, de_data, workdir, padj_cutoff, lfc_cutoff)
    print(f'  结果行数: {len(result_df)}')

    # ---- Step 4: 邻近基因注释 ----
    print('\n[Step 4] 邻近基因注释 ...')
    if not Path(gtf_path).exists():
        print(f'  [WARNING] GTF 文件不存在: {gtf_path}, 跳过注释')
    else:
        result_df = annotate_nearest_genes(result_df, gtf_path, max_dist)
        n_annotated = (result_df['nearest_gene_id'] != 'NA').sum()
        print(f'  成功注释: {n_annotated}/{len(result_df)} sRNAs')

    # ---- Step 5: 保存结果 ----
    out_path = args.output or str(Path(workdir) / 'results' / 'target_sRNA_annotated.csv')
    result_df.to_csv(out_path, index=False)
    print(f'\n[完成] 结果保存至: {out_path}')
    print(f'  共 {len(result_df)} 个目标 sRNA')

    # 打印简要统计
    if len(result_df) > 0:
        print('\n--- 结果摘要 ---')
        print(f'  染色体分布:')
        for chrom, cnt in result_df['chrom'].value_counts().head(10).items():
            print(f'    {chrom}: {cnt}')
        if 'sRNA_type' in result_df.columns:
            print(f'  sRNA 类型:')
            for t, cnt in result_df['sRNA_type'].value_counts().items():
                print(f'    {t}: {cnt}')
        if 'gene_direction' in result_df.columns:
            print(f'  基因方向分布:')
            for d, cnt in result_df['gene_direction'].value_counts().items():
                print(f'    {d}: {cnt}')


if __name__ == '__main__':
    main()
