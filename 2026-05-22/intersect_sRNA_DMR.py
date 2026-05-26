#!/usr/bin/env python3
"""
intersect_sRNA_DMR.py - sRNA 目标基因与 Nanopore DMR 基因联合分析

功能: 找到 sRNA 邻近基因与 DMR 注释基因的交集
输入:
  - sRNA: target_sRNA_annotated.csv (filter_target_sRNA.py 输出)
  - DMR:  4种修饰类型 (4mc/5hmc/5mc/6ma) 的 genic DMR 文件

输出:
  - 交集基因汇总表 (含 sRNA 信息 + DMR 信息)
  - 各修饰类型交集统计

用法:
  python3 intersect_sRNA_DMR.py [--srna ...] [--dmr-dir ...] [--output ...]
"""

import argparse
import sys
import pandas as pd
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description='sRNA 与 DMR 交集基因分析')
    p.add_argument('--srna',
                   default='/public/home/h14166/fang/Heyufei/smrna_seq/sRNA_diff_analysis/results/target_sRNA_annotated.csv',
                   help='sRNA 筛选结果文件')
    p.add_argument('--dmr-dir',
                   default='/public/home/h14166/fang/Heyufei/Nanopore/filter_location',
                   help='DMR 数据目录 (含 4mc/5hmc/5mc/6ma 子目录)')
    p.add_argument('--output',
                   default='/public/home/h14166/fang/Heyufei/smrna_seq/sRNA_diff_analysis/results/sRNA_DMR_intersect.csv',
                   help='输出文件路径')
    return p.parse_args()


# =============================================================================
# 数据加载
# =============================================================================

# DMR 文件命名规则: {mod}_H_annotated_DMR_afterselect_genic.tsv
DMR_MODS = ['4mc', '5hmc', '5mc', '6ma']


def load_srna(srna_path):
    """加载 sRNA 筛选结果"""
    df = pd.read_csv(srna_path)
    # 只保留有基因注释的行
    df = df[df['nearest_gene_id'] != 'NA'].copy()
    print(f'  sRNA: {len(df)} 条有基因注释')
    print(f'  涉及基因: {df["nearest_gene_id"].nunique()} 个')
    return df


def load_dmr(dmr_dir, mod):
    """加载某种修饰的 genic DMR 文件"""
    # 尝试两种文件名格式
    genic_file = Path(dmr_dir) / mod / f'{mod}_H_annotated_DMR_afterselect_genic.tsv'
    if not genic_file.exists():
        print(f'  [WARNING] 文件不存在: {genic_file}')
        return None

    df = pd.read_csv(genic_file, sep='\t', comment='#',
                     header=None, low_memory=False)

    # 读取表头（第一行带 #）
    with open(genic_file) as f:
        header_line = f.readline().lstrip('#').strip()
    cols = header_line.split('\t')
    df.columns = cols[:len(df.columns)]

    print(f'  {mod} genic DMR: {len(df)} 条, 涉及 {df["Gene"].nunique()} 个基因')
    return df


# =============================================================================
# 交集分析
# =============================================================================

def find_intersect(srna_df, dmr_df, mod_name):
    """找到 sRNA 基因与 DMR 基因的交集"""
    srna_genes = set(srna_df['nearest_gene_id'].unique())
    dmr_genes = set(dmr_df['Gene'].unique())

    common_genes = srna_genes & dmr_genes
    print(f'    交集基因: {len(common_genes)}')

    if len(common_genes) == 0:
        return pd.DataFrame()

    # 提取交集基因的 sRNA 信息
    srna_hit = srna_df[srna_df['nearest_gene_id'].isin(common_genes)].copy()

    # 提取交集基因的 DMR 信息（每个基因取代表性 DMR：选 Delta 绝对值最大的）
    dmr_hit = dmr_df[dmr_df['Gene'].isin(common_genes)].copy()

    # DMR 汇总: 每个基因的 DMR 数目 + 各组状态
    dmr_summary = dmr_hit.groupby('Gene').agg(
        n_DMR=(dmr_hit.columns[0], 'count'),  # DMR 区域数
    ).reset_index()

    # 获取各组状态列 (如果有)
    state_cols = [c for c in dmr_df.columns if c.endswith('_state')]
    if state_cols:
        # 每个基因取最常见状态
        for col in state_cols:
            mode_vals = dmr_hit.groupby('Gene')[col].agg(
                lambda x: x.mode().iloc[0] if len(x.mode()) > 0 else 3
            ).reset_index()
            mode_vals.columns = ['Gene', f'{mod_name}_{col}']
            dmr_summary = dmr_summary.merge(mode_vals, on='Gene', how='left')

    # 合并 DMR 数目到 sRNA 结果
    dmr_summary = dmr_summary.rename(columns={'n_DMR': f'{mod_name}_n_DMR'})
    result = srna_hit.merge(dmr_summary, left_on='nearest_gene_id', right_on='Gene', how='left')
    if 'Gene' in result.columns:
        result = result.drop(columns=['Gene'])

    return result


# =============================================================================
# 主流程
# =============================================================================

def main():
    args = parse_args()

    print('=' * 60)
    print('  sRNA 与 Nanopore DMR 交集基因分析')
    print('=' * 60)
    print(f'  sRNA 文件: {args.srna}')
    print(f'  DMR 目录:  {args.dmr_dir}')
    print(f'  输出文件:  {args.output}')
    print('=' * 60)

    # 1. 加载 sRNA
    print('\n[1] 加载 sRNA 结果 ...')
    srna_df = load_srna(args.srna)
    srna_genes = set(srna_df['nearest_gene_id'].unique())

    if len(srna_df) == 0:
        print('[ERROR] sRNA 无有效基因注释')
        sys.exit(1)

    # 2. 逐个修饰类型分析
    print('\n[2] 加载 DMR 数据并找交集 ...')
    all_intersect_genes = {}  # mod -> set of genes
    merged_results = srna_df.copy()

    for mod in DMR_MODS:
        print(f'\n  --- {mod.upper()} ---')
        dmr_df = load_dmr(args.dmr_dir, mod)
        if dmr_df is None or len(dmr_df) == 0:
            continue

        dmr_genes = set(dmr_df['Gene'].unique())
        common = srna_genes & dmr_genes
        all_intersect_genes[mod] = common
        print(f'    sRNA基因: {len(srna_genes)}, DMR基因: {len(dmr_genes)}, 交集: {len(common)}')

        if len(common) > 0:
            # 合并 DMR 数目信息
            dmr_count = dmr_df[dmr_df['Gene'].isin(common)].groupby('Gene').size().reset_index(name=f'{mod}_n_DMR')

            # 获取状态列
            state_cols = [c for c in dmr_df.columns if c.endswith('_state')]
            if state_cols:
                state_summary = dmr_df[dmr_df['Gene'].isin(common)].groupby('Gene')[state_cols].agg(
                    lambda x: x.mode().iloc[0] if len(x.mode()) > 0 else 3
                ).reset_index()
                state_summary.columns = ['Gene'] + [f'{mod}_{c}' for c in state_cols]
                dmr_count = dmr_count.merge(state_summary, on='Gene', how='left')

            merged_results = merged_results.merge(
                dmr_count, left_on='nearest_gene_id', right_on='Gene', how='left'
            )
            if 'Gene' in merged_results.columns:
                merged_results = merged_results.drop(columns=['Gene'])

    # 3. 标记哪些基因在各修饰中有 DMR
    for mod in DMR_MODS:
        col = f'{mod}_n_DMR'
        if col in merged_results.columns:
            merged_results[f'{mod}_has_DMR'] = merged_results[col].notna().astype(int)
        else:
            merged_results[f'{mod}_has_DMR'] = 0

    # 4. 总交集: 至少在一种修饰中有 DMR 的基因
    dmr_flag_cols = [f'{mod}_has_DMR' for mod in DMR_MODS]
    merged_results['any_DMR'] = merged_results[dmr_flag_cols].sum(axis=1)
    has_any_dmr = merged_results[merged_results['any_DMR'] > 0]

    # 5. 保存
    # 完整表
    merged_results.to_csv(args.output, index=False)

    # 只保留有 DMR 交集的子集
    out_intersect = args.output.replace('.csv', '_filtered.csv')
    has_any_dmr.to_csv(out_intersect, index=False)

    # 6. 汇总统计
    print('\n' + '=' * 60)
    print('  结果汇总')
    print('=' * 60)
    print(f'  sRNA 目标基因总数: {len(srna_genes)}')
    print(f'  与任意 DMR 有交集的 sRNA: {len(has_any_dmr)} 条')
    print(f'  与任意 DMR 有交集的基因: {has_any_dmr["nearest_gene_id"].nunique()} 个')
    print()
    print('  各修饰交集基因数:')
    for mod in DMR_MODS:
        n = len(all_intersect_genes.get(mod, set()))
        print(f'    {mod.upper():>5}: {n} 个基因')

    # 多重修饰交集
    if len(all_intersect_genes) > 1:
        all_mod_genes = list(all_intersect_genes.values())
        multi_intersect = set.intersection(*all_mod_genes) if all_mod_genes else set()
        print(f'\n  所有修饰共同交集: {len(multi_intersect)} 个基因')
        if multi_intersect:
            print(f'    基因列表: {sorted(multi_intersect)[:20]}')

    print(f'\n[完成]')
    print(f'  全部结果: {args.output}')
    print(f'  交集子集: {out_intersect}')


if __name__ == '__main__':
    main()
