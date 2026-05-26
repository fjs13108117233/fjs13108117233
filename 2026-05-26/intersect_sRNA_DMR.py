#!/usr/bin/env python3
"""
intersect_sRNA_DMR.py - sRNA 与 Nanopore DMR 趋势一致性交集分析

筛选逻辑:
  1. sRNA 邻近基因 ∩ DMR 注释基因 → 交集基因
  2. 趋势一致性过滤:
     sRNA: OE上调(0) + KO下调(1)
     DMR:  OE甲基化升高(0) + KO甲基化降低(1)
     → 同一基因的 sRNA 和 DMR 变化趋势必须一致

输入:
  - sRNA: target_sRNA_annotated.csv
  - DMR: filter_location/{4mc,5hmc,5mc,6ma}/*_genic.tsv

输出:
  - sRNA_DMR_intersect.csv           (全部交集，含趋势标记)
  - sRNA_DMR_intersect_consistent.csv (仅趋势一致的)

state 编码: 0=升高/上调, 1=降低/下调, 3=不显著

用法:
  python3 intersect_sRNA_DMR.py [--srna ...] [--dmr-dir ...] [--output ...]
"""

import argparse
import sys
import pandas as pd
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description='sRNA 与 DMR 趋势一致性交集分析')
    p.add_argument('--srna',
                   default='/public/home/h14166/fang/Heyufei/smrna_seq/sRNA_diff_analysis/results/target_sRNA_annotated.csv')
    p.add_argument('--dmr-dir',
                   default='/public/home/h14166/fang/Heyufei/Nanopore/filter_location')
    p.add_argument('--output',
                   default='/public/home/h14166/fang/Heyufei/smrna_seq/sRNA_diff_analysis/results/sRNA_DMR_intersect.csv')
    p.add_argument('--strict', action='store_true',
                   help='严格模式: DMR 四组状态必须全部与 sRNA 一致 (默认: 宽松模式，OE至少1个=0且KO至少1个=1)')
    return p.parse_args()


# =============================================================================
# 配置
# =============================================================================

# 修饰类型 -> 实际文件中间标识符
# 4mc: 4mc_annotated_DMR_afterselect_genic.tsv (无中间字母)
# 5hmc: 5hmc_H_annotated_DMR_afterselect_genic.tsv
# 5mc: 5mc_M_annotated_DMR_afterselect_genic.tsv
# 6ma: 6ma_A_annotated_DMR_afterselect_genic.tsv
DMR_MODS = {
    '4mc':  '4mc_annotated_DMR_afterselect_genic.tsv',
    '5hmc': '5hmc_H_annotated_DMR_afterselect_genic.tsv',
    '5mc':  '5mc_M_annotated_DMR_afterselect_genic.tsv',
    '6ma':  '6ma_A_annotated_DMR_afterselect_genic.tsv',
}

# DMR state 列名 -> 对应 sRNA 组
# sRNA 趋势: OE1=0(上调), OE2=0(上调), KO6bp=1(下调), KO8bp=1(下调)
# DMR 趋势应一致: OE1=0(hyper), OE2=0(hyper), KO6bp=1(hypo), KO8bp=1(hypo)
DMR_STATE_MAP = {
    'OE1_state': 0,         # 期望 OE1 甲基化升高
    'RDM_OE2_state': 0,     # 期望 OE2 甲基化升高
    'RDM_-6bp_state': 1,    # 期望 KO6bp 甲基化降低
    'RDM_-8bp_state': 1,    # 期望 KO8bp 甲基化降低
}


# =============================================================================
# 数据加载
# =============================================================================

def load_srna(path):
    """加载 sRNA 结果，仅保留有基因注释的"""
    df = pd.read_csv(path)
    df = df[df['nearest_gene_id'] != 'NA'].copy()
    print(f'  sRNA: {len(df)} 条 (有基因注释), 涉及 {df["nearest_gene_id"].nunique()} 个基因')
    return df


def load_dmr(dmr_dir, mod, filename):
    """加载某种修饰的 genic DMR"""
    genic_file = Path(dmr_dir) / mod / filename
    if not genic_file.exists():
        print(f'  [SKIP] 文件不存在: {genic_file}')
        return None

    # 读表头（带 # 号）
    with open(genic_file) as f:
        header = f.readline().lstrip('#').strip().split('\t')

    df = pd.read_csv(genic_file, sep='\t', comment='#', header=None, low_memory=False)
    df.columns = header[:len(df.columns)]

    print(f'  {mod}: {len(df)} DMRs, {df["Gene"].nunique()} genes')
    return df


# =============================================================================
# 趋势一致性判断
# =============================================================================

def check_trend_consistent(dmr_gene_df, strict=False):
    """
    判断某基因的 DMR 是否与 sRNA 趋势一致

    strict=True:  OE1=0 & OE2=0 & KO6bp=1 & KO8bp=1 (全部一致)
    strict=False: OE中至少1个=0, KO中至少1个=1 (宽松)

    返回: True/False
    """
    # 每个基因可能有多个 DMR 区域，取各 state 的众数(mode)
    available_cols = [c for c in DMR_STATE_MAP.keys() if c in dmr_gene_df.columns]
    if not available_cols:
        return False

    # 取该基因所有 DMR 区域中各 state 的众数
    states = {}
    for col in available_cols:
        vals = dmr_gene_df[col].dropna()
        if len(vals) > 0:
            mode = vals.mode()
            states[col] = int(mode.iloc[0]) if len(mode) > 0 else 3
        else:
            states[col] = 3

    if strict:
        # 严格: 每个 state 都必须等于期望值
        for col, expected in DMR_STATE_MAP.items():
            if col in states and states[col] != expected:
                return False
            elif col not in states:
                return False
        return True
    else:
        # 宽松: OE 中至少一个=0, KO 中至少一个=1
        oe_cols = [c for c in available_cols if 'OE' in c]
        ko_cols = [c for c in available_cols if 'bp' in c]

        oe_ok = any(states.get(c, 3) == 0 for c in oe_cols)
        ko_ok = any(states.get(c, 3) == 1 for c in ko_cols)
        return oe_ok and ko_ok


# =============================================================================
# 主流程
# =============================================================================

def main():
    args = parse_args()
    mode_str = '严格' if args.strict else '宽松'

    print('=' * 60)
    print('  sRNA × DMR 趋势一致性交集分析')
    print('=' * 60)
    print(f'  sRNA:   {args.srna}')
    print(f'  DMR:    {args.dmr_dir}')
    print(f'  模式:   {mode_str}')
    print(f'  输出:   {args.output}')
    print('=' * 60)

    # 1. 加载 sRNA
    print('\n[1] 加载 sRNA ...')
    srna_df = load_srna(args.srna)
    srna_genes = set(srna_df['nearest_gene_id'].unique())

    if not srna_genes:
        sys.exit('[ERROR] 无有效 sRNA 基因')

    # 2. 逐修饰分析
    print('\n[2] 逐修饰类型分析 ...')
    results_by_mod = {}  # mod -> DataFrame of consistent genes

    for mod, filename in DMR_MODS.items():
        print(f'\n  --- {mod.upper()} ---')
        dmr_df = load_dmr(args.dmr_dir, mod, filename)
        if dmr_df is None:
            continue

        # 基因交集
        dmr_genes = set(dmr_df['Gene'].unique())
        common = srna_genes & dmr_genes
        print(f'    基因交集: {len(common)}')

        if not common:
            continue

        # 趋势一致性过滤
        consistent_genes = []
        gene_dmr_info = {}

        for gene in common:
            gene_dmr = dmr_df[dmr_df['Gene'] == gene]
            if check_trend_consistent(gene_dmr, strict=args.strict):
                consistent_genes.append(gene)

                # 记录该基因 DMR 摘要
                state_cols = [c for c in DMR_STATE_MAP.keys() if c in gene_dmr.columns]
                info = {'Gene': gene, f'{mod}_n_DMR': len(gene_dmr)}
                for col in state_cols:
                    vals = gene_dmr[col].dropna()
                    info[f'{mod}_{col}'] = int(vals.mode().iloc[0]) if len(vals) > 0 else 3
                gene_dmr_info[gene] = info

        print(f'    趋势一致: {len(consistent_genes)} 个基因')
        results_by_mod[mod] = {
            'genes': set(consistent_genes),
            'info': gene_dmr_info,
        }

    # 3. 合并结果
    print('\n[3] 合并结果 ...')
    merged = srna_df.copy()

    for mod in DMR_MODS:
        if mod not in results_by_mod:
            merged[f'{mod}_consistent'] = 0
            continue

        info_dict = results_by_mod[mod]['info']
        info_df = pd.DataFrame(info_dict.values())

        if len(info_df) > 0:
            merged = merged.merge(info_df, left_on='nearest_gene_id', right_on='Gene', how='left')
            if 'Gene' in merged.columns:
                merged = merged.drop(columns=['Gene'])

        # 标记是否趋势一致
        col = f'{mod}_n_DMR'
        merged[f'{mod}_consistent'] = merged[col].notna().astype(int) if col in merged.columns else 0

    # 总计一致的修饰数
    consist_cols = [f'{mod}_consistent' for mod in DMR_MODS]
    merged['n_consistent_mods'] = merged[consist_cols].sum(axis=1)

    # 4. 输出
    # 全部
    merged.to_csv(args.output, index=False)

    # 仅趋势一致的子集
    consistent_df = merged[merged['n_consistent_mods'] > 0].copy()
    out_consistent = args.output.replace('.csv', '_consistent.csv')
    consistent_df.to_csv(out_consistent, index=False)

    # 5. 汇总
    print('\n' + '=' * 60)
    print('  结果汇总')
    print('=' * 60)
    print(f'  sRNA 目标基因总数: {len(srna_genes)}')
    print(f'  趋势一致的 sRNA 条数: {len(consistent_df)}')
    print(f'  趋势一致的基因数: {consistent_df["nearest_gene_id"].nunique() if len(consistent_df) > 0 else 0}')
    print()
    print(f'  各修饰趋势一致基因数 ({mode_str}模式):')
    for mod in DMR_MODS:
        n = len(results_by_mod.get(mod, {}).get('genes', set()))
        print(f'    {mod.upper():>5}: {n}')

    # 多修饰共同一致
    all_consistent = [results_by_mod[m]['genes'] for m in DMR_MODS if m in results_by_mod and results_by_mod[m]['genes']]
    if len(all_consistent) > 1:
        multi = set.intersection(*all_consistent)
        print(f'\n  所有修饰共同一致: {len(multi)} 个基因')
        if multi:
            for g in sorted(multi)[:20]:
                print(f'    {g}')

    print(f'\n[完成]')
    print(f'  全部: {args.output}')
    print(f'  一致: {out_consistent}')


if __name__ == '__main__':
    main()
