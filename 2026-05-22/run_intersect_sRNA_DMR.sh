#!/bin/bash
#SBATCH -J sRNA_DMR
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -o sRNA_DMR_%j.out
#SBATCH -e sRNA_DMR_%j.err

# ============================================================
# sRNA 与 Nanopore DMR 交集分析
# 提交: sbatch run_intersect_sRNA_DMR.sh
# ============================================================

source ~/.bashrc
conda activate DL

SRNA="/public/home/h14166/fang/Heyufei/smrna_seq/sRNA_diff_analysis/results/target_sRNA_annotated.csv"
DMR_DIR="/public/home/h14166/fang/Heyufei/Nanopore/filter_location"
OUTPUT="/public/home/h14166/fang/Heyufei/smrna_seq/sRNA_diff_analysis/results/sRNA_DMR_intersect.csv"
SCRIPT="/public/home/h14166/fang/Heyufei/smrna_seq/tools/intersect_sRNA_DMR.py"

echo "开始: $(date)"

python3 ${SCRIPT} \
    --srna ${SRNA} \
    --dmr-dir ${DMR_DIR} \
    --output ${OUTPUT}

echo "完成: $(date)"
