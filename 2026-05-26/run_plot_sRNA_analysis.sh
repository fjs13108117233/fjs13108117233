#!/bin/bash
#SBATCH -J plot_sRNA
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -o plot_sRNA_%j.out
#SBATCH -e plot_sRNA_%j.err

# ============================================================
# sRNA × DMR 联合分析可视化
# 提交: sbatch run_plot_sRNA_analysis.sh
# ============================================================

source ~/.bashrc
conda activate DL

SRNA="/public/home/h14166/fang/Heyufei/smrna_seq/sRNA_diff_analysis/results/target_sRNA_annotated.csv"
RPM="/public/home/h14166/fang/Heyufei/smrna_seq/sRNA_diff_analysis/counts/sRNA_rpm.csv"
META="/public/home/h14166/fang/Heyufei/smrna_seq/sRNA_diff_analysis/counts/sample_metadata.csv"
INTERSECT="/public/home/h14166/fang/Heyufei/Unite/sRNA_DMR_intersect.csv"
CONSISTENT="/public/home/h14166/fang/Heyufei/Unite/sRNA_DMR_intersect_consistent.csv"
OUTDIR="/public/home/h14166/fang/Heyufei/Unite/plots"
SCRIPT="/public/home/h14166/fang/Heyufei/Unite/tools/plot_sRNA_analysis.R"

echo "开始: $(date)"

Rscript ${SCRIPT} \
    --srna ${SRNA} \
    --rpm ${RPM} \
    --meta ${META} \
    --intersect ${INTERSECT} \
    --consistent ${CONSISTENT} \
    --outdir ${OUTDIR}

echo "完成: $(date)"
