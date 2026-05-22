#!/bin/bash
#SBATCH -p com300
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -J sRNA_diff
#SBATCH --output logs/sRNA_diff_%j.out
#SBATCH --error  logs/sRNA_diff_%j.err

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"
mkdir -p logs

source "$HOME/anaconda3/etc/profile.d/conda.sh"
conda activate RNA

# =========================
# 用户参数
# =========================
WORKDIR="/public/home/h14166/fang/Heyufei/smrna_seq"
TOOLS="${WORKDIR}/tools"
OUTDIR="${WORKDIR}/sRNA_diff_analysis"

# =========================
# Step 1: 构建 count 矩阵
# =========================
echo ">>> Step 1: 构建 count 矩阵 ..."

python3 ${TOOLS}/sRNA_diff_pipeline.py \
    --workdir ${WORKDIR} \
    --bed ${WORKDIR}/all_sRNA.bed \
    --list ${WORKDIR}/list \
    --outdir ${OUTDIR} \
    --min-samples 2 \
    --rscript ${TOOLS}/sRNA_deseq2_analysis.R \
    --skip-r

# =========================
# Step 2: DESeq2 差异分析
# =========================
echo ""
echo ">>> Step 2: DESeq2 差异分析 ..."

Rscript ${TOOLS}/sRNA_deseq2_analysis.R ${OUTDIR}

echo ""
echo "=== 全部完成: $(date) ==="
echo "输出目录: ${OUTDIR}"
