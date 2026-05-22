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
OUTDIR="${WORKDIR}/sRNA_diff_analysis"
BED="${WORKDIR}/all_sRNA.bed"
LIST="${WORKDIR}/list"
MIN_SAMPLES=2

# 脚本路径 (与本脚本同目录)
SCRIPT_DIR="${SLURM_SUBMIT_DIR}"
BUILD_PY="${SCRIPT_DIR}/build_sRNA_counts.py"
DESEQ2_R="${SCRIPT_DIR}/sRNA_deseq2_analysis.R"

# =========================
# 检查
# =========================
echo "=== sRNA 差异分析流水线 ==="
echo "开始时间: $(date)"
echo "工作目录: $WORKDIR"
echo "输出目录: $OUTDIR"
echo "BED 文件: $BED"
echo "最小样本数: $MIN_SAMPLES"
echo ""

[ -f "$BED" ]       || { echo "[ERROR] BED 文件不存在: $BED" >&2; exit 1; }
[ -f "$LIST" ]      || { echo "[ERROR] 样本列表不存在: $LIST" >&2; exit 1; }
[ -f "$BUILD_PY" ]  || { echo "[ERROR] Python 脚本不存在: $BUILD_PY" >&2; exit 1; }
[ -f "$DESEQ2_R" ]  || { echo "[ERROR] R 脚本不存在: $DESEQ2_R" >&2; exit 1; }

command -v python3 >/dev/null 2>&1 || { echo "[ERROR] python3 not found" >&2; exit 1; }
command -v Rscript >/dev/null 2>&1 || { echo "[ERROR] Rscript not found" >&2; exit 1; }

# =========================
# Step 1: 构建 count 矩阵
# =========================
echo "============================================"
echo ">>> Step 1: 构建 sRNA count 矩阵"
echo "============================================"
echo ""

python3 "$BUILD_PY" \
    --bed "$BED" \
    --workdir "$WORKDIR" \
    --list "$LIST" \
    --outdir "$OUTDIR" \
    --min-samples "$MIN_SAMPLES"

echo ""

# 检查输出
if [ ! -f "${OUTDIR}/counts/sRNA_counts.csv" ]; then
    echo "[ERROR] count 矩阵未生成" >&2
    exit 1
fi

# =========================
# Step 2: DESeq2 差异分析
# =========================
echo "============================================"
echo ">>> Step 2: DESeq2 差异分析"
echo "============================================"
echo ""

Rscript "$DESEQ2_R" "$OUTDIR"

echo ""

# =========================
# 完成
# =========================
echo "============================================"
echo "=== 全部完成: $(date) ==="
echo "============================================"
echo ""
echo "输出文件:"
echo "  Counts:  ${OUTDIR}/counts/"
echo "  Results: ${OUTDIR}/results/"
echo "  Plots:   ${OUTDIR}/plots/"
echo ""
echo "结果文件列表:"
ls -lh "${OUTDIR}/results/"*.csv 2>/dev/null || echo "  (无 csv)"
echo ""
echo "图片文件列表:"
ls -lh "${OUTDIR}/plots/"* 2>/dev/null || echo "  (无图片)"
