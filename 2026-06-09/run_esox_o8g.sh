#!/usr/bin/env bash
###############################################################################
# esox 8-oxo-dG (O8G) detection pipeline for 5 samples
#
# Inputs already prepared:
#   - raw fast5 per sample : ${FAST5_ROOT}/<sample>/*.fast5
#   - Dorado fastq per sample (big): ${FQ_ROOT}/<sample>.fastq
#
# Steps per sample:
#   0) split the big fastq -> per-fast5 matching fastq (single streaming pass)
#   1) basecall.py : raw fast5 + matching fastq -> esox .fastq + .npz
#   2) modcall.py  : .npz/.fastq -> per-G 8-oxo-dG scores (.txt)
#   3) summarize.py: aggregate O8G calls per sample
#
# NOTE: data is large (~100GB fastq/sample). Run on a GPU node. Steps 1-2 are
#       the heavy GPU steps; consider submitting one job per sample.
###############################################################################
set -euo pipefail

# ---------------------------------------------------------------------------
ESOX_DIR="/public/home/h14166/fang/Heyufei/O8G/tools/esox"
FAST5_ROOT="/public/home/h14166/fang/Heyufei/O8G/fast5"     # contains the 5 sample subdirs
FQ_ROOT="/public/home/h14166/fang/Heyufei/O8G/fq"           # contains <sample>.fastq
OUT_ROOT="/public/home/h14166/fang/Heyufei/O8G/esox_out"
DEVICE="cuda:0"            # or "cpu"
THRESHOLD="0.95"
MAX_OPEN="512"            # max simultaneously open files during fastq split
# ---------------------------------------------------------------------------

# Run all 5 samples by default, or a single sample if passed as $1:
#   bash run_esox_o8g.sh            # all samples
#   bash run_esox_o8g.sh OE1        # just OE1 (good for per-sample job submission)
if [ "$#" -ge 1 ]; then
    SAMPLES=("$1")
else
    SAMPLES=(OE1 RDM_-6bp_SAM_DNA RDM_-8bp_SAM_DNA RDM_OE2_SAM_DNA WT_SAM_DNA)
fi

BONITO_MODEL="${ESOX_DIR}/static/models/bonito.pt"
REMORA_MODEL="${ESOX_DIR}/static/models/remora.pt"

cd "${ESOX_DIR}"
mkdir -p "${OUT_ROOT}"

for S in "${SAMPLES[@]}"; do
    echo "==================================================================="
    echo "Sample: ${S}"
    echo "==================================================================="

    FAST5_DIR="${FAST5_ROOT}/${S}"
    BIG_FQ="${FQ_ROOT}/${S}.fastq"
    MATCHED_FQ="${OUT_ROOT}/${S}/matched_fastq"
    BC_OUT="${OUT_ROOT}/${S}/basecall_out"
    MC_OUT="${OUT_ROOT}/${S}/modcall_out"

    mkdir -p "${MATCHED_FQ}" "${BC_OUT}" "${MC_OUT}"

    # 0) split big fastq -> per-fast5 fastq (streaming, low memory)
    echo "[0/3] splitting fastq to match fast5..."
    python3 scripts/split_fastq_by_fast5.py \
        --fast5-path  "${FAST5_DIR}" \
        --guppy-fastq "${BIG_FQ}" \
        --output-path "${MATCHED_FQ}" \
        --max-open "${MAX_OPEN}"

    # 1) esox basecalling -> .fastq + .npz
    echo "[1/3] basecalling (esox bonito)..."
    python3 scripts/basecall.py \
        --fast5-path  "${FAST5_DIR}" \
        --fastq-path  "${MATCHED_FQ}" \
        --output-path "${BC_OUT}" \
        --model-file  "${BONITO_MODEL}" \
        --progress-bar \
        --device "${DEVICE}"

    # 2) modification calling -> per-G oxog_score .txt
    echo "[2/3] modification calling (esox remora)..."
    python3 scripts/modcall.py \
        --input-path  "${BC_OUT}" \
        --output-path "${MC_OUT}" \
        --model-file  "${REMORA_MODEL}" \
        --progress-bar \
        --device "${DEVICE}"

    # 3) summarize this sample
    echo "[3/3] summarizing..."
    python3 scripts/summarize.py \
        --modcall-path "${MC_OUT}" \
        --threshold "${THRESHOLD}" \
        --sample "${S}" \
        --output "${OUT_ROOT}/summary.tsv"
done

echo "ALL DONE. Per-sample summary -> ${OUT_ROOT}/summary.tsv"
