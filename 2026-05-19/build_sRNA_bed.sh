#!/bin/bash
#SBATCH -p com300
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G
#SBATCH -J build_bed
#SBATCH --output logs/build_bed_%j.out
#SBATCH --error  logs/build_bed_%j.err

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"

mkdir -p logs

source "$HOME/anaconda3/etc/profile.d/conda.sh"
conda activate RNA

# =========================
# 用户参数
# =========================
WORKDIR="/public/home/h14166/fang/Heyufei/smrna_seq"
LIST="${WORKDIR}/list"
REF_MIRNA_GFF="/public/home/h14166/refer/Zea_mays.miRNA.gff3"
OUTBED="${WORKDIR}/all_sRNA.bed"

# Python 分类脚本路径 (与本脚本同目录)
CLASSIFY_PY="${SLURM_SUBMIT_DIR}/classify_sRNA.py"

# 可选注释文件
REF_TE_BED="${WORKDIR}/Zea_mays.TE.bed"
REF_TAS_BED="${WORKDIR}/Zea_mays.TAS_loci.bed"

# PhaseScore 阈值
PHASE_CUTOFF_21=1.31
PHASE_CUTOFF_24=2.0
MIN_PHASE_SUPPORT=1

# =========================
# 检查
# =========================
command -v bedtools >/dev/null 2>&1 || { echo "[ERROR] bedtools not found" >&2; exit 1; }
command -v python3 >/dev/null 2>&1 || { echo "[ERROR] python3 not found" >&2; exit 1; }

[ -f "$LIST" ]          || { echo "[ERROR] 样本列表不存在: $LIST" >&2; exit 1; }
[ -f "$REF_MIRNA_GFF" ] || { echo "[ERROR] miRNA GFF3 不存在: $REF_MIRNA_GFF" >&2; exit 1; }
[ -f "$CLASSIFY_PY" ]   || { echo "[ERROR] 分类脚本不存在: $CLASSIFY_PY" >&2; exit 1; }

HAS_TE="no"
[ -f "$REF_TE_BED" ] && HAS_TE="yes" && echo "[INFO] TE 注释: $REF_TE_BED"

HAS_TAS="no"
[ -f "$REF_TAS_BED" ] && HAS_TAS="yes" && echo "[INFO] TAS 注释: $REF_TAS_BED"

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

echo "=== build_sRNA_bed.sh 开始: $(date) ==="
echo "工作目录: $WORKDIR"
echo "输出文件: $OUTBED"
echo "PHASE_CUTOFF_21=$PHASE_CUTOFF_21  PHASE_CUTOFF_24=$PHASE_CUTOFF_24  MIN_PHASE_SUPPORT=$MIN_PHASE_SUPPORT"
echo ""

# =========================
# 1. 解析 ShortStack Results.txt
# =========================
echo ">>> [1/6] 解析 ShortStack Results.txt ..."

ALL_BED="${TMP}/all_loci.bed"
> "$ALL_BED"

while IFS= read -r SAMPLE; do
    SAMPLE="$(echo "$SAMPLE" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')"
    [ -z "$SAMPLE" ] && continue
    [[ "$SAMPLE" =~ ^# ]] && continue

    R="${WORKDIR}/${SAMPLE}/shortstack_out/Results.txt"
    if [ ! -f "$R" ]; then
        echo "  [WARNING] 缺少: $R" >&2
        continue
    fi
    echo "  样本: $SAMPLE"

    awk -v S="$SAMPLE" '
    BEGIN { FS = OFS = "\t" }
    function get_col(name) { return (name in idx) ? $(idx[name]) : "" }
    NR == 1 {
        for (i = 1; i <= NF; i++) idx[$i] = i
        next
    }
    {
        loc = get_col("#Locus")
        if (loc == "") loc = get_col("Locus")
        if (loc == "") next

        chrom = loc; sub(/:[0-9]+-[0-9]+$/, "", chrom)
        coord = loc; sub(/^.*:/, "", coord)
        if (split(coord, a, "-") != 2) next
        start = a[1] - 1; end = a[2]
        if (start < 0) start = 0
        if (end <= start) next

        dicer  = get_col("DicerCall")
        phase  = get_col("PhaseScore")
        mirna  = get_col("MIRNA")
        strand = get_col("Strand")

        if (dicer == "" || dicer == "." || dicer == "NA") dicer = "N"
        if (phase == "" || phase == "." || phase == "NA") phase = "0"
        if (mirna == "" || mirna == "." || mirna == "NA") mirna = "N"
        if (strand == "" || strand == "." || strand == "NA") strand = "."

        print chrom, start, end, S, dicer, phase, mirna, strand
    }
    ' "$R" >> "$ALL_BED"
done < "$LIST"

RAW_N=$(wc -l < "$ALL_BED")
echo "  原始 locus 数: $RAW_N"
[ "$RAW_N" -eq 0 ] && { echo "[ERROR] 无 locus" >&2; exit 1; }

# =========================
# 2. merge
# =========================
echo ""
echo ">>> [2/6] merge ..."

sort -k1,1 -k2,2n "$ALL_BED" > "${TMP}/sorted.bed"
bedtools merge -i "${TMP}/sorted.bed" -d 0 \
    -c 4,5,6,7,8 -o collapse,collapse,collapse,collapse,collapse \
    > "${TMP}/merged.bed"

awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"sRNA_"NR,$4,$5,$6,$7,$8}' \
    "${TMP}/merged.bed" > "${TMP}/with_id.bed"

echo "  合并后: $(wc -l < "${TMP}/with_id.bed")"

# =========================
# 3. intersect miRNA
# =========================
echo ""
echo ">>> [3/6] intersect miRNA ..."

awk 'BEGIN{OFS="\t"} !/^#/ && $3=="miRNA_primary_transcript" {print}' "$REF_MIRNA_GFF" \
    | sort -k1,1 -k4,4n > "${TMP}/ref.gff3"

if [ "$(wc -l < "${TMP}/ref.gff3")" -eq 0 ]; then
    awk 'BEGIN{OFS="\t"} !/^#/ && NF>=9 {print}' "$REF_MIRNA_GFF" \
        | sort -k1,1 -k4,4n > "${TMP}/ref.gff3"
fi

awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4}' "${TMP}/with_id.bed" > "${TMP}/bed4.bed"

bedtools intersect -a "${TMP}/bed4.bed" -b "${TMP}/ref.gff3" -wa -wb -loj \
    > "${TMP}/inter_mirna.tsv"

echo "  miRNA 记录数: $(wc -l < "${TMP}/ref.gff3")"

# =========================
# 4. intersect TE / TAS
# =========================
echo ""
echo ">>> [4/6] intersect TE / TAS ..."

if [ "$HAS_TE" = "yes" ]; then
    bedtools intersect -a "${TMP}/bed4.bed" -b "$REF_TE_BED" -wa -wb -loj \
        > "${TMP}/inter_te.tsv"
else
    > "${TMP}/inter_te.tsv"
    echo "  [SKIP] 无 TE 文件"
fi

if [ "$HAS_TAS" = "yes" ]; then
    bedtools intersect -a "${TMP}/bed4.bed" -b "$REF_TAS_BED" -wa -wb -loj \
        > "${TMP}/inter_tas.tsv"
else
    > "${TMP}/inter_tas.tsv"
    echo "  [SKIP] 无 TAS 文件"
fi

# =========================
# 5. Python 分类
# =========================
echo ""
echo ">>> [5/6] 分类注释 ..."

python3 "$CLASSIFY_PY" \
    "${TMP}/with_id.bed" \
    "${TMP}/inter_mirna.tsv" \
    "${TMP}/inter_te.tsv" \
    "${TMP}/inter_tas.tsv" \
    "$OUTBED" \
    "$PHASE_CUTOFF_21" \
    "$PHASE_CUTOFF_24" \
    "$MIN_PHASE_SUPPORT"

# =========================
# 6. 检查
# =========================
echo ""
echo ">>> [6/6] 检查输出 ..."

if [ ! -s "$OUTBED" ]; then
    echo "[ERROR] 输出为空: $OUTBED" >&2
    exit 1
fi

echo "  总行数: $(wc -l < "$OUTBED")"
echo ""
head -8 "$OUTBED" | column -t -s "$(printf '\t')" || head -8 "$OUTBED"

echo ""
echo "=== 完成: $(date) ==="
