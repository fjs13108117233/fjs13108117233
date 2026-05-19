#!/bin/bash
#SBATCH -p com300
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH -J dmr_srna
#SBATCH --output logs/dmr_srna_%j.out
#SBATCH --error  logs/dmr_srna_%j.err

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
mkdir -p logs

source "$HOME/anaconda3/etc/profile.d/conda.sh"
conda activate RNA

# =========================
# 路径
# =========================
SRNA_BED="/public/home/h14166/fang/Heyufei/smrna_seq/all_sRNA.bed"
DMR_DIR="/public/home/h14166/fang/Heyufei/Nanopore/filter_location"
OUTDIR="/public/home/h14166/fang/Heyufei/Nanopore/srna_overlap"

mkdir -p "$OUTDIR"

# 4 种修饰
MODS=("4mc" "5hmc" "5mc" "6ma")

# =========================
# 准备 sRNA BED
# =========================
echo "=== DMR vs sRNA 交集分析: $(date) ==="
echo "sRNA BED: $SRNA_BED"
echo ""

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT


# 去表头, 提取 chrom start end sRNA_name sRNA_type
tail -n +2 "$SRNA_BED" | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5}' \
    | sort -k1,1 -k2,2n > "${TMP}/srna.bed"

SRNA_TOTAL=$(wc -l < "${TMP}/srna.bed")
echo "sRNA loci 总数: $SRNA_TOTAL"
echo ""

# =========================
# 汇总表头
# =========================
SUMMARY="${OUTDIR}/overlap_summary.tsv"
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "mod_type" "region" "DMR_count" "overlap_DMR" "pct_DMR" \
    "overlap_sRNA" "pct_sRNA" "sRNA_type_top3" \
    > "$SUMMARY"

# =========================
# 逐修饰逐文件处理
# =========================
for MOD in "${MODS[@]}"; do
    echo "============================================"
    echo ">>> 修饰: $MOD"

    MOD_DIR="${DMR_DIR}/${MOD}"
    if [ ! -d "$MOD_DIR" ]; then
        echo "  [WARNING] 目录不存在: $MOD_DIR"
        continue
    fi

    for DMR_FILE in "${MOD_DIR}"/*_afterselect_*.tsv; do
        [ -f "$DMR_FILE" ] || continue

        # 判断是 genic 还是 intergenic
        FNAME=$(basename "$DMR_FILE")
        if echo "$FNAME" | grep -q "intergenic"; then
            REGION="intergenic"
        else
            REGION="genic"
        fi

        echo ""
        echo "  --- $MOD / $REGION ---"
        echo "  文件: $FNAME"


        # DMR 转 BED (去表头, 取 chrom start end gene)
        tail -n +2 "$DMR_FILE" | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$5}' \
            | sort -k1,1 -k2,2n > "${TMP}/dmr.bed"

        DMR_COUNT=$(wc -l < "${TMP}/dmr.bed")
        echo "  DMR 位点数: $DMR_COUNT"

        if [ "$DMR_COUNT" -eq 0 ]; then
            echo "  [SKIP] 无 DMR"
            printf "%s\t%s\t%d\t%d\t%.2f\t%d\t%.2f\t%s\n" \
                "$MOD" "$REGION" 0 0 0 0 0 "NA" >> "$SUMMARY"
            continue
        fi

        # --- 1. 有多少 DMR 与 sRNA 重叠 ---
        OVERLAP_DMR=$(bedtools intersect \
            -a "${TMP}/dmr.bed" -b "${TMP}/srna.bed" -u | wc -l)
        PCT_DMR=$(awk "BEGIN {printf \"%.2f\", $OVERLAP_DMR/$DMR_COUNT*100}")

        # --- 2. 有多少 sRNA 被 DMR 命中 ---
        OVERLAP_SRNA=$(bedtools intersect \
            -a "${TMP}/srna.bed" -b "${TMP}/dmr.bed" -u | wc -l)
        PCT_SRNA=$(awk "BEGIN {printf \"%.2f\", $OVERLAP_SRNA/$SRNA_TOTAL*100}")

        echo "  与 sRNA 重叠的 DMR: $OVERLAP_DMR / $DMR_COUNT ($PCT_DMR%)"
        echo "  被 DMR 命中的 sRNA: $OVERLAP_SRNA / $SRNA_TOTAL ($PCT_SRNA%)"


        # --- 3. 详细交集输出 (DMR + sRNA 完整信息) ---
        OUT_DETAIL="${OUTDIR}/${MOD}_${REGION}_x_sRNA_detail.tsv"

        # 表头
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "DMR_chrom" "DMR_start" "DMR_end" "DMR_gene" \
            "sRNA_chrom" "sRNA_start" "sRNA_end" "sRNA_name" "sRNA_type" \
            > "$OUT_DETAIL"

        bedtools intersect \
            -a "${TMP}/dmr.bed" \
            -b "${TMP}/srna.bed" \
            -wa -wb \
            >> "$OUT_DETAIL"

        echo "  详细结果: $OUT_DETAIL"

        # --- 4. sRNA 类型分布 ---
        echo "  sRNA 类型分布:"
        TYPE_DIST=$(tail -n +2 "$OUT_DETAIL" \
            | awk '{print $9}' | sort | uniq -c | sort -rn)
        echo "$TYPE_DIST" | head -5 | while read cnt tp; do
            printf "      %-20s %s\n" "$tp" "$cnt"
        done

        # 取 top3 类型写入汇总
        TOP3=$(echo "$TYPE_DIST" | head -3 \
            | awk '{printf "%s(%s)", $2, $1; if(NR<3) printf ","}')
        [ -z "$TOP3" ] && TOP3="NA"


        # --- 5. 写入汇总 ---
        printf "%s\t%s\t%d\t%d\t%s\t%d\t%s\t%s\n" \
            "$MOD" "$REGION" "$DMR_COUNT" "$OVERLAP_DMR" "$PCT_DMR" \
            "$OVERLAP_SRNA" "$PCT_SRNA" "$TOP3" \
            >> "$SUMMARY"

    done
done

# =========================
# 最终汇总
# =========================
echo ""
echo "============================================"
echo "============================================"
echo ""
echo ">>> 汇总表: $SUMMARY"
echo ""
column -t -s "$(printf '\t')" "$SUMMARY"
echo ""
echo ">>> 所有详细结果在: $OUTDIR/"
ls -lh "$OUTDIR"/*.tsv
echo ""
echo "=== 完成: $(date) ==="
