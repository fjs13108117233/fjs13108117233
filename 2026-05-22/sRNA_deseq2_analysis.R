#!/usr/bin/env Rscript
# =============================================================================
# sRNA_deseq2_analysis.R
# =============================================================================
# DESeq2 差异分析：所有 sRNA（基于 build_sRNA_counts.py 输出）
#
# 用法:
#   Rscript sRNA_deseq2_analysis.R /path/to/sRNA_diff_analysis
#
# 输入 (在 WORKDIR/counts/ 下):
#   - sRNA_counts.csv        : count 矩阵 (sRNA x samples)
#   - sample_metadata.csv    : 样本分组
#   - sRNA_annotation.csv    : sRNA 注释 (类型、坐标等)
#
# 输出:
#   results/  : DE 结果表 (全部 + 显著)
#   plots/    : 火山图、PCA、热图、MA 图
# =============================================================================

suppressMessages({
    library(DESeq2)
    library(ggplot2)
    library(pheatmap)
    library(RColorBrewer)
})

# ===== 参数 =====
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("用法: Rscript sRNA_deseq2_analysis.R <WORKDIR>")
}
WORKDIR <- args[1]

# 阈值
PADJ_CUTOFF <- 0.05
LFC_CUTOFF  <- 1.0
MIN_COUNT   <- 10   # 行过滤: rowSums >= MIN_COUNT

# ===== 路径 =====
counts_file <- file.path(WORKDIR, "counts", "sRNA_counts.csv")
meta_file   <- file.path(WORKDIR, "counts", "sample_metadata.csv")
anno_file   <- file.path(WORKDIR, "counts", "sRNA_annotation.csv")
results_dir <- file.path(WORKDIR, "results")
plots_dir   <- file.path(WORKDIR, "plots")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir,   showWarnings = FALSE, recursive = TRUE)

# ===== 1. 读取数据 =====
message(">>> [1/6] 读取数据 ...")

counts  <- read.csv(counts_file, row.names = 1, check.names = FALSE)
coldata <- read.csv(meta_file, row.names = 1)
anno    <- read.csv(anno_file)

counts <- as.matrix(counts)
mode(counts) <- "integer"

# 确保样本顺序一致
common_samples <- intersect(colnames(counts), rownames(coldata))
if (length(common_samples) == 0) {
    stop("[ERROR] counts 与 metadata 无共同样本")
}
counts  <- counts[, common_samples, drop = FALSE]
coldata <- coldata[common_samples, , drop = FALSE]

message("  样本数: ", ncol(counts))
message("  sRNA 数: ", nrow(counts))
message("  分组: ", paste(unique(coldata$group), collapse = ", "))

# ===== 2. DESeq2 构建 =====
message(">>> [2/6] DESeq2 建模 ...")

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData   = coldata,
                              design    = ~ group)

# 预过滤低表达
keep <- rowSums(counts(dds)) >= MIN_COUNT
dds  <- dds[keep, ]
message("  过滤后 sRNA 数 (rowSums >= ", MIN_COUNT, "): ", nrow(dds))

# 设置 WT 为参考水平
dds$group <- relevel(dds$group, ref = "WT")

# 运行 DESeq2
dds <- DESeq(dds)

# ===== 3. 提取差异结果 =====
message(">>> [3/6] 提取差异结果 ...")

groups    <- levels(dds$group)
wt_group  <- "WT"
mt_groups <- setdiff(groups, wt_group)

summary_stats <- data.frame()
all_results   <- list()

for (mt in mt_groups) {
    message("  比较: ", mt, " vs ", wt_group)

    res <- results(dds, contrast = c("group", mt, wt_group),
                   alpha = PADJ_CUTOFF)
    res <- res[order(res$padj), ]

    df <- as.data.frame(res)
    df$sRNA_name   <- rownames(df)
    df$comparison  <- paste(mt, "vs", wt_group)

    # 合并注释
    if (nrow(anno) > 0 && "sRNA_name" %in% colnames(anno)) {
        df <- merge(df, anno[, c("sRNA_name", "sRNA_type", "chrom", "start", "end")],
                    by = "sRNA_name", all.x = TRUE)
    }

    # 显著基因
    sig <- subset(df, !is.na(padj) & padj < PADJ_CUTOFF & abs(log2FoldChange) > LFC_CUTOFF)
    up   <- sum(sig$log2FoldChange > 0, na.rm = TRUE)
    down <- sum(sig$log2FoldChange < 0, na.rm = TRUE)

    summary_stats <- rbind(summary_stats, data.frame(
        Comparison    = paste(mt, "vs", wt_group),
        Upregulated   = up,
        Downregulated = down,
        Total_DE      = up + down,
        stringsAsFactors = FALSE
    ))

    all_results[[mt]] <- df

    # 保存
    write.csv(df,  file.path(results_dir, paste0("DE_", mt, "_vs_WT.csv")),
              row.names = FALSE)
    write.csv(sig, file.path(results_dir, paste0("DE_", mt, "_vs_WT_significant.csv")),
              row.names = FALSE)

    message("    Up=", up, "  Down=", down, "  Total_DE=", up + down)
}

# 汇总表
write.csv(summary_stats, file.path(results_dir, "DE_summary_statistics.csv"),
          row.names = FALSE)
write.csv(do.call(rbind, all_results),
          file.path(results_dir, "all_comparisons.csv"), row.names = FALSE)

# ===== 4. 火山图 =====
message(">>> [4/6] 绘制火山图 ...")

for (mt in mt_groups) {
    df <- all_results[[mt]]
    df$sig <- ifelse(!is.na(df$padj) & df$padj < PADJ_CUTOFF & abs(df$log2FoldChange) > LFC_CUTOFF,
                     ifelse(df$log2FoldChange > 0, "Up", "Down"), "NS")
    df$sig <- factor(df$sig, levels = c("Up", "Down", "NS"))

    n_up   <- sum(df$sig == "Up", na.rm = TRUE)
    n_down <- sum(df$sig == "Down", na.rm = TRUE)

    p <- ggplot(df, aes(log2FoldChange, -log10(padj), color = sig)) +
        geom_point(alpha = 0.6, size = 1.2) +
        scale_color_manual(values = c(Up = "firebrick", Down = "dodgerblue", NS = "grey70"),
                           labels = c(paste0("Up (", n_up, ")"),
                                      paste0("Down (", n_down, ")"),
                                      "NS")) +
        geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = "dashed", color = "grey40") +
        geom_hline(yintercept = -log10(PADJ_CUTOFF), linetype = "dashed", color = "grey40") +
        labs(title = paste0(mt, " vs WT"),
             subtitle = paste0("padj < ", PADJ_CUTOFF, ", |log2FC| > ", LFC_CUTOFF),
             x = "log2 Fold Change", y = "-log10(padj)", color = "Status") +
        theme_bw(base_size = 12) +
        theme(legend.position = "top")

    ggsave(file.path(plots_dir, paste0("volcano_", mt, "_vs_WT.pdf")),
           p, width = 8, height = 6)
    ggsave(file.path(plots_dir, paste0("volcano_", mt, "_vs_WT.png")),
           p, width = 8, height = 6, dpi = 150)
}

# ===== 5. PCA =====
message(">>> [5/6] PCA 分析 ...")

vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

pcaData <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p <- ggplot(pcaData, aes(PC1, PC2, color = group, shape = group)) +
    geom_point(size = 4) +
    geom_text(aes(label = name), size = 2.5, vjust = -1.2, show.legend = FALSE) +
    labs(x = paste0("PC1: ", percentVar[1], "% variance"),
         y = paste0("PC2: ", percentVar[2], "% variance"),
         title = "PCA - sRNA expression",
         color = "Group", shape = "Group") +
    theme_bw(base_size = 12) +
    theme(legend.position = "right")

ggsave(file.path(plots_dir, "PCA_plot.pdf"), p, width = 9, height = 7)
ggsave(file.path(plots_dir, "PCA_plot.png"), p, width = 9, height = 7, dpi = 150)

# ===== 6. 样本距离热图 + DE 热图 =====
message(">>> [6/6] 热图 ...")

# 样本距离热图
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(vsd)
colnames(sampleDistMatrix) <- colnames(vsd)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

pdf(file.path(plots_dir, "sample_distance_heatmap.pdf"), width = 8, height = 7)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         main = "Sample-to-Sample Distance")
dev.off()

# 差异 sRNA 表达热图 (取所有比较中显著的 union)
all_sig_names <- unique(unlist(lapply(all_results, function(df) {
    subset(df, !is.na(padj) & padj < PADJ_CUTOFF & abs(log2FoldChange) > LFC_CUTOFF)$sRNA_name
})))

if (length(all_sig_names) > 0 && length(all_sig_names) <= 500) {
    sig_vsd <- assay(vsd)[rownames(assay(vsd)) %in% all_sig_names, , drop = FALSE]
    if (nrow(sig_vsd) > 1) {
        # z-score 标准化
        sig_scaled <- t(scale(t(sig_vsd)))

        # 注释列
        anno_col <- data.frame(Group = coldata[colnames(sig_vsd), "group"],
                               row.names = colnames(sig_vsd))

        pdf(file.path(plots_dir, "DE_sRNA_heatmap.pdf"),
            width = 10, height = min(20, 4 + nrow(sig_scaled) * 0.15))
        pheatmap(sig_scaled,
                 annotation_col = anno_col,
                 show_rownames = (nrow(sig_scaled) <= 80),
                 clustering_method = "ward.D2",
                 main = paste0("Differentially Expressed sRNAs (n=", nrow(sig_scaled), ")"),
                 fontsize_row = 6)
        dev.off()
        message("  DE 热图: ", nrow(sig_scaled), " sRNAs")
    }
} else if (length(all_sig_names) > 500) {
    message("  [SKIP] DE sRNA > 500, 跳过热图 (可取 top 子集)")
} else {
    message("  [SKIP] 无显著 DE sRNA")
}

# ===== 汇总 =====
message("\n=== sRNA 差异分析完成 ===")
message("结果目录: ", results_dir)
message("图片目录: ", plots_dir)
message("")
print(summary_stats)
