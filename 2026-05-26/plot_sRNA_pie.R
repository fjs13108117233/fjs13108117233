#!/usr/bin/env Rscript
# =============================================================================
# plot_sRNA_pie.R - 各类 sRNA 数量饼图 / 环形图 (PPT 用)
# =============================================================================
# 风格: Nature (NPG) 配色 + Arial + 大字体, 适合 PPT 展示
#
# 输入: all_sRNA.bed (build_sRNA_bed.sh 输出, 首行 # 开头, 含 sRNA_type 列)
# 输出:
#   sRNA_pie.pdf/png     - 饼图 (扇区标百分比, 图例标类型+数量)
#   sRNA_donut.pdf/png   - 环形图 (中心显示总数)
#   sRNA_type_bar.pdf/png- 横向柱状图 (备选)
#
# 用法:
#   Rscript plot_sRNA_pie.R --bed all_sRNA.bed --outdir plots
# =============================================================================

suppressMessages({
    library(optparse)
    library(ggplot2)
})

option_list <- list(
    make_option("--bed", type = "character",
                default = "/public/home/h14166/fang/Heyufei/smrna_seq/all_sRNA.bed",
                help = "all_sRNA.bed 路径"),
    make_option("--outdir", type = "character", default = "plots",
                help = "输出目录"),
    make_option("--min-samples", type = "integer", default = 1,
                help = "仅统计 sample_count >= N 的 sRNA (默认 1=全部)")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

FONT <- "Arial"

# NPG 配色 (Nature Publishing Group)
NPG <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
         "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85")

save_plot <- function(base, plot, width, height) {
    ggsave(file.path(opt$outdir, paste0(base, ".pdf")), plot,
           width = width, height = height, device = cairo_pdf, limitsize = FALSE)
    ggsave(file.path(opt$outdir, paste0(base, ".png")), plot,
           width = width, height = height, dpi = 300, limitsize = FALSE)
}

# ===== 读取 all_sRNA.bed =====
message(">>> 读取 ", opt$bed)
hdr <- readLines(opt$bed, n = 1)
cols <- strsplit(sub("^#", "", hdr), "\t")[[1]]
df <- read.delim(opt$bed, header = FALSE, comment.char = "#", stringsAsFactors = FALSE)
colnames(df) <- cols[seq_len(ncol(df))]

if (!"sRNA_type" %in% colnames(df)) {
    stop("[ERROR] 未找到 sRNA_type 列。实际列: ", paste(colnames(df), collapse = ", "))
}

# 可选: 按 sample_count 过滤
if (opt$`min-samples` > 1 && "sample_count" %in% colnames(df)) {
    df <- df[as.integer(df$sample_count) >= opt$`min-samples`, ]
    message("  按 sample_count >= ", opt$`min-samples`, " 过滤后: ", nrow(df))
}

message("  sRNA 总数: ", nrow(df))

# ===== 统计各类型数量 =====
tab <- as.data.frame(table(df$sRNA_type), stringsAsFactors = FALSE)
colnames(tab) <- c("type", "count")
tab <- tab[order(-tab$count), ]
tab$pct <- tab$count / sum(tab$count) * 100
tab$label <- paste0(tab$type, "\n", tab$count, " (", sprintf("%.1f%%", tab$pct), ")")
tab$legend_label <- paste0(tab$type, "  (n=", tab$count, ", ", sprintf("%.1f%%", tab$pct), ")")
tab$type <- factor(tab$type, levels = tab$type)  # 按数量排序固定

total_n <- sum(tab$count)
n_type <- nrow(tab)
pal <- setNames(NPG[seq_len(n_type)], levels(tab$type))

message("  类型数: ", n_type)
for (i in seq_len(nrow(tab))) {
    message(sprintf("    %-18s %5d (%.1f%%)", tab$type[i], tab$count[i], tab$pct[i]))
}

# 大字体配置 (PPT 用)
PPT_TITLE <- 30
PPT_TEXT  <- 24

# =============================================================================
# 图1: 饼图
# =============================================================================
message(">>> [1] 饼图 ...")
# 计算标签位置 (累积角度中点)
tab2 <- tab[order(tab$type), ]
tab2$ymax <- cumsum(tab2$count)
tab2$ymin <- c(0, head(tab2$ymax, -1))
tab2$ymid <- (tab2$ymax + tab2$ymin) / 2
# 仅对占比 >= 3% 的扇区在图内标百分比 (避免重叠)
tab2$inlab <- ifelse(tab2$pct >= 3, sprintf("%.1f%%", tab2$pct), "")

p_pie <- ggplot(tab2, aes(x = "", y = count, fill = type)) +
    geom_col(width = 1, colour = "white", linewidth = 1.2) +
    coord_polar(theta = "y", start = 0) +
    geom_text(aes(y = ymid, label = inlab),
              colour = "white", fontface = "bold", family = FONT, size = 7) +
    scale_fill_manual(values = pal, labels = tab$legend_label, name = "sRNA Type") +
    labs(title = paste0("sRNA Classification (n = ", total_n, ")")) +
    theme_void(base_family = FONT) +
    theme(
        plot.title   = element_text(family = FONT, face = "bold", size = PPT_TITLE,
                                    hjust = 0.5, margin = margin(b = 10)),
        legend.title = element_text(family = FONT, face = "bold", size = PPT_TEXT),
        legend.text  = element_text(family = FONT, face = "bold", size = PPT_TEXT - 4),
        legend.key.size = unit(1, "cm"),
        plot.margin  = margin(15, 15, 15, 15)
    )
save_plot("sRNA_pie", p_pie, 14, 9)

# =============================================================================
# 图2: 环形图 (donut, 中心显示总数)
# =============================================================================
message(">>> [2] 环形图 ...")
p_donut <- ggplot(tab2, aes(x = 2, y = count, fill = type)) +
    geom_col(width = 1, colour = "white", linewidth = 1.2) +
    coord_polar(theta = "y", start = 0) +
    xlim(0.5, 2.5) +
    geom_text(aes(y = ymid, label = inlab),
              colour = "white", fontface = "bold", family = FONT, size = 7) +
    annotate("text", x = 0.5, y = 0, label = paste0("Total\n", total_n),
             family = FONT, fontface = "bold", size = 11) +
    scale_fill_manual(values = pal, labels = tab$legend_label, name = "sRNA Type") +
    labs(title = paste0("sRNA Classification (n = ", total_n, ")")) +
    theme_void(base_family = FONT) +
    theme(
        plot.title   = element_text(family = FONT, face = "bold", size = PPT_TITLE,
                                    hjust = 0.5, margin = margin(b = 10)),
        legend.title = element_text(family = FONT, face = "bold", size = PPT_TEXT),
        legend.text  = element_text(family = FONT, face = "bold", size = PPT_TEXT - 4),
        legend.key.size = unit(1, "cm"),
        plot.margin  = margin(15, 15, 15, 15)
    )
save_plot("sRNA_donut", p_donut, 14, 9)

# =============================================================================
# 图3: 横向柱状图 (备选, 数量一目了然)
# =============================================================================
message(">>> [3] 横向柱状图 ...")
p_bar <- ggplot(tab, aes(x = reorder(type, count), y = count, fill = type)) +
    geom_col(width = 0.75, colour = "black", linewidth = 0.8) +
    geom_text(aes(label = paste0(count, " (", sprintf("%.1f%%", pct), ")")),
              hjust = -0.1, family = FONT, fontface = "bold", size = 7) +
    coord_flip(clip = "off") +
    scale_fill_manual(values = pal) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(title = paste0("sRNA Classification (n = ", total_n, ")"),
         x = NULL, y = "Count") +
    theme_bw(base_size = PPT_TEXT, base_family = FONT) +
    theme(
        text         = element_text(family = FONT, face = "bold", colour = "black"),
        plot.title   = element_text(family = FONT, face = "bold", size = PPT_TITLE, hjust = 0.5),
        axis.title   = element_text(family = FONT, face = "bold", size = PPT_TEXT),
        axis.text    = element_text(family = FONT, face = "bold", size = PPT_TEXT - 2,
                                    colour = "black"),
        axis.ticks   = element_line(colour = "black", linewidth = 1),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
        panel.grid   = element_blank(),
        legend.position = "none",
        plot.margin  = margin(15, 30, 15, 15)
    )
save_plot("sRNA_type_bar", p_bar, 13, 8)

# 同时保存统计表
write.csv(tab[, c("type", "count", "pct")],
          file.path(opt$outdir, "sRNA_type_counts.csv"), row.names = FALSE)

message("\n=== 完成 ===")
message("输出: ", opt$outdir)
message("  sRNA_pie         - 饼图")
message("  sRNA_donut       - 环形图")
message("  sRNA_type_bar    - 横向柱状图")
message("  sRNA_type_counts.csv - 统计表")
