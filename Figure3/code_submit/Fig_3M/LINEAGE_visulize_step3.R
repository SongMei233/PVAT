library(Seurat)
library(ggplot2)
library(gplots)
library(Rtsne)
library(ape)
library(tidytree)
library(ggtree)

setwd('~/Documents/Code/R/06_JunyeChen/code_submit/Fig_3M/')

# 1. 读入 Seurat，对象合并
seu <- readRDS("example.rds")

marker_cols <- c(
  PECAM1 = "#619cff",
  CXCL14 = "#00ba38",
  CD55   = "#f8766d"
)

tissue_cols <- c(
  PLA  = "#66C5CC",
  PVAT = "#F6CF71"
)

label_palette_base <- c(
  "#F89C74", "#DCB0F2", "#87C55F","#9EB9F3", "#FE88B1", "#C9DB74", "#8BE0A4"
)

seu$marker <- ifelse(
  seu$tissue == "PLA", "PECAM1",
  ifelse(seu$RNA_snn_res.1.2 == "9", "CXCL14", "CD55")
)
table(seu$marker)
table(seu$patient)
seu$marker <- factor(seu$marker, levels = names(marker_cols))
patient_list <- c("patient_1", "patient_3", "patient_4")
setwd("example_step3/")

for (i in patient_list) {
  message("Processing: ", i)
  
  
  df <- readRDS(paste0("../", i, "_master.rds"))
  colnames(df) <- gsub("cell_", "", colnames(df))
  sds <- apply(df, 1, sd)
  df  <- df[sds > 0.3, ]
  rm(sds)
  
  
  ty <- read.csv(
    paste0("../example_step2/", i, "_tsne_label.csv"),
    row.names = 1
  )
  label_levels <- sort(unique(ty$label))
  label_cols   <- label_palette_base[seq_along(label_levels)]
  names(label_cols) <- label_levels
  ty$label <- factor(ty$label, levels = label_levels)
  seu <- AddMetaData(seu, metadata = ty)
  cells_use <- rownames(ty)
  plot_meta <- subset(seu@meta.data, rownames(seu@meta.data) %in% cells_use)
  plot_meta$marker <- factor(plot_meta$marker, levels = names(marker_cols))
  plot_meta$label  <- factor(plot_meta$label,  levels = label_levels)
  p_tsne_label <- ggplot(plot_meta, aes(x = TSNE_1, y = TSNE_2)) +
    ggforce::geom_mark_hull(
      aes(
        fill  = label,
        label = label
      ),
      alpha  = 0.15,                    
      color  = "grey40",                 
      expand = unit(0.8, "mm"),
      radius = unit(0.1, "mm")
    ) +
    geom_point(aes(color = marker), size = 0.6) +
    scale_fill_manual(values = label_cols, name = "label") +
    scale_color_manual(values = marker_cols, name = "marker") +
    theme_bw() +
    theme(
      panel.border     = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line        = element_line(color = "black")
    )
  ggsave(
    paste0(i, "_tsne_sd_0.3_label.pdf"),
    p_tsne_label,
    width  = 7,
    height = 6
  )
  
  
  df1 <- df[, rownames(ty)]
  colnames(df1) <- paste0(rownames(ty), "__", ty$label)
  plot_color <- label_cols[as.character(ty$label)]
  names(plot_color) <- rownames(ty)
  group_labels <- label_levels
  group_colors <- label_cols
  pdf(file = paste0(i, "_heatmap_sd_0.3_label.pdf"), width = 50, height = 10)
  heatmap.2(
    as.matrix(df1),
    key           = TRUE,
    ColSideColors = plot_color,
    col           = colorRampPalette(c("blue","white","red"))(20),
    breaks        = seq(-1, 1, 0.1),
    symkey        = TRUE,
    density.info  = "none",
    trace         = "none",
    scale         = "row",
    cexRow        = 0.4,
    cexCol        = 0.8,
    margins       = c(10, 10)
  )
  legend(
    "topright",
    legend = group_labels,
    fill   = group_colors,
    border = "black",
    bty    = "n",
    title  = "Groups",
    cex    = 2
  )
  dev.off()
  
  
  rm(df1, plot_color, group_labels, group_colors,
     df, ty, cells_use, plot_meta, p_tsne_label, label_levels, label_cols)
}
