library(Seurat)
library(ggplot2)
library(gplots)
library(Rtsne)
library(ape)
library(tidytree)
library(ggtree)

setwd('~/Documents/Code/R/06_JunyeChen/code_submit/Fig_3M/')

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

table(seu$anno)
seu$marker <- ifelse(
  seu$tissue == "PLA", "PECAM1",
  ifelse(seu$RNA_snn_res.1.2 == "9", "CXCL14", "CD55")
)
table(seu$marker)
table(seu$patient)


patient_list <- c("patient_1", "patient_3", "patient_4")

setwd("example_step1/")

for (i in patient_list) {
  message("Processing: ", i)
  
  
  df <- readRDS(paste0("../", i, "_master.rds"))
  colnames(df) <- gsub("cell_", "", colnames(df))
  sds <- apply(df, 1, sd)
  df  <- df[sds > 0.3, ]
  rm(sds)
  
  
  pdf(file = paste0(i, "_heatmap_sd_0.3.pdf"), width = 30, height = 10)
  heatmap.2(
    as.matrix(df),
    col    = colorRampPalette(c("blue", "white", "red"))(256),
    key    = TRUE,
    symkey = FALSE,
    density.info = "none",
    trace        = "none",
    scale        = "row",
    cexRow       = 0.4,
    cexCol       = 0.8
  )
  dev.off()
  
  
  set.seed(123)
  t  <- Rtsne(t(df), perplexity = 5, check_duplicates = FALSE, normalize = FALSE)
  ty <- as.data.frame(t$Y)
  rownames(ty) <- colnames(df)
  colnames(ty) <- c("TSNE_1", "TSNE_2")
  write.csv(ty, paste0(i, "_tsne.csv"))
  seu <- AddMetaData(seu, metadata = ty)
  p_tsne <- ggplot(
    data = subset(seu@meta.data, rownames(seu@meta.data) %in% rownames(ty)),
    aes(x = TSNE_1, y = TSNE_2, color = marker)
  ) +
    geom_point(size = 0.6) +
    scale_color_manual(values = marker_cols) +
    theme_bw() +
    theme(
      panel.border     = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line        = element_line(color = "black")
    )
  ggsave(paste0(i, "_tsne_sd_0.3.pdf"), p_tsne, width = 6, height = 5)
  rm(t, p_tsne)
  
  
  hc        <- hclust(dist(t(df)))
  tree_root <- as.phylo(hc)
  tip_names <- tree_root$tip.label
  meta_idx  <- match(tip_names, seu$name)
  group <- data.frame(
    name   = tip_names,
    tissue = seu$tissue[meta_idx],
    marker = seu$marker[meta_idx],
    row.names = tip_names,
    stringsAsFactors = FALSE
  )
  p_base_tissue <- ggtree(
    tree_root,
    layout    = "circular",
    ladderize = FALSE,
    size      = 0.1,
    color     = "grey80"
  )
  ann_tissue <- data.frame(tissue = group$tissue)
  rownames(ann_tissue) <- rownames(group)  

  ann_marker <- data.frame(marker = group$marker)
  rownames(ann_marker) <- rownames(group)

  p_double <- gheatmap(
    p_base_tissue,
    ann_tissue,
    width  = 0.05,
    offset = 0,
    colnames = FALSE
  )
  p_double <- gheatmap(
    p_double,
    ann_marker,
    width  = 0.05,
    offset = 0.12,  
    colnames = FALSE
  ) +
    scale_fill_manual(values = c(tissue_cols, marker_cols), name = NULL)

  ggsave(paste0(i, "_tree_tissue_marker_double_ring.pdf"),
         p_double, width = 7, height = 7)
  
  rm(df, ty, hc, tree_raw, tree_root,
     tip_names, meta_idx, group,
     p_base_tissue, p_tree_tissue,
     p_base_marker, p_tree_marker,
     ann_tissue, ann_marker)
}
