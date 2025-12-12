library(monocle)
library(Seurat)
library(dplyr)
library(ggplot2)
library(scales)
library(tidyr)
library(ggpubr)
library(ggsci)
library(scCustomize)
library(tidyverse)
library(RColorBrewer)
library(harmony)
library(showtext)
library(plyr)
library(qs)
library(SeuratWrappers)

pal = DiscretePalette_scCustomize(num_colors = 12, palette = "ditto_seq")[-c(8)]


total_ADSC <- readRDS("total_ADSC.rds")
total_pla <- qread("total_pla.qs")
total_endo <- subset(total_pla,subset = cluster_name_slim == "Endothelial cell")

total_ADSC_sub <- subset(total_ADSC,subset = orig.ident == "CAS6")
total_endo_sub <- subset(total_endo,subset = sample == "Patient2")

data <- merge(total_ADSC_sub,total_endo_sub)
data$cluster_name_slim <- factor(data$cluster_name_slim,levels = c("ADSC_CD55","ADSC_CXCL14","Endothelial cell"))

Idents(data) <- "cluster_name_slim"
data <- subset(data,downsample=200)

exp.matrix <- as.sparse(data@assays$RNA@counts)
feature_ann <- data.frame(gene_id=rownames(exp.matrix),gene_short_name=rownames(exp.matrix))
rownames(feature_ann) <- rownames(data@assays$RNA@counts)
fd <- new("AnnotatedDataFrame", data = feature_ann)
sample_ann <- data@meta.data
rownames(sample_ann) <- colnames(exp.matrix)
pd <- new("AnnotatedDataFrame", data = sample_ann)

cds <- newCellDataSet(exp.matrix, phenoData =pd, featureData =fd, expressionFamily=negbinomial.size())

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(cds),  num_cells_expressed >= 10))
disp_table <- dispersionTable(cds)
ordering_genes <- as.character(subset(disp_table,mean_expression >= 0.3 & dispersion_empirical >= dispersion_fit)$gene_id)
cds <- setOrderingFilter(cds, ordering_genes)
cds <- reduceDimension(cds, max_components = 2, reduction_method = 'DDRTree',residualModelFormulaStr = "~sample")
cds1 <- orderCells(cds1,reverse = F)

pdf("FigureS5_Pseudotime_cell_name_P3.pdf",width = 5,height = 5)
plot_cell_trajectory(cds3, color_by = "cluster_name_slim",show_branch_points = F,cell_size = 1)+
  scale_color_manual(values = pal[c(3,1,5)])
dev.off()

pdf("FigureS5_Pseudotime_time_P3.pdf",width = 5,height = 5)
plot_cell_trajectory(cds3, color_by = "Pseudotime",show_branch_points = F,cell_size = 1)
dev.off()

cds_list <- list(cds1,cds2,cds3)
names(cds_list) <- c("P1","P2","P3")

qsave(cds_list,file = "cds_list.qs")
