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

#define function of dimensionality reduction and cell clustering

Reduct <- function(total){
  DefaultAssay(total) <- "integrated"
  total <- RunPCA(total, verbose = FALSE)
  ElbowPlot(total, ndims = 50)
  total <- RunUMAP(total, dims = 1:30)
  total <- FindNeighbors(total, dims = 1:30)
  for(i in c(0.5,0.8,1,1.5)){
    total <- FindClusters(total, resolution = i)
  }
  return(total)
}

#define color palette

pal = DiscretePalette_scCustomize(num_colors = 12, palette = "ditto_seq")[-c(8)]

######################Macropahge##############################

total_immune <- readRDS("total_innate_immune.rds")

total_macro <- subset(total_immune,subset = cluster_name_slim == "Macrophage")

DefaultAssay(total_macro) <- "RNA"
total_macro[['SCT']] <- NULL
total_macro[['integrated']] <- NULL
total_macro <- NormalizeData(total_macro)
total_macro <- FindVariableFeatures(total_macro,selection.method = "vst", nfeatures = 1000)
total_macro <- ScaleData(total_macro)
total_macro <- RunPCA(total_macro, verbose = FALSE)

total_macro@meta.data$sample <- as.factor(total_macro@meta.data$sample)
total_macro <- RunHarmony(total_macro,group.by.vars="sample",plot_convergence = TRUE)
ElbowPlot(total_macro,ndims = 50)
total_macro <- RunUMAP(total_macro,reduction = "harmony", dims = 1:30)
total_macro <- FindNeighbors(total_macro,reduction = "harmony", dims = 1:30)

total_macro <- FindClusters(total_macro, resolution = 0.5)

DimPlot(total_macro, reduction = "umap", label=T,group.by = "RNA_snn_res.0.5")

DimPlot(total_macro,group.by = "sample")

markers <- c("PTPRC","MS4A1","CD19","BANK1","JCHAIN","SDC1","IGKC","CD3D","CD4","CCR7","SELL","CD8A","GZMK","CCL5","KLRD1","GNLY","NKG7","MKI67","STMN1","TUBA1B","IL1B","S100A8","S100A9","CD68","C1QA","LYVE1","CD1E","LAMP3","FSCN1","PTGDS","IL3RA","LILRA4","F3","CXCL14","CD55")

DotPlot(total_macro, features = markers, dot.scale = 6)  + theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) + RotatedAxis()

total_macro <- subset(total_macro,subset = RNA_snn_res.0.5 != "1")

markers <- FindAllMarkers(total_macro,logfc.threshold = 0.25,test.use = "wilcox",min.pct = 0.5,min.diff.pct = 0.2,only.pos = T,recorrect_umi = FALSE)

top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

clu14 <- FindMarkers(total_macro,ident.1 = c(1,4),logfc.threshold = 0.25,test.use = "wilcox",min.pct = 0.5,min.diff.pct = 0.2,only.pos = T,recorrect_umi = FALSE)

total_macro$cluster_name_slim <- "Macro_C1QA"
total_macro$cluster_name_slim[total_macro$RNA_snn_res.0.5 == "1"] <- "Macro_APOE"
total_macro$cluster_name_slim[total_macro$RNA_snn_res.0.5 == "2"] <- "Macro_CXCL8"
total_macro$cluster_name_slim[total_macro$RNA_snn_res.0.5 == "3"] <- "Macro_CCL3"

total_macro$cluster_name_slim <- factor(total_macro$cluster_name_slim,
                                        levels = c("Macro_C1QA","Macro_APOE","Macro_CXCL8","Macro_CCL3"))
Idents(total_macro) <- 'cluster_name_slim'

p <- DimPlot_scCustom(total_macro,colors_use = pal[5:8],figure_plot = T,label = F,label.size = 5,pt.size = 1);p

ggsave("Figure6G_Macrophage_marker_UMAP.pdf",p,width = 5,height = 4)

markers <- c("CD68","CD163","C1QA","APOE","CXCL8","CCL3")

p <- VlnPlot(total_macro,features = markers,stack = F,
             group.by = "cluster_name_slim",cols = pal[5:8],
             pt.size = 0,ncol = 3);p

ggsave("Figure6H_Macrophage_marker_VlnPlot.pdf",p,width = 9,height = 6)

#cell proportion

cell.prop <- table(total_macro$cluster_name_slim,total_macro$disease) %>% as.data.frame.matrix()
#cell.prop <- t(t(cell.prop)/colSums(cell.prop)) %>% as.data.frame()
cell.prop$cell <- row.names(cell.prop)
cell.prop <- pivot_longer(cell.prop,cols = c("Ctl","CAS"),names_to = "disease",values_to = "num")
cell.prop = ddply(cell.prop, 'cell', transform, percent=num/sum(num)*100)

cell.prop$cell <- factor(cell.prop$cell,levels = rev(levels(total_macro$cluster_name_slim)))
cell.prop$disease <- factor(cell.prop$disease,levels = c("CAS","Ctl"))

index <- table(total_macro$disease)
gap <- 100*index["Ctl"]/sum(index)

p <- ggplot(cell.prop,aes(x=percent,y=cell,fill=disease))+
  geom_bar(stat = 'identity',width = 0.8,colour='black')+
  geom_vline(aes(xintercept = gap), colour="black", linetype="dashed",linewidth = 1)+
  scale_fill_manual(values = c("#f26d5b","#2b90d9"))+
  theme_classic()+
  theme(axis.text = element_text(size = 12,colour = "black"),
        axis.title.x = element_text(size = 14,colour = "black"))+
  ylab("")+xlab("Cell ratio (%)");p
ggsave("Figure6I_Macrophage_cell_proportion_bar.pdf",p,width = 6,height = 3)
