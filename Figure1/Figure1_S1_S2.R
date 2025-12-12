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

#######################Pre-treatment##########################################

#samples in list1 were prepared using the 10X library construction method
list1 <- paste0("PVAT",1:6)
#samples in list2 were prepared using the BGI's library construction method
list2 <- c(paste0("CAS",5:7),"CBT3")

data.list <- list()
i=1
for(x in list1){
  raw <- Read10X(data.dir=x)
  data <- CreateSeuratObject(raw,project=x,min.cells=10)
  data <- RenameCells(object=data,add.cell.id=x)
  data$sample <- x
  data <- PercentageFeatureSet(data, pattern = "^MT-", col.name = "percent.mt")
  data <- subset(data, subset = nFeature_RNA > 400 & nFeature_RNA < 8000 & percent.mt < 15)
  data <- SCTransform(data)
  data.list[[i]] <- data
  names(data.list)[i] <- x
  i=i+1
}

for(x in list2){
  raw1 <- Read10X(data.dir=paste0(x,"_1"),gene.column=1)
  raw2 <- Read10X(data.dir=paste0(x,"_2"),gene.column=1)
  count1 <- raw1 %>% as.data.frame()
  count2 <- raw2 %>% as.data.frame()
  count <- merge(count1,count2,by="row.names",all=T)
  count <- column_to_rownames(count,var = "Row.names")
  index <- strsplit(colnames(count),"\\.") %>% lapply(.,function(x){x[1]}) %>% unlist()
  count <- t(count) %>% as.data.frame()
  count <- aggregate(count,by = list(index),FUN = sum)
  count <- column_to_rownames(count,var = "Group.1")
  count <- t(count) %>% as.data.frame()
  data <- CreateSeuratObject(counts = as.sparse(count),project=x,min.cells=10)
  data <- RenameCells(object=data,add.cell.id=x)
  data$sample <- x
  data <- PercentageFeatureSet(data, pattern = "^MT-", col.name = "percent.mt")
  data <- subset(data, subset = nFeature_RNA > 400 & nFeature_RNA < 8000 & percent.mt < 15)
  data <- SCTransform(data)
  data.list[[i]] <- data
  names(data.list)[i] <- x
  i=i+1
}

saveRDS(data.list,file="../data.list.rds")


x <- "CBT4"
count3 <- Read10X(data.dir=paste0(x,"_3"),gene.column=1)
count3 <- CreateSeuratObject(counts = as.sparse(count3),project=x)
count4 <- Read10X(data.dir=paste0(x,"_4"),gene.column=1)
count4 <- CreateSeuratObject(counts = as.sparse(count4),project=x)

count <- merge(count3,count4)
count <- count@assays$RNA$counts

index <- paste(strsplit(colnames(count),"_") %>% lapply(.,function(x){x[1]}) %>% unlist(),
               strsplit(colnames(count),"_") %>% lapply(.,function(x){x[2]}) %>% unlist(),
               sep = "_")
count <- as.data.frame(count) %>% t() %>% as.data.frame()
count <- as.sparse(count)
count <- aggregate(count,by = list(index),FUN = sum)
count <- column_to_rownames(count,var = "Group.1")
count <- t(count)
data <- CreateSeuratObject(counts = as.sparse(count),project=x,min.cells=10)
data <- RenameCells(object=data,add.cell.id=x)
data$sample <- x
data <- PercentageFeatureSet(data, pattern = "^MT-", col.name = "percent.mt")
data <- subset(data, subset = nFeature_RNA > 400 & nFeature_RNA < 8000 & percent.mt < 15)
data <- SCTransform(data)

saveRDS(data,file = "CBT4.rds")

data.list <- readRDS("../data.list.rds")
data.list[[11]] <- data
names(data.list)[11] <- "CBT4"
data.list[[10]] <- NULL
saveRDS(data.list,file = "../data.list.2.rds")

options(future.globals.maxSize= 100000000000000)

data.features <- SelectIntegrationFeatures(object.list = data.list)
data.list <- PrepSCTIntegration(object.list = data.list,assay = "SCT", anchor.features = data.features,verbose = FALSE)
data.anchors <- FindIntegrationAnchors(object.list = data.list, normalization.method = "SCT", anchor.features = data.features, verbose = FALSE)
total <- IntegrateData(anchorset = data.anchors, normalization.method = "SCT", verbose = FALSE)

saveRDS(total,file="../total.2.rds")

##########################Clusering and annotation##############################################

total <- readRDS("total.2.rds")

total <- Reduct(total)

DimPlot(total, reduction = "umap", label=T,group.by = "integrated_snn_res.1")

DefaultAssay(total) <- "SCT"
Idents(total) <- "integrated_snn_res.1"
markers <- c("PTPRC","MS4A1","CD19","BANK1","JCHAIN","SDC1","IGKC","CD3D","CD4","CCR7","SELL","CD8A","GZMK","CCL5","KLRD1","GNLY","NKG7","MKI67","STMN1","TUBA1B","IL1B","S100A8","S100A9","CD68","C1QA","LYVE1","CD1E","LAMP3","FSCN1","PTGDS","IL3RA","LILRA4","F3","CXCL14","CD55","HBB","HBA1","VWF")
DotPlot(total, features = markers, dot.scale = 6)  + theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) + RotatedAxis()

total <- PrepSCTFindMarkers(total)
clu17.markers <- FindMarkers(total,ident.1 = "17",min.pct = 0.5, min.diff.pct = 0.2, logfc.threshold =0.25, only.pos = TRUE, test.use = "wilcox")

Cluster_Highlight_Plot(total, cluster_name = "27", highlight_color = "navy",
                       background_color = "lightgray")

total$cluster_name_slim <- "CD4T cell"
total$cluster_name_slim[total$integrated_snn_res.0.5 %in% c(3,4,6)] <- "B cell"
total$cluster_name_slim[total$integrated_snn_res.0.5 %in% c(16)] <- "Plasma cell"
total$cluster_name_slim[total$integrated_snn_res.0.5 %in% c(5)] <- "CD8T cell"
total$cluster_name_slim[total$integrated_snn_res.0.5 %in% c(11)] <- "NK"
total$cluster_name_slim[total$integrated_snn_res.0.5 %in% c(12)] <- "Proliferating Lymphocyte"
total$cluster_name_slim[total$integrated_snn_res.1 %in% c(22,23)] <- "Monocyte/Neutrophil"
total$cluster_name_slim[total$integrated_snn_res.1 %in% c(30)] <- "Macrophage"
total$cluster_name_slim[total$integrated_snn_res.1 %in% c(33)] <- "DC"
total$cluster_name_slim[total$integrated_snn_res.0.5 %in% c(17)] <- "pDC"
total$cluster_name_slim[total$integrated_snn_res.0.5 %in% c(7)] <- "ADSC"
total$cluster_name_slim[total$integrated_snn_res.0.5 %in% c(20)] <- "Doublets"
total$cluster_name_slim[total$integrated_snn_res.0.5 %in% c(23)] <- "Endothelial"
total <- subset(total,subset = cluster_name_slim != "Doublets")
total <- subset(total,subset = cluster_name_slim != "Endothelial")
total$cluster_name_slim <- factor(total$cluster_name_slim,
                                  levels = c("B cell","Plasma cell","CD4T cell","CD8T cell","NK","Proliferating Lymphocyte","Monocyte/Neutrophil","Macrophage","DC","pDC","ADSC"))
Idents(total) <- "cluster_name_slim"
DimPlot(total)

##########################Remove doublets#######################

library(DoubletFinder)

setwd("~/Project/5.PUMCH_Blood/doublets")

data.list <- SplitObject(total,split.by = "sample")
name <- names(data.list)

i=1
for(data in data.list){
  DefaultAssay(data) <- "RNA"
  data <- SCTransform(data)
  data <- RunPCA(data)
  data <- RunUMAP(data,dims = 1:30)
  
  sweep.res.list <- paramSweep_v3(data, PCs = 1:30, sct = T)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  mpK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  
  annotations <- data@meta.data$cluster_name_slim
  homotypic.prop <- modelHomotypic(annotations)
  DoubletRate = ncol(data)*8*1e-6
  nExp_poi <- round(DoubletRate*length(data$cluster_name_slim))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  data <- doubletFinder_v3(data, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE,sct = T)
  data <- doubletFinder_v3(data, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE,sct = T)
  
  data@meta.data[,"DF_hi.lo"] <- data@meta.data[,ncol(data@meta.data)-2]
  data@meta.data$DF_hi.lo[which(data@meta.data$DF_hi.lo == "Doublet" & data@meta.data[,ncol(data@meta.data)-1] == "Singlet")] <- "Doublet-Low Confidience"
  data@meta.data$DF_hi.lo[which(data@meta.data$DF_hi.lo == "Doublet")] <- "Doublet-High Confidience"
  
  pdf(paste0(names(data.list)[i],".pdf"),width = 8,height = 6)
  p <- DimPlot(data, reduction = "umap", group.by ="DF_hi.lo",cols =c("black","red","gold"))
  print(p)
  dev.off()
  
  data <- subset(data,subset = DF_hi.lo!="Doublet-High Confidience")
  saveRDS(data,file = paste0(names(data.list)[i],".rds"))
  i = i+1
}


data.list <- list()
i=1
for(x in name){
  data.list[i] <- readRDS(paste0(x,".rds"))
  i=i+1
}

names(data.list) <- name

##########################Re-merge data#################################

total <- merge(data.list[[1]],y = data.list[-1])

meta <- total@meta.data
meta <- meta[,1:13]
total@meta.data <- meta

DefaultAssay(total) <- "integrated"
VariableFeatures(total) <- row.names(total)

total <- Reduct(total)
DimPlot(total, reduction = "umap", label=T,group.by = "cluster_name_slim")

Idents(total) <- "integrated_snn_res.1.5"
DimPlot(total, reduction = "umap", label=T,group.by = "integrated_snn_res.1")
DefaultAssay(total) <- "SCT"
DotPlot(total, features = markers, dot.scale = 6)  + theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) + RotatedAxis()

total$sample[total$sample %in% c("PVAT4")] <- "Ctl1"
total$sample[total$sample %in% c("PVAT6")] <- "Ctl2"
total$sample[total$sample %in% c("CBT4")] <- "Ctl3"
total$sample[total$sample %in% c("PVAT1")] <- "CAS1"
total$sample[total$sample %in% c("PVAT2")] <- "CAS2"
total$sample[total$sample %in% c("PVAT3")] <- "CAS3"
total$sample[total$sample %in% c("PVAT5")] <- "CAS4"
total$sample <- factor(total$sample,
                       levels = c("Ctl1","Ctl2","Ctl3","CAS1","CAS2","CAS3","CAS4","CAS5","CAS6","CAS7"))

total$disease <- "CAS"
total$disease[total$sample %in% c("Ctl1","Ctl2","Ctl3")] <- "Ctl"
total$disease <- factor(total$disease,levels = c("Ctl","CAS"))

total$condition <- as.character(total$disease)
total$condition[total$sample %in% c("CAS1","CAS4")] <- "CAS_NS"
total$condition[total$sample %in% c("CAS2","CAS3","CAS5","CAS6","CAS7")] <- "CAS_Smo"
total$condition <- factor(total$condition,levels = c("Ctl","CAS_NS","CAS_Smo"))

#######################Figure S1A QC########################
pdf(file = "FigureS1A_QualityControl.pdf",width = 15,height = 5)
VlnPlot(total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,cols = pal[1:12],pt.size = 0)
dev.off()

######################Figure S1B Pseudo-bulk and PCA###########################
library(SingleCellExperiment)
library(scuttle)
library(GenomicFeatures)
library(rtracklayer)

pesudo_counts <- total@assays$RNA@counts
pesudo_metadata <- total@meta.data


sce_bulk <- SingleCellExperiment(assay = list(counts = pesudo_counts),
                                 colData = pesudo_metadata)

sce_QC <- perCellQCMetrics(sce_bulk)

sce_bulk$is_outlier <- isOutlier(
  metric = sce_QC$total,
  nmads = 3, type = "both", log = TRUE)

sce_bulk <- sce_bulk[,sce_bulk$is_outlier==FALSE]
sce_bulk <- sce_bulk[rowSums(counts(sce_bulk) > 1) >= 10, ]

group <- colData(sce_bulk)[,c("sample")]

pb <- aggregate(t(counts(sce_bulk)),by = list(group), FUN = "sum")
pb <- column_to_rownames(pb,var = "Group.1")

pb_t <- t(pb)

gtf_file <- "genes.gtf"
gtf <- readGFF(gtf_file)
gene <- gtf[gtf$type=="gene",]
gene$length <- gene$end-gene$start+1
gene <- gene[gene$gene_name %in% row.names(pb_t),]
gene <- gene[!duplicated(gene$gene_name),]
row.names(gene) <- gene$gene_name

pb_t <- as.data.frame(pb_t)
pb_t <- pb_t[row.names(pb_t) %in% gene$gene_name,]

gene <- gene[row.names(pb_t),]

kb <- gene$length/1000
rpk <- pb_t/kb
fpkm <- t(t(rpk)*10^6/colSums(pb_t))
fpkm <- as.data.frame(fpkm)
saveRDS(fpkm,file = "fpkm.rds")

#batch effect
library(limma)

#fpkm <- fpkm[,c(9:11,8,1:7)]

group <- factor(c(rep("Ctl",3),rep("CAS",7)),levels = c("Ctl","CAS"))
design <- model.matrix(~0+factor(group))
colnames(design)=levels(factor(group))
rownames(design)=colnames(fpkm)

batch <- c("Batch1","Batch1","Batch3","Batch1","Batch1","Batch1","Batch1","Batch2","Batch2","Batch2")

fpkm <- removeBatchEffect(log2(fpkm+1),batch=batch,design = design)

#PCA
library(ropls)

pca1 <- prcomp(t(fpkm),center = TRUE,scale. = TRUE)
df1 <- pca1$x %>% as.data.frame()
summ1 <- summary(pca1)
xlab1 <- paste0("PC1(",round(summ1$importance[2,1]*100,2),"%)")
ylab1 <- paste0("PC2(",round(summ1$importance[2,2]*100,2),"%)")

df1$Group <- factor(c(rep("Ctl",3),rep("CAS",7)),levels = c("Ctl","CAS"))
df1$sample <- row.names(df1)

p <- ggplot(df1, aes(PC1,PC2,color=Group))+
  geom_point(shape =16,size = 3)+
  ggrepel::geom_text_repel(aes(PC1,PC2,label=sample),size=3)+
  stat_ellipse(aes(group=Group,color = Group))+
  scale_color_manual(values = c("#2b90d9","#f26d5b"))+
  xlab(xlab1)+
  ylab(ylab1)+
  theme_bw();p

ggsave("FigureS1B_PCA.pdf",width = 4,height = 3)

######################Figure 1 and related Figure S1, S2##############################

###Figure 1B_C_D dimplot
Idents(total) <- "integrated_snn_res.0.5"

pdf(file = "FigureS1C_dimplot_unclassified.pdf",width = 8,height = 6)
DimPlot_scCustom(total,colors_use = DiscretePalette_scCustomize(num_colors = 26, palette = "ditto_seq"),figure_plot = T,label = T,label.size = 5,pt.size = 0.5)
dev.off()

total$cluster_name_slim <- factor(total$cluster_name_slim,
                                  levels = c("B cell","Plasma cell","CD4T cell","CD8T cell","NK","Proliferating Lymphocyte","Monocyte/Neutrophil","Macrophage","DC","pDC","ADSC"))
Idents(total) <- "cluster_name_slim"

pal = DiscretePalette_scCustomize(num_colors = 13, palette = "ditto_seq")[-c(8)]

pdf(file = "Figure1B_dimplot_classified.pdf",width = 8,height = 6)
DimPlot_scCustom(total,colors_use = pal[1:12],figure_plot = T,label = T,label.size = 5,pt.size = 0.5)
dev.off()

Idents(total) <- "disease"
pdf(file = "Figure1C_dimplot_disease.pdf",width = 8,height = 6)
DimPlot_scCustom(total,reduction = "umap",label = F,pt.size = 0.1,colors_use = c("#2b90d9","#f26d5b"),figure_plot = T)
dev.off()

Idents(total) <- "sample"
pdf(file = "Figure1D_dimplot_sample.pdf",width = 8,height = 6)
DimPlot_scCustom(total,reduction = "umap",label = F,pt.size = 0.1,colors_use = brewer.pal(10,"Set3"),figure_plot = T)
dev.off()

###Figure 1E dotplot and markers

markers <- c("PTPRC","MS4A1","CD19","BANK1","JCHAIN","SDC1","IGKC","CD3D","CD4","CCR7","SELL","CD8A","GZMK","CCL5","KLRD1","GNLY","NKG7","MKI67","STMN1","TUBA1B","IL1B","S100A8","S100A9","CD68","C1QA","LYVE1","CD1E","LAMP3","FSCN1","PTGDS","IL3RA","LILRA4","F3","CXCL14","CD55")
Idents(total) <- "cluster_name_slim"
DefaultAssay(total) <- "SCT"
pdf(file = "Figure1E_DotPlot.pdf",width = 12,height = 5)
DotPlot(total, features = markers, dot.scale = 6)  + theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) + RotatedAxis()
dev.off()

Idents(total) <- "integrated_snn_res.0.5"
DefaultAssay(total) <- "SCT"
pdf(file = "FigureS1D_DotPlot.pdf",width = 12,height = 6)
DotPlot(total, features = markers, dot.scale = 6)  + theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) + RotatedAxis()
dev.off()

DefaultAssay(total) <- "SCT"
markers <- FindAllMarkers(total,logfc.threshold = 0.25,test.use = "wilcox",min.pct = 0.5,min.diff.pct = 0.2,only.pos = T,recorrect_umi = FALSE)
top30 <- markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
write.csv(markers,file = "total_markers.csv")

#Density plot

library(viridisLite)
library(viridis)
library(ggpointdensity)

Idents(total) <- "disease"
total_sampled <- subset(total,downsample=30000)

Idents(total_sampled) <- "cluster_name_slim"
coord = Embeddings(object = total_sampled, reduction = "umap")
coord = coord[,c(1,2)]
colnames(coord) = c("UMAP_1", "UMAP_2")
coord = data.frame(ID = rownames(coord), coord)
coord$disease <- total_sampled$disease[row.names(coord)]

p2 <- ggplot(data = coord, mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_pointdensity(size = 0.5) +
  stat_density_2d(color=alpha("grey",0.8),bins = 10,adjust = 1.5)+
  #scale_color_viridis(option="D", alpha = 1, begin = 0, end = 1) +
  scale_color_gradientn(colors = rev(brewer.pal(9,"RdYlBu")),oob = scales::squish)+
  theme_classic()+
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())+
  facet_wrap(~disease, ncol = 2);p2
ggsave("FigureS1E_DensityPlot.pdf", p2, width = 10, height = 5)

Idents(total) <- "cluster_name_slim"
pdf(file = "FigureS1E_dimplot_disease_split.pdf",width = 11,height = 5)
DimPlot(total,reduction = "umap",label = T,pt.size = 0.1,split.by = "disease",cols = pal[1:12],raster=F)
dev.off()

###Figure 1F Cell Proportion Barplot
cell_num <- table(total$cluster_name_slim,total$sample) %>% as.data.frame.matrix()
write.csv(cell_num,file = "cell_num.csv")

cell_num <- table(total$cluster_name_slim,total$disease) %>% as.data.frame.table()
colnames(cell_num)<-c("cluster","group","proportion")
new <- cell_num %>% group_by(group) %>% # 按group分组
  mutate(perc = proportion / sum(proportion)) # 添加新列perc，表示value所占比例
new$cluster <- factor(new$cluster,levels = levels(total$cluster_name_slim))
p <- ggplot(new,aes(proportion,group,fill=cluster))+
  geom_bar(stat="identity",position="fill",width = 0.7)+
  scale_fill_manual(values = pal[1:12])+
  theme_classic()+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size=16,colour = "black"),
        axis.text = element_text(size=16,colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size=16))+
  xlab("Cell proportion")
ggsave("Figure1F_Cell_proportion.pdf",width = 12,height = 4)

cell_num <- table(total$cluster_name_slim,total$sample) %>% as.data.frame.table()
colnames(cell_num)<-c("cluster","group","proportion")
new <- cell_num %>% group_by(group) %>% # 按group分组
  mutate(perc = proportion / sum(proportion)) # 添加新列perc，表示value所占比例
new$cluster <- factor(new$cluster,levels = levels(total$cluster_name_slim))
p <- ggplot(new,aes(group,proportion,fill=cluster))+
  geom_bar(stat="identity",position="fill",width = 0.8)+
  scale_fill_manual(values = pal[1:12])+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=16,colour = "black"),
        axis.text = element_text(size=14,colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size=16))+
  ylab("Cell proportion");p
ggsave("FigureS2A_Cell_proportion_sample.pdf",width = 10,height = 5)

saveRDS(total,file = "total_final.2.rds")

### Figure S2 cell proportion box plot
cell.prop <- as.data.frame.matrix((table(total$cluster_name_slim, total$sample)))
cell.prop <- t(cell.prop)/colSums(cell.prop)
cell.prop <- as.data.frame(cell.prop)
cell.prop$Patient <- row.names(cell.prop)
cell.prop <- pivot_longer(cell.prop,cols = setdiff(colnames(cell.prop),"Patient"),names_to = "Cell",values_to = "Freq")

meta <- total@meta.data[,c("sample",'disease')] %>% distinct()
row.names(meta) <- meta$sample
cell.prop$disease <- meta[cell.prop$Patient,"disease"]

cell.prop$Cell[cell.prop$Cell=="Monocyte/Neutrophil"] <- "Monocyte"
index <- c("ADSC","CD4T cell","CD8T cell","B cell","Plasma cell","NK", "Monocyte","Macrophage")
cell.prop <- cell.prop[cell.prop$Cell %in% index,]

cell.prop$Cell <- factor(cell.prop$Cell, levels = index)

my_comparisons <- list( c("Ctl", "CAS") )
p <- ggplot(cell.prop,aes(x=disease,y=Freq))+
  geom_boxplot(aes(fill=disease),outlier.shape = NA)+
  geom_jitter(aes(color=disease),size=1.5)+
  scale_fill_manual(values = alpha(c("#2b90d9","#f26d5b"),0.2))+
  scale_color_manual(values = c("#2b90d9","#f26d5b"))+
  #stat_compare_means(comparisons = my_comparisons,label = "p.format",method = "t.test")+
  scale_y_continuous(labels = scales::percent)+
  theme_pubr()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 16),
        strip.text.x = element_text(size = 12),
        strip.background = element_rect(fill = "white"),
        panel.border = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.line.y = element_line(colour = "black"))+
  ylab("Cell Proportion")+xlab("")+
  facet_wrap(~Cell,ncol = 4,strip.position="bottom",scales = "free_y");p

ggsave("FigureS2B_Cell proportion between Ctl and CAS.pdf",p,width = 6,height = 4)

#######################Ro/e#################################
library("plyr")

do.tissueDist <- function(cellInfo.tb = cellInfo.tb,
                          meta.cluster = cellInfo.tb$meta.cluster,
                          colname.patient = "patient",
                          loc = cellInfo.tb$loc,
                          out.prefix,
                          pdf.width=3,
                          pdf.height=5,
                          verbose=0){
  ##input data 
  library(data.table)
  dir.create(dirname(out.prefix),F,T)
  
  cellInfo.tb = data.table(cellInfo.tb)
  cellInfo.tb$meta.cluster = as.character(meta.cluster)
  
  if(is.factor(loc)){
    cellInfo.tb$loc = loc
  }else{cellInfo.tb$loc = as.factor(loc)}
  
  loc.avai.vec <- levels(cellInfo.tb[["loc"]])
  count.dist <- unclass(cellInfo.tb[,table(meta.cluster,loc)])[,loc.avai.vec]
  freq.dist <- sweep(count.dist,1,rowSums(count.dist),"/")
  freq.dist.bin <- floor(freq.dist * 100 / 10)
  print(freq.dist.bin)
  
  {
    count.dist.melt.ext.tb <- test.dist.table(count.dist)
    p.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="p.value")
    OR.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="OR")
    OR.dist.mtx <- as.matrix(OR.dist.tb[,-1])
    rownames(OR.dist.mtx) <- OR.dist.tb[[1]]
  }
  
  sscVis::plotMatrix.simple(OR.dist.mtx,
                            out.prefix=sprintf("%s.OR.dist",out.prefix),
                            show.number=F,
                            waterfall.row=T,par.warterfall = list(score.alpha = 2,do.norm=T),
                            exp.name=expression(italic(OR)),
                            z.hi=4,
                            palatte=viridis::viridis(7),
                            pdf.width = 4, pdf.height = pdf.height)
  if(verbose==1){
    return(list("count.dist.melt.ext.tb"=count.dist.melt.ext.tb,
                "p.dist.tb"=p.dist.tb,
                "OR.dist.tb"=OR.dist.tb,
                "OR.dist.mtx"=OR.dist.mtx))
  }else{
    return(OR.dist.mtx)
  }
}

test.dist.table <- function(count.dist,min.rowSum=0)
{
  count.dist <- count.dist[rowSums(count.dist)>=min.rowSum,,drop=F]
  sum.col <- colSums(count.dist)
  sum.row <- rowSums(count.dist)
  count.dist.tb <- as.data.frame(count.dist)
  setDT(count.dist.tb,keep.rownames=T)
  count.dist.melt.tb <- melt(count.dist.tb,id.vars="rn")
  colnames(count.dist.melt.tb) <- c("rid","cid","count")
  count.dist.melt.ext.tb <- as.data.table(ldply(seq_len(nrow(count.dist.melt.tb)), function(i){
    this.row <- count.dist.melt.tb$rid[i]
    this.col <- count.dist.melt.tb$cid[i]
    this.c <- count.dist.melt.tb$count[i]
    other.col.c <- sum.col[this.col]-this.c
    this.m <- matrix(c(this.c,
                       sum.row[this.row]-this.c,
                       other.col.c,
                       sum(sum.col)-sum.row[this.row]-other.col.c),
                     ncol=2)
    res.test <- fisher.test(this.m)
    data.frame(rid=this.row,
               cid=this.col,
               p.value=res.test$p.value,
               OR=res.test$estimate)
  }))
  count.dist.melt.ext.tb <- merge(count.dist.melt.tb,count.dist.melt.ext.tb,
                                  by=c("rid","cid"))
  count.dist.melt.ext.tb[,adj.p.value:=p.adjust(p.value,"BH")]
  return(count.dist.melt.ext.tb)
}


meta.tb <- total@meta.data
meta.tb$disease <- as.character(meta.tb$disease) %>% factor(.,levels=c("Ctl","CAS"))
OR.list <- do.tissueDist(cellInfo.tb=meta.tb,
                         out.prefix="./Roe",
                         meta.cluster=meta.tb$cluster_name_slim,
                         loc=meta.tb$disease,
                         pdf.width=4,pdf.height=6,verbose=1)

matrix <- OR.list[["OR.dist.tb"]]
matrix <- column_to_rownames(matrix,var = "rid")
matrix <- matrix[levels(total$cluster_name_slim),]
matrix <- t(matrix)
write.csv(matrix,file = "Figure1F_Roe_res.csv")

label <- matrix
label[label>1] <- "+++"
label[label<1&label>0.8] <- "++"
label[label<0.8&label>0.2] <- "+"
label[label<0.2&label>0] <- "+/-"