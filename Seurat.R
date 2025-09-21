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

pal = DiscretePalette_scCustomize(num_colors = 12, palette = "ditto_seq")[-c(8)]

#######################Pre-treatment##########################################
list1 <- paste0("PVAT",1:6)
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
#count1 <- Read10X(data.dir=paste0(x,"_1"),gene.column=1) 
#count1 <- CreateSeuratObject(counts = as.sparse(count1),project=x)
#count2 <- Read10X(data.dir=paste0(x,"_2"),gene.column=1)
#count2 <- CreateSeuratObject(counts = as.sparse(count2),project=x)
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

##########################Clusering##############################################
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
#total$cluster_name_slim[total$integrated_snn_res.0.5 %in% c(20)] <- "Doublets"
#total$cluster_name_slim[total$integrated_snn_res.0.5 %in% c(23)] <- "Endothelial"
#total <- subset(total,subset = cluster_name_slim != "Doublets")
total$cluster_name_slim <- factor(total$cluster_name_slim,
                                       levels = c("B cell","Plasma cell","CD4T cell","CD8T cell","NK","Proliferating Lymphocyte","Monocyte/Neutrophil","Macrophage","DC","pDC","ADSC"))
Idents(total) <- "cluster_name_slim"
DimPlot(total)

#total_clu22 <- subset(total,subset = integrated_snn_res.1.5 == "22")
#total_clu22 <- Reduct(total_clu22)
#DimPlot(total_clu22, reduction = "umap", label=T,group.by = "integrated_snn_res.1.5")
#Idents(total_clu22) <- "integrated_snn_res.1.5"
#DefaultAssay(total_clu22) <- "SCT"
#DotPlot(total_clu22, features = markers, dot.scale = 6)  + theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) + RotatedAxis()

#total$cluster_name_slim <- as.character(total$cluster_name_slim)
#index <- colnames(total_clu22)[total_clu22$integrated_snn_res.1.5 %in% c(5)]
#total$cluster_name_slim[index] <- "Doublets"

#total_clu25 <- subset(total,subset = integrated_snn_res.1.5 %in% c("25"))
#total_clu25 <- Reduct(total_clu25)
#DimPlot(total_clu25, reduction = "umap", label=T,group.by = "integrated_snn_res.1.5")
#Idents(total_clu25) <- "integrated_snn_res.1.5"
#DefaultAssay(total_clu25) <- "SCT"
#DotPlot(total_clu25, features = markers, dot.scale = 6)  + theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) + RotatedAxis()

#index <- colnames(total_clu27)[total_clu27$integrated_snn_res.1.5 %in% c(8)]
#total$cluster_name_slim[index] <- "Doublets"

#total$cluster_name_slim[total$integrated_snn_res.1.5 %in% c(39)] <- "B cell"

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

######################Project to public data####################################
total_public <- readRDS("./PublicData/total_BAT_WAT.rds")
total <- readRDS("total_final.2.rds")

total_public <- RunUMAP(total_public,reduction = "harmony", dims = 1:30,return.model = T)
total_public <- SCTransform(total_public)
DefaultAssay(total) <- "SCT"

DefaultAssay(total) <- "RNA"
total <- NormalizeData(total)
total <- FindVariableFeatures(total)
total <- ScaleData(total)

anchors <- FindTransferAnchors(
  reference = total_public,
  query = total,
  features = row.names(total@assays$integrated$data),
  normalization.method = "LogNormalize",
  reference.reduction = "pca",
  reference.assay = "RNA",
  query.assay = "RNA",
  dims = 1:50
)

total <- MapQuery(
  anchorset = anchors,
  query = total,
  reference = total_public,
  refdata = list(
    cluster_name_slim = "cluster_name_slim",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "pca", 
  reduction.model = "umap"
)

total$predicted.cluster_name_slim <- factor(total$predicted.cluster_name_slim,
                                                levels = c("ADSC","B","T","Myeloid","Endothelial"))

DimPlot(total, reduction = "ref.umap", group.by = "predicted.cluster_name_slim", label = F, pt.size = 0.8, label.size = 3, repel = TRUE)
DimPlot(total, reduction = "ref.umap", group.by = "cluster_name_slim", label = F, pt.size = 0.8, label.size = 3, repel = TRUE)


######################Figure 1##############################
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


###Figure 1B_C_D dimplot
Idents(total) <- "integrated_snn_res.0.5"

pdf(file = "Figure1B_dimplot_unclassified.pdf",width = 8,height = 6)
DimPlot_scCustom(total,colors_use = DiscretePalette_scCustomize(num_colors = 26, palette = "ditto_seq"),figure_plot = T,label = T,label.size = 5,pt.size = 0.5)
dev.off()

total$cluster_name_slim <- factor(total$cluster_name_slim,
                                  levels = c("B cell","Plasma cell","CD4T cell","CD8T cell","NK","Proliferating Lymphocyte","Monocyte/Neutrophil","Macrophage","DC","pDC","ADSC"))
Idents(total) <- "cluster_name_slim"

pal = DiscretePalette_scCustomize(num_colors = 13, palette = "ditto_seq")[-c(8)]

pdf(file = "Figure1C_dimplot_classified.pdf",width = 8,height = 6)
DimPlot_scCustom(total,colors_use = pal[1:12],figure_plot = T,label = T,label.size = 5,pt.size = 0.5)
dev.off()

Idents(total) <- "disease"
pdf(file = "Figure1D_dimplot_disease.pdf",width = 8,height = 6)
DimPlot_scCustom(total,reduction = "umap",label = F,pt.size = 0.1,colors_use = c("#2b90d9","#f26d5b"),figure_plot = T)
dev.off()

Idents(total) <- "sample"
pdf(file = "FigureS1_dimplot_sample.pdf",width = 8,height = 6)
DimPlot_scCustom(total,reduction = "umap",label = F,pt.size = 0.1,colors_use = brewer.pal(10,"Set3"),figure_plot = T)
dev.off()

###Figure 1E dotplot and markers

markers <- c("PTPRC","MS4A1","CD19","BANK1","JCHAIN","SDC1","IGKC","CD3D","CD4","CCR7","SELL","CD8A","GZMK","CCL5","KLRD1","GNLY","NKG7","MKI67","STMN1","TUBA1B","IL1B","S100A8","S100A9","CD68","C1QA","LYVE1","CD1E","LAMP3","FSCN1","PTGDS","IL3RA","LILRA4","F3","CXCL14","CD55")
Idents(total) <- "cluster_name_slim"
DefaultAssay(total) <- "SCT"
pdf(file = "Figure1E_DotPlot.pdf",width = 12,height = 5)
DotPlot(total, features = markers, dot.scale = 6)  + theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) + RotatedAxis()
dev.off()

DefaultAssay(total) <- "SCT"
markers <- FindAllMarkers(total,logfc.threshold = 0.25,test.use = "wilcox",min.pct = 0.5,min.diff.pct = 0.2,only.pos = T,recorrect_umi = FALSE)
top30 <- markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
write.csv(markers,file = "total_markers.csv")
write.csv(top30,file = "TableS_markers_top30.csv")

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
ggsave("Figure1F_Cell_proportion_sample.pdf",width = 10,height = 5)

###Figure S1 QC
pdf(file = "SupplementaryFigure1_QualityControl.pdf",width = 15,height = 5)
VlnPlot(total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,cols = pal[1:12],pt.size = 0)
dev.off()

saveRDS(total,file = "total_final.2.rds")

#cell proportion box plot
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

ggsave("Figure1F_Cell proportion between Ctl and CAS.pdf",p,width = 6,height = 4)



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

######################Pseudo-bulk###########################
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

ggsave("Figure1S_PCA.pdf",width = 4,height = 3)

######################Figure 2 ADSC##############################
total <- readRDS("total_final.2.rds")

total_ADSC <- subset(total,subset = cluster_name_slim == "ADSC")

total_ADSC <- Reduct(total_ADSC)

DimPlot(total_ADSC, reduction = "umap", label=T,group.by = "integrated_snn_res.1.5")

markers <- c("PTPRC","MS4A1","CD19","BANK1","JCHAIN","SDC1","IGKC","CD3D","CD4","CCR7","SELL","CD8A","GZMK","CCL5","KLRD1","GNLY","NKG7","MKI67","STMN1","TUBA1B","IL1B","S100A8","S100A9","CD68","C1QA","LYVE1","CD1E","LAMP3","FSCN1","PTGDS","IL3RA","LILRA4","F3","CXCL14","CD55")
Idents(total_ADSC) <- "integrated_snn_res.1.5"
DefaultAssay(total_ADSC) <- "SCT"
DotPlot(total_ADSC, features = markers, dot.scale = 6)  + theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) + RotatedAxis()


total_ADSC <- subset(total_ADSC,subset = integrated_snn_res.1.5 %in% c(0:7,9,10,12,13))
total_ADSC <- Reduct(total_ADSC)

DimPlot(total_ADSC, reduction = "umap", label=T,group.by = "integrated_snn_res.1.5")

#re-integrate
data.list <- SplitObject(total_ADSC,split.by = "sample")
data.list[["CAS1"]] <- NULL
data.list[["Ctl1"]] <- NULL

data.list <- lapply(X=data.list,FUN=function(x){
  DefaultAssay(x) <- "RNA"
  x[['integrated']] <- NULL
  x <- SCTransform(x)
})


options(future.globals.maxSize= 100000000000000)

data.features <- SelectIntegrationFeatures(object.list = data.list)
data.list <- PrepSCTIntegration(object.list = data.list,assay = "SCT", anchor.features = data.features,verbose = FALSE)
data.anchors <- FindIntegrationAnchors(object.list = data.list, normalization.method = "SCT", anchor.features = data.features, verbose = FALSE,k.filter = 50)
total_ADSC <- IntegrateData(anchorset = data.anchors, normalization.method = "SCT", verbose = FALSE, k.weight = 50)

total_ADSC <- Reduct(total_ADSC)

DimPlot(total_ADSC, reduction = "umap", label=T,group.by = "integrated_snn_res.1.5")

markers <- c("PTPRC","MS4A1","CD19","BANK1","JCHAIN","SDC1","IGKC","CD3D","CD4","CCR7","SELL","CD8A","GZMK","CCL5","KLRD1","GNLY","NKG7","MKI67","STMN1","TUBA1B","IL1B","S100A8","S100A9","CD68","C1QA","LYVE1","CD1E","LAMP3","FSCN1","PTGDS","IL3RA","LILRA4","F3","CXCL14","CD55")
Idents(total_ADSC) <- "integrated_snn_res.1.5"
DefaultAssay(total_ADSC) <- "SCT"
DotPlot(total_ADSC, features = markers, dot.scale = 6)  + theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) + RotatedAxis()

total_ADSC <- subset(total_ADSC,subset = integrated_snn_res.1.5 %in% c(0:8,10,12,16))
total_ADSC <- Reduct(total_ADSC)

DimPlot(total_ADSC, reduction = "umap", label=T,group.by = "integrated_snn_res.1.5")

total_ADSC$cluster_name_slim <- "ADSC_CXCL14"
total_ADSC$cluster_name_slim[total_ADSC$integrated_snn_res.1.5 %in% c(1,5,6,8,10)] <- "ADSC_CD55"
total_ADSC$cluster_name_slim <- factor(total_ADSC$cluster_name_slim,
                                       levels = c("ADSC_CXCL14","ADSC_CD55"))

DimPlot(total_ADSC, reduction = "umap", label=T,group.by = "cluster_name_slim")

Idents(total_ADSC) <- "cluster_name_slim"

pdf(file = "Figure2A_dimplot_ADSC.pdf",width = 6,height = 5)
DimPlot(total_ADSC,reduction = "umap",label = F,pt.size = 0.15,cols = pal[c(1,3)])
dev.off()


markers <- FindAllMarkers(total_ADSC,logfc.threshold = 0.25,test.use = "wilcox",min.pct = 0.5,min.diff.pct = 0.2,only.pos = T,recorrect_umi = FALSE)
write.csv(markers,file="total_markers_ADSC.csv")

top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf(file = "Figure2B_Top10geneHeatmap.pdf",width = 8,height = 6)
DoHeatmap(total_ADSC, features = top10$gene,slot = "data",group.colors = pal[c(1,3)]) + 
  scale_fill_gradientn(colors = c("grey98","gold","orange","red","red3")) + 
  theme(axis.text=element_text(size=12))+
  NoLegend()
dev.off()

pdf(file = "Figure2C_featureplot.pdf",width = 10,height = 5)
Plot_Density_Custom(seurat_object = total_ADSC, features = c("CXCL14","CD55"))
dev.off()

p <- FeaturePlot(total_ADSC,features = c("CD34","ITGB1","F3"),order = T,cols = c("lightgrey","firebrick"),ncol = 3);p
ggsave("FigureS2B_feature_ADSC.pdf",p,width = 12,height = 4)

p <- FeaturePlot(total_ADSC,features = c("PDGFRA","PDGFRB"),order = T,cols = c("lightgrey","firebrick"),ncol = 2);p
ggsave("Feature_ADSC_PDGFR.pdf",p,width = 8,height = 4)

p <- FeaturePlot(total_ADSC,features = c("GPC3"),order = T,cols = c("lightgrey","firebrick"),ncol = 1);p
ggsave("Feature_ADSC_GPC3.pdf",p,width = 4,height = 4)

saveRDS(total_ADSC,file = "total_ADSC.rds")

marker <- list(gene1=c("PI16","ANXA3","FN1","SEMA3C","SMPD3","PTGS2","DPP4","COL14A1","IGFBP4","CD55"),
               gene2=c("COL4A1","COL4A2","CXCL14","SPARCL1","COL5A3","FABP4","COL15A1","LPL","COL6A3","IGFBP7"))

DefaultAssay(total_ADSC)

expr <- as.matrix(total_ADSC@assays$SCT@data)
expr <- expr[rowSums(expr)>0,]
gsea_mat <- gsva(expr, marker, method="ssgsea", verbose=T, parallel.sz=0, ssgsea.norm=T)

gsea_mat <- t(gsea_mat)

total_ADSC <- AddMetaData(total_ADSC,gsea_mat)

p <- VlnPlot(total_ADSC,features = c("gene1","gene2"),cols = pal[c(1,3)],pt.size = 0)

ggsave("Signature_ADSC.pdf",p,width = 6,height = 4)

#cell proportion box plot
cell.prop <- as.data.frame.matrix((table(total_ADSC$cluster_name_slim, total_ADSC$sample)))
cell.prop <- t(cell.prop)/colSums(cell.prop)
cell.prop <- as.data.frame(cell.prop)
cell.prop$Patient <- row.names(cell.prop)
cell.prop <- pivot_longer(cell.prop,cols = setdiff(colnames(cell.prop),"Patient"),names_to = "Cell",values_to = "Freq")

meta <- total_ADSC@meta.data[,c("sample",'disease')] %>% distinct()
row.names(meta) <- meta$sample
cell.prop$disease <- factor(meta[cell.prop$Patient,"disease"],
                            levels = c("Ctl","CAS"))
cell.prop$Cell <- factor(cell.prop$Cell,
                         levels = levels(total_ADSC$cluster_name_slim))

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
  facet_wrap(~Cell,ncol = 2,strip.position="bottom",scales = "free_y");p

ggsave("Figure2_Cell proportion between Ctl and CAS_ADSC.pdf",p,width = 5,height = 3)

#Roe
meta.tb <- total_ADSC@meta.data
meta.tb$disease <- as.character(meta.tb$disease) %>% factor(.,levels=c("Ctl","CAS"))
OR.list <- do.tissueDist(cellInfo.tb=meta.tb,
                         out.prefix="./Roe",
                         meta.cluster=meta.tb$cluster_name_slim,
                         loc=meta.tb$disease,
                         pdf.width=4,pdf.height=6,verbose=1)

matrix <- OR.list[["OR.dist.tb"]]
matrix <- column_to_rownames(matrix,var = "rid")
matrix <- matrix[levels(total_ADSC$cluster_name_slim),]
matrix <- t(matrix)
write.csv(matrix,file = "Figure2_Roe_res_ADSC.csv")

label <- matrix
label[label>1] <- "+++"
label[label<1&label>0.8] <- "++"
label[label<0.8&label>0.2] <- "+"
label[label<0.2&label>0] <- "+/-"


#secreted angiogenesis

library(msigdbr)

m_df_gobp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
m_df_gobp <- m_df_gobp %>% split(x = .$gene_symbol, f = .$gs_name)

m_df_gocc <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "CC")
m_df_gocc <- m_df_gocc %>% split(x = .$gene_symbol, f = .$gs_name)

angio <- m_df_gobp$GOBP_SPROUTING_ANGIOGENESIS
vascular <- m_df_gobp$GOBP_VASCULATURE_DEVELOPMENT
vascular <- m_df_gobp$GOBP_BLOOD_VESSEL_REMODELING

secreted <- read.table("SecretedProteins.txt",header = F)

gene <- intersect(vascular,secreted$V1)


total_ADSC <- readRDS("total_ADSC.rds")
Idents(total_ADSC)
#test1
DefaultAssay(total_ADSC)
total_ADSC <- PrepSCTFindMarkers(total_ADSC)
markers <- FindMarkers(total_ADSC,ident.1 = "ADSC_CD55",logfc.threshold = 0.2,test.use = "wilcox",only.pos = T)

#test2
total_ADSC_CD55 <- subset(total_ADSC,subset = cluster_name_slim == "ADSC_CD55")
total_ADSC_CD55$disease <- factor(total_ADSC_CD55$disease,levels = c("Ctl","CAS"))
Idents(total_ADSC_CD55) <- "disease"
DefaultAssay(total_ADSC_CD55)
total_ADSC_CD55 <- PrepSCTFindMarkers(total_ADSC_CD55)
markers <- FindMarkers(total_ADSC_CD55,ident.1 = "CAS",logfc.threshold = 0.2,test.use = "wilcox",only.pos = T)


res <- markers[gene,] %>% na.omit() %>% row.names()


VlnPlot(total_ADSC,features = gene,pt.size = 0,ncol = 4)&
  stat_compare_means(method = "wilcox.test", 
                      label = "p.signif",##星号设置
                      label.x.npc = "middle",#位置
                     label.y.npc = 0.9,hide.ns = T,
                      size = 5)&
  theme(axis.title.x = element_blank())
ggsave("VlnPlot_remodeling_secreted.pdf",width = 12,height = 9)

#CD55-CXCL14

DefaultAssay(total_ADSC) <- "RNA"
total_ADSC <- NormalizeData(total_ADSC)

exp <- FetchData(total_ADSC,vars = c("CD55","CXCL14"))

p <- ggplot(exp,aes(x=CD55,y=CXCL14))+
  geom_point(color=alpha("black",0.7))+
  theme_bw()+
  theme(axis.text = element_text(color = "black",size = 14),
        axis.title = element_text(size = 16));p

ggsave("Point_CD55_CXCL14.pdf",p,width = 4,height = 4)

######################Figure 3 ADSC's function##############################
library(WebGestaltR)
library(Hmisc)

markers <- read.csv("total_markers_ADSC.csv",row.names = 1)

top20 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

#volcano plot

res <- markers
res$avg_log2FC[res$cluster=="ADSC_CD55"] <- -res$avg_log2FC[res$cluster=="ADSC_CD55"]

res$p_val_adj[res$p_val_adj==0] <- 1*10e-300

limit=1
data_vol <- res %>% subset(.,select=c(2,5))
colnames(data_vol) <- c("log2FoldChange","pvalue")
data_vol$sig[((data_vol$log2FoldChange<limit)&(data_vol$log2FoldChange>-limit))|(data_vol$pvalue>0.05)] <- "no"
data_vol$sig[(data_vol$log2FoldChange>=limit)&(data_vol$pvalue<=0.05)] <- "up"
data_vol$sig[(data_vol$log2FoldChange<=-limit)&(data_vol$pvalue<=0.05)] <- "down"
data_vol$label <- row.names(data_vol)
data_vol$label[data_vol$sig=="no"] <- ""

#saveRDS(data_vol,file = "data_vol.rds")

p <- ggplot(data_vol,aes(log2FoldChange,-1*log10(pvalue),
                         color = sig))+geom_point()+
  geom_text_repel(aes(log2FoldChange,-1*log10(pvalue),label=label),size=3)+
  labs(x="log2(FoldChange)",y="-log10(P-value)")
p <- p + scale_color_manual(values =c("#0072B5","grey","#BC3C28"))+
  geom_hline(yintercept=-log10(0.05),linetype=4)+
  geom_vline(xintercept=c(-limit,limit),linetype=4)
p <- p + theme(panel.grid =element_blank())+
  theme(axis.line = element_line(size=0))+theme_classic()
p <- p  +guides(colour = FALSE)
p <- p +theme(axis.text=element_text(size=20,color = "black"),axis.title=element_text(size=20))
p

ggsave("FigureS2A_VolcanoPlot.pdf",p,width = 8,height = 6)



gene <- markers[markers$cluster=="ADSC_CXCL14","gene"]
gene <- markers[markers$cluster=="ADSC_CD55","gene"]

GOBP <- WebGestaltR(enrichMethod = "ORA", organism = "hsapiens",
                    enrichDatabase = "geneontology_Biological_Process_noRedundant",
                    enrichDatabaseType = "genesymbol",
                    interestGene = gene, interestGeneType = "genesymbol",
                    referenceSet = "genome",
                    fdrThr = 0.2,
                    isOutput = F)
GOBP <- GOBP[GOBP$FDR<0.05&GOBP$overlap>3,]

gobp_plot <- GOBP
gobp_plot$FDR <- (-log10(gobp_plot$pValue))
gobp_plot <- gobp_plot[order(gobp_plot$enrichmentRatio,decreasing = T),]
gobp_plot <- gobp_plot[1:10,] %>% na.omit()
gobp_plot$description <- capitalize(gobp_plot$description)
gobp_plot$description <- factor(gobp_plot$description,levels = rev(gobp_plot$description))

p <- ggplot(gobp_plot,aes(x=enrichmentRatio,y=description))+
  geom_bar(stat="identity",fill=alpha(pal[3],0.5))+
  geom_text(aes(x=0.1,y=description,label=description),hjust=0)+
  #scale_fill_gradientn(colours = rev(brewer.pal(9,"Reds")))+
  theme_classic()+
  scale_x_continuous(expand = c(0,0))+
  theme(legend.position = "right",
        legend.text = element_text(color = "black",size = 14),
        legend.title = element_blank(),
        axis.text = element_text(color = "black",size = 14),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(color = "black",size = 16),
        axis.title.y = element_blank())+
  xlab("Enrichment Ratio");p

ggsave("Figure3A_GeneEnrichment_ADSC_CD55.pdf",p,width = 10,height = 4)

#GSEA
library("GSVA")
library("msigdbr")

m_df_gobp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
m_df_gobp <- m_df_gobp %>% split(x = .$gene_symbol, f = .$gs_name)

features <- c("GOBP_RESPONSE_TO_TUMOR_NECROSIS_FACTOR",
              "GOBP_ENDOTHELIAL_CELL_PROLIFERATION",
              "GOBP_ARTERIAL_ENDOTHELIAL_CELL_DIFFERENTIATION",
              "GOBP_VASCULOGENESIS",
              "GOBP_ARTERY_MORPHOGENESIS",
              "GOBP_SPROUTING_ANGIOGENESIS",
              "GOBP_POSITIVE_REGULATION_OF_SPROUTING_ANGIOGENESIS",
              "GOBP_POSITIVE_REGULATION_OF_ADIPOSE_TISSUE_DEVELOPMENT",
              "GOBP_EXTRACELLULAR_MATRIX_ASSEMBLY",
              "GOBP_FAT_CELL_DIFFERENTIATION",
              "GOBP_REGULATION_OF_MESENCHYMAL_STEM_CELL_DIFFERENTIATION",
              "GOBP_MESENCHYMAL_CELL_PROLIFERATION",
              "GOBP_MESENCHYMAL_STEM_CELL_PROLIFERATION",
              "GOBP_POSITIVE_REGULATION_OF_MESENCHYMAL_STEM_CELL_PROLIFERATION",
              "GOBP_FIBROBLAST_PROLIFERATION",
              "GOBP_CELL_CHEMOTAXIS",
              "GOBP_POSITIVE_CHEMOTAXIS")

m_df_gobp_2 <- m_df_gobp[features]

expr <- as.matrix(total_ADSC@assays$SCT@data)
expr <- expr[rowSums(expr)>0,]
gsea_mat <- gsva(expr, m_df_gobp_2, method="ssgsea", verbose=T, parallel.sz=0, ssgsea.norm=T)
total_ADSC[["geneList"]] <- CreateAssayObject(gsea_mat)
DefaultAssay(total_ADSC) <- "geneList"
diff_gobp <- FindMarkers(total_ADSC,ident.1 = "ADSC_CXCL14",ident.2 = "ADSC_CD55",logfc.threshold = 0,min.pct = 0,min.diff.pct = 0,test.use = "wilcox")


features <-  c("GOBP-CELL-CHEMOTAXIS","GOBP-SPROUTING-ANGIOGENESIS")
data <- FetchData(total_ADSC,vars = features)
data$group <- total_ADSC$cluster_name_slim[row.names(data)]
colnames(data)[1:2] <- c("BP1","BP2")

p <- ggviolin(data, x = "group", # 分组列
              y = "BP2", # 基因列
              fill = "group", #按分组填充颜色
              alpha = 1,#透明图 0-1
              width = 0.5, #宽度
              title = features[2],
              legend = "top",legend.title = "",#legend及位置
              font.legend = c(12, "plain", "black"),
              ylab="Enrichment score", xlab=FALSE, #xy轴标签，去掉X轴标签
              font.y = 16,#xy轴标题大小
              x.text.angle = 45, y.text.angle = 0,#xy轴标题角度
              font.tickslab = c(15,"plain","black"), #xy轴刻度大小/样式/颜色
              add = "boxplot", #添加图 "dotplot", "jitter", "boxplot", "point"
              add.params = list(fill = "white", #填充颜色 白色
                                width = 0.1,#宽度
                                linetype = 1)#线型
)+stat_compare_means(method = "wilcox.test", 
                     #label = "p.signif",##星号设置
                     aes(label = "p.format"), #显示方式
                     label.x.npc ="middle",#位置
                     size = 5)+
  scale_fill_manual(values = pal[c(1,3)]);p
ggsave(paste0("Figure3B_",features[2],".pdf"),width = 4,height = 5)


#highly-expressed genes
gene1 <- m_df_gobp$GOBP_CELL_CHEMOTAXIS
gene2 <- read.table("angiogenesis.txt")$V1

gene1 <- intersect(gene1,markers[markers$cluster=="ADSC_CXCL14","gene"])
gene2 <- intersect(gene2,markers[markers$cluster=="ADSC_CD55","gene"])
gene <- c(gene1,gene2)

gene <- c("CXCL14","CXCL12","VCAM1","MYOC",
          "CD248","SCARA5","FN1","ECM1")

DefaultAssay(total_ADSC) <- "SCT"
DotPlot(total_ADSC,features = gene)+
  RotatedAxis()+
  scale_color_gradientn(colors = rev(brewer_pal(palette = "RdYlBu")(10)))+
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 14))
ggsave("Figure3C_Genes_dotplot.pdf",width = 7,height = 4)

######################Figure 3 monocle##############################
library(monocle)

exp.matrix <- as.sparse(total_ADSC@assays$RNA@counts)
feature_ann <- data.frame(gene_id=rownames(exp.matrix),gene_short_name=rownames(exp.matrix))
rownames(feature_ann) <- rownames(total_ADSC@assays$RNA@counts)
fd <- new("AnnotatedDataFrame", data = feature_ann)
sample_ann <- total_ADSC@meta.data
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
cds <- orderCells(cds,reverse = T)

pdf("Figure3J_Pseudotime_cell_name.pdf",width = 5,height = 5)
plot_cell_trajectory(cds, color_by = "cluster_name_slim",show_branch_points = F,cell_size = 1)+
  scale_color_manual(values = pal[c(1,3)])
dev.off()

pdf("Figure3J_Pseudotime_time.pdf",width = 5,height = 5)
plot_cell_trajectory(cds, color_by = "Pseudotime",show_branch_points = F,cell_size = 1)
dev.off()

pdf("Figure3J_Pseudotime_disease.pdf",width = 5,height = 5)
plot_cell_trajectory(cds, color_by = "disease",show_branch_points = F,cell_size = 1)+
  scale_color_manual(values = c("#f26d5b","#2b90d9"))
dev.off()

saveRDS(cds,file = "monocle_ADSC.rds")

pData(cds)$KLF3 = log2(exprs(cds)['KLF3',]+1)
pdf("Figure3J_Pseudotime_KLF3.pdf",width = 5,height = 5)
plot_cell_trajectory(cds, color_by = "KLF3",show_branch_points = F,cell_size = 1) + scale_color_gsea()
dev.off()

pData(cds)$JDP2 = log2(exprs(cds)['JDP2',]+1)
pdf("Figure3J_Pseudotime_JDP2.pdf",width = 5,height = 5)
plot_cell_trajectory(cds, color_by = "JDP2",show_branch_points = F,cell_size = 1) + scale_color_gsea()
dev.off()

pData(cds)$JDP2 = log2(exprs(cds)['JDP2',]+1)
pdf("Figure3J_Pseudotime_JDP2.pdf",width = 5,height = 5)
plot_cell_trajectory(cds, color_by = "JDP2",show_branch_points = F,cell_size = 1) + scale_color_gsea()
dev.off()

pData(cds)$NR2F2 = log2(exprs(cds)['NR2F2',]+1)
pdf("Figure3J_Pseudotime_NR2F2.pdf",width = 5,height = 5)
plot_cell_trajectory(cds, color_by = "NR2F2",show_branch_points = F,cell_size = 1) + scale_color_gsea()
dev.off()

for(x in c("KLF3","FOSB","LMX1A","JDP2","NR2F2","TCF4","CEBPA","MEOX2","EPAS1")){
  pData(cds)$TF = log2(exprs(cds)[x,]+1)
  pdf(paste0("Figure3J_Pseudotime_",x,".pdf"),width = 5,height = 5)
  p <- plot_cell_trajectory(cds, color_by = "TF",show_branch_points = F,cell_size = 1) + scale_color_gsea()
  print(p)
  dev.off()
}

######################Figure 3 cytoTrace##############################
library(CytoTRACE)
library(reticulate)
library(RColorBrewer)

matrix <- as.matrix(total_ADSC@assays$RNA@counts)
phe <- as.character(total_ADSC$cluster_name_slim)
names(phe) <- row.names(total_ADSC@meta.data)

results <- CytoTRACE(matrix)
plotCytoTRACE(results, phenotype = phe,colors = pal[c(1,3)])

data <- data.frame(order = results$CytoTRACE,
                   cluster = total_ADSC$cluster_name_slim[names(results$CytoTRACE)],
                   row.names = names(results$CytoTRACE))
test <- t.test(data[data$cluster=="ADSC_CXCL14","order"],
               data[data$cluster=="ADSC_CD55","order"],
               paired = F)

######################Figure 3 SCENIC##############################
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(tidyr)
library(data.table)

scenic_matrix <- t(as.matrix(total_ADSC@assays$RNA@counts))
scenic_matrix <- scenic_matrix[,!str_detect(colnames(scenic_matrix),"\\.")]
scenic_matrix <- as.data.frame(scenic_matrix)
fwrite(scenic_matrix,file = "./scenic_count.csv",row.names = T,col.names = T)

sce_SCENIC <- open_loom("./sce_SCENIC.loom")
exprMat <- get_dgem(sce_SCENIC)
exprMat_log <- log2(exprMat+1)
regulons_incidMat <- get_regulons(sce_SCENIC, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(sce_SCENIC, column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(sce_SCENIC)

cellTypes <- data.frame(type=total_ADSC$cluster_name_slim,
                        index=row.names(total_ADSC@meta.data),
                        row.names=row.names(total_ADSC@meta.data))

rss <- calcRSS(AUC=getAUC(regulonAUC),
               cellAnnotation=cellTypes[colnames(regulonAUC),"type"])

rss <- as.data.frame(rss)
rss$label <- row.names(rss)

rss_sub <- rss[order(rss$ADSC_CXCL14,decreasing = T),]
rss_sub$order <- as.numeric(1:nrow(rss_sub))
rss_sub$label[!rss_sub$order %in% c(1:5)] <- ""
rss_sub$color[rss_sub$order %in% c(1:5)] <- "show"
rss_sub$color[!rss_sub$order %in% c(1:5)] <- "no"

p <- ggplot(rss_sub,aes(x=order,y=ADSC_CXCL14,color=color))+
  geom_point(size=1)+
  geom_text_repel(aes(order,ADSC_CXCL14,label=label),max.overlaps = 50)+
  scale_color_manual(values=c("blue","red"))+
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",fill = "transparent"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(family = "sans",colour = "black",size=14),
        axis.title = element_text(family = "sans",colour = "black",size=14),
        legend.position = "none")+
  xlab("Regulon")+
  ylab("Regulon specificity score (RSS)")

ggsave("Figure3I_Regulon_ADSC_CXCL14.pdf",p,width = 3,height = 4)

AUC <- getAUC(regulonAUC)
AUC <- AUC[,colnames(total_ADSC)]

total_ADSC[["Regulon"]] <- CreateAssayObject(AUC)
DefaultAssay(total_ADSC) <- "Regulon"

VariableFeatures(total_ADSC) <- row.names(total_ADSC)
total_ADSC <- ScaleData(total_ADSC,features = row.names(total_ADSC))

p <- DoHeatmap(total_ADSC,features = c("FOSB(-)","LMX1A(+)","KLF3(+)","JDP2(+)","LMX1A(-)",
                                  "NR2F2(+)","TCF4(+)","CEBPA(+)","MEOX2(+)","EPAS1(+)"),
          slot = "scale.data",group.colors = pal[c(1,3)])+
  scale_fill_gradientn(colors = brewer.pal(9,"YlGn"))+
  theme(axis.text=element_text(size=16,color = "black"))
ggsave("FigureS2E_Heatmap_Regulon.pdf",p,width = 8,height = 6)

######################Figure 4 innate immune cell##############################

total <- readRDS("total_final.2.rds")

total_immune <- subset(total,idents=c("Macrophage","Monocyte/Neutrophil","NK","DC","pDC"))

data.list <- SplitObject(total_immune,split.by = "sample")

data.list <- lapply(X=data.list,FUN=function(x){
  DefaultAssay(x) <- "RNA"
  x[['integrated']] <- NULL
  x <- SCTransform(x)
})


options(future.globals.maxSize= 100000000000000)

data.features <- SelectIntegrationFeatures(object.list = data.list)
data.list <- PrepSCTIntegration(object.list = data.list,assay = "SCT", anchor.features = data.features,verbose = FALSE)
data.anchors <- FindIntegrationAnchors(object.list = data.list, normalization.method = "SCT", anchor.features = data.features, verbose = FALSE)
total_immune <- IntegrateData(anchorset = data.anchors, normalization.method = "SCT", verbose = FALSE)

total_immune <- Reduct(total_immune)

DimPlot(total_immune, reduction = "umap", label=T,group.by = "integrated_snn_res.1.5")

markers <- c("PTPRC","MS4A1","CD19","BANK1","JCHAIN","SDC1","IGKC","CD3D","CD4","CCR7","SELL","CD8A","GZMK","CCL5","KLRD1","GNLY","NKG7","MKI67","STMN1","TUBA1B","IL1B","S100A8","S100A9","CD68","C1QA","LYVE1","CD1E","LAMP3","FSCN1","PTGDS","IL3RA","LILRA4","F3","CXCL14","CD55")
Idents(total_immune) <- "integrated_snn_res.1.5"
DefaultAssay(total_immune) <- "SCT"
DotPlot(total_immune, features = markers, dot.scale = 6)  + theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) + RotatedAxis()

markers <- FindMarkers(total_immune,ident.1 = c(4,10),logfc.threshold = 0.25,test.use = "wilcox",min.pct = 0.5,min.diff.pct = 0.2,only.pos = T,recorrect_umi = FALSE)

total_immune <- subset(total_immune,subset = integrated_snn_res.1.5 %in% c(0:6,8,10:14,16:18,21))
total_immune <- Reduct(total_immune)

DimPlot(total_immune, reduction = "umap", label=T,group.by = "integrated_snn_res.1.5")

markers <- c("PTPRC","NCAM1","FCGR3A","ITGA1","B3GAT1","KLRD1","GNLY","NKG7","IL1B","S100A8","S100A9","CSF3R","CD68","C1QA","LYVE1","XCR1","CLEC9A","CD1C","CD1E","LAMP3","FSCN1","PTGDS","IL3RA","LILRA4")
Idents(total_immune) <- "integrated_snn_res.1.5"
DefaultAssay(total_immune) <- "SCT"
DotPlot(total_immune, features = markers, dot.scale = 6)  + theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) + RotatedAxis()

list <- list(bright = c("NCAM1","CCR7","CSF2","CXCR3","IFNG","IL2RB","IL7R","KIT","KLRC1","KLRD1","NCR1","SELL"),
             dim = c("FCGR3A","NCAM1","CX3CR1","CXCR1","ITGB2","KIR3DL1","KLRC2","KLRG1","PRF1"),
             resident = c("CCR7","CD69","EOMES","ICAM1","IL7R","ITGA1","KIR3DL1","NCR1","NCR2","NCR3","SELL","TBX21"),
             adaptive = c("B3GAT1","KLRC2","LILRB1"))

total_immune <- AddModuleScore(total_immune,features = list,assay = "SCT")

VlnPlot(total_immune,features = c("Cluster1","Cluster2","Cluster3","Cluster4"),idents = c(1,3,4,8,10))

total_immune$cluster_name_slim <- "NK_CD56bright"
total_immune$cluster_name_slim[total_immune$integrated_snn_res.1.5 == 1] <- "NK_CD56dim"
total_immune$cluster_name_slim[total_immune$integrated_snn_res.1.5 == 8] <- "NK_tr"
total_immune$cluster_name_slim[total_immune$integrated_snn_res.1.5 %in% c(0,12)] <- "cDC2"
total_immune$cluster_name_slim[total_immune$integrated_snn_res.1.5 %in% c(13)] <- "cDC1"
total_immune$cluster_name_slim[total_immune$integrated_snn_res.1.5 %in% c(15)] <- "DC_LAMP3+"
total_immune$cluster_name_slim[total_immune$integrated_snn_res.1.5 %in% c(2,5,11)] <- "Neutrophil"
total_immune$cluster_name_slim[total_immune$integrated_snn_res.1.5 %in% c(7,14)] <- "Macrophage"
total_immune$cluster_name_slim[total_immune$integrated_snn_res.1.5 %in% c(16)] <- "Monocyte"
total_immune$cluster_name_slim[total_immune$integrated_snn_res.1.5 %in% c(6,9)] <- "pDC"
total_immune$cluster_name_slim <- factor(total_immune$cluster_name_slim,
                                         levels = c("NK_CD56bright","NK_CD56dim","NK_tr","cDC1","cDC2","DC_LAMP3+","Monocyte","Neutrophil","Macrophage","pDC"))
Idents(total_immune) <- "cluster_name_slim"
DimPlot(total_immune, reduction = "umap", label=T,group.by = "cluster_name_slim")

saveRDS(total_immune,file = "total_innate_immune.rds")

pdf(file = "Figure4A_dimplot_innate_immune.pdf",width = 8,height = 6)
DimPlot_scCustom(total_immune,colors_use = pal[1:10],figure_plot = T,label = T,label.size = 5,pt.size = 0.5)
dev.off()

Idents(total_immune) <- "disease"
pdf(file = "Figure4B_dimplot_disease.pdf",width = 7,height = 6)
DimPlot_scCustom(total_immune,reduction = "umap",label = F,pt.size = 0.1,colors_use = c("#2b90d9","#f26d5b"),figure_plot = T)
dev.off()

markers <- FindAllMarkers(total_immune,logfc.threshold = 0.25,test.use = "wilcox",min.pct = 0.5,min.diff.pct = 0.2,only.pos = T,recorrect_umi = FALSE)
write.csv(markers,file = "total_markers_innate_immune.csv")

markers <- c("NCAM1",'FCGR3A',"CD69","NKG7","XCR1","CLEC9A","CD1C","CLEC10A","LAMP3","CCL19","CDKN1C",'LILRB2',"FCN1","S100A8","S100A9","CSF3R","CD163","C1QA","MRC1","IL3RA","LILRA4","TCF4")
Idents(total_immune) <- "cluster_name_slim"
DefaultAssay(total_immune) <- "SCT"
pdf(file = "FigureS3A_DotPlot_innate_immune.pdf",width = 12,height = 5)
DotPlot(total_immune, features = markers, dot.scale = 6)  + theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) + RotatedAxis()
dev.off()

#cell proportion in bar plot

cell.prop <- table(total_immune$cluster_name_slim,total_immune$disease) %>% as.data.frame.matrix()
#cell.prop <- t(t(cell.prop)/colSums(cell.prop)) %>% as.data.frame()
cell.prop$cell <- row.names(cell.prop)
cell.prop <- pivot_longer(cell.prop,cols = c("Ctl","CAS"),names_to = "disease",values_to = "num")
cell.prop = ddply(cell.prop, 'cell', transform, percent=num/sum(num)*100)

cell.prop$cell <- factor(cell.prop$cell,levels = rev(c("NK_CD56bright","NK_CD56dim","NK_tr","cDC1","cDC2","DC_LAMP3+","Monocyte","Neutrophil","Macrophage","pDC")))
cell.prop$disease <- factor(cell.prop$disease,levels = c("CAS","Ctl"))

p <- ggplot(cell.prop,aes(x=percent,y=cell,fill=disease))+
  geom_bar(stat = 'identity',width = 0.8,colour='black')+
  geom_vline(aes(xintercept = 28.84), colour="black", linetype="dashed",linewidth = 1)+
  scale_fill_manual(values = c("#f26d5b","#2b90d9"))+
  theme_classic()+
  theme(axis.text = element_text(size = 12,colour = "black"),
        axis.title.x = element_text(size = 14,colour = "black"))+
  ylab("")+xlab("Cell ratio (%)");p
ggsave("Figure4C_cell_proportion_bar.pdf",p,width = 5,height = 3.5)



#cell proportion box plot
cell.prop <- as.data.frame.matrix((table(total_immune$cluster_name_slim, total_immune$sample)))
cell.prop <- t(cell.prop)/colSums(cell.prop)
cell.prop <- as.data.frame(cell.prop)
cell.prop$Patient <- row.names(cell.prop)
cell.prop <- pivot_longer(cell.prop,cols = setdiff(colnames(cell.prop),"Patient"),names_to = "Cell",values_to = "Freq")

meta <- total_immune@meta.data[,c("sample",'disease')] %>% distinct()
row.names(meta) <- meta$sample
cell.prop$disease <- factor(meta[cell.prop$Patient,"disease"],
                            levels = c("Ctl","CAS"))
cell.prop$Cell <- factor(cell.prop$Cell,
                         levels = levels(total_immune$cluster_name_slim))

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

ggsave("Figure4C_Cell proportion between Ctl and CAS_innate.pdf",p,width = 8,height = 6)

#Roe
meta.tb <- total_immune@meta.data
meta.tb$disease <- as.character(meta.tb$disease) %>% factor(.,levels=c("Ctl","CAS"))
OR.list <- do.tissueDist(cellInfo.tb=meta.tb,
                         out.prefix="./Roe",
                         meta.cluster=meta.tb$cluster_name_slim,
                         loc=meta.tb$disease,
                         pdf.width=4,pdf.height=6,verbose=1)

matrix <- OR.list[["OR.dist.tb"]]
matrix <- column_to_rownames(matrix,var = "rid")
matrix <- matrix[levels(total_immune$cluster_name_slim),]
write.csv(matrix,file = "Figure4_Roe_res_innate.csv")

label <- matrix
label[label>1] <- "+++"
label[label<1&label>0.8] <- "++"
label[label<0.8&label>0.2] <- "+"
label[label<0.2&label>0] <- "+/-"

#Correlation heatmap
library(psych)

total_SCT_mat <- as.matrix(total_immune@assays$SCT@data)
total_SCT_mat_ClusterMeanSplit <- apply(total_SCT_mat,1,function(input){aggregate(input,list(total_immune@meta.data$cluster_name_slim),mean)[,2]})
rownames(total_SCT_mat_ClusterMeanSplit) <- levels(total_immune@meta.data$cluster_name_slim)
total_SCT_mat_ClusterMeanSplit <- t(total_SCT_mat_ClusterMeanSplit)
total_SCT_mat_ClusterMeanSplit <- total_SCT_mat_ClusterMeanSplit[rownames(total_SCT_mat),]

res <- corr.test(total_SCT_mat_ClusterMeanSplit)
r <- res$r
#bk <- seq(0.7,1,by=0.01)
pheatmap::pheatmap(r,cellwidth = 16,cellheight = 16,
                   width = 5,height = 4,
                   fontsize_col = 8,fontsize_row = 8,
                   #color = colorRampPalette(c("steelblue","white","chocolate1"))(length(bk)),
                   #legend_breaks=seq(0.75,1,0.05),breaks=bk,
                   treeheight_col = 20,treeheight_row = 20,
                   filename = "Figure4D_correlation.pdf")

#GSVA
library(msigdbr)
library(GSVA)

total_immune$disease <- factor(total_immune$disease,levels = c("Ctl","CAS"))
total_immune$cell_disease <- paste(total_immune$cluster_name_slim,total_immune$disease,sep = "-")
total_immune$cell_disease <- factor(total_immune$cell_disease,
                                    levels = paste(rep(levels(total_immune$cluster_name_slim),each=2),rep(levels(total_immune$disease),10),sep = "-"))

m_df_hallmark = msigdbr(species = "Homo sapiens", category = "H")
m_df_hallmark <- m_df_hallmark %>% split(x = .$gene_symbol, f = .$gs_name)

expr <- as.matrix(total_immune@assays$SCT@data)
gsea_mat <- gsva(expr, m_df_hallmark, method="ssgsea", verbose=T, parallel.sz=0, ssgsea.norm=T)

total_immune[["hallmark"]] <- CreateAssayObject(gsea_mat)

matrix <- AverageExpression(total_immune,assays = "hallmark",group.by = "cell_disease") %>% as.data.frame()
 
test <- matrix(nrow = nrow(matrix),ncol = ncol(matrix)) %>%
  as.data.frame(row.names = row.names(matrix))
colnames(test) <- colnames(matrix)
for(x in unique(total_immune$cluster_name_slim)){
  print(x)
  DefaultAssay(total_immune) <- "hallmark"
  Idents(total_immune) <- "cell_disease"
  res <- FindMarkers(total_immune,ident.1 = paste(x,"CAS",sep = "-"),ident.2 = paste(x,"Ctl",sep = "-"), min.pct = 0, logfc.threshold =0.2, only.pos = TRUE, test.use = "wilcox")
  hallmark <- res[res$p_val<0.05,] %>% row.names()
  x <- gsub("_",".",x)
  x <- gsub("\\+",".",x)
  test[hallmark,paste("hallmark",x,"CAS",sep = ".")] <- 1
}
test[is.na(test)] <- 0

test <- test[rowSums(test)>0,]
matrix <- matrix[row.names(test),]
test[test==1] <- "*"
test[test==0] <- ""

row.names(matrix) <- gsub("HALLMARK-","",row.names(matrix))
row.names(matrix) <- gsub("-"," ",row.names(matrix))
row.names(matrix) <- str_to_title(row.names(matrix))
row.names(test) <- row.names(matrix)

anno_col <- data.frame(Group = rep(levels(total_immune$disease),10),
                       Cell = rep(levels(total_immune$cluster_name_slim),each=2),
                       row.names = colnames(matrix))
anno_color <- list(Cell = c(NK_CD56bright = pal[1],
                            NK_CD56dim = pal[2],
                            NK_tr = pal[3],
                            cDC1 = pal[4],
                            cDC2 = pal[5],
                            "DC_LAMP3+" = pal[6],
                            Monocyte = pal[7],
                            Neutrophil = pal[8],
                            Macrophage = pal[9],
                            pDC = pal[10]),
                   Group = c(Ctl = "#2b90d9", CAS = "#f26d5b"))

pheatmap::pheatmap(matrix,scale = "row",cluster_cols = F,
                   display_numbers = test,
                   annotation_col = anno_col,
                   annotation_colors = anno_color,
                   show_colnames = F,
                   gaps_col = seq(2,18,2),
                   width = 12,
                   fontsize_row = 14,
                   filename = "Figure4E_hallmark_heatmap.pdf")

#DEGs
Idents(total_immune) <- "cell_disease"

diffExp <- function(total,ident1,ident2){
  DefaultAssay(total) <- "SCT"
  Cell_SCT_diff <- FindMarkers(total, ident.1 = ident1, ident.2 = ident2, only.pos = F, test.use="wilcox",logfc.threshold = 0.2,min.pct = 0.4,min.diff.pct = 0.1,recorrect_umi = FALSE)
  return(Cell_SCT_diff)
}

cluster1 <- levels(total_immune$cluster_name_slim)
DefaultAssay(total_immune) <- "SCT"
data.features <- row.names(total_immune)
RB_genes = c(grep("^RPS", data.features, v=T),grep("^RPL", data.features, v=T))
MT_genes = grep("^MT-", data.features, v=T)
remove <- c(RB_genes,MT_genes)
total_remove <- subset(total_immune,features = setdiff(data.features,remove))


dir.create("DEG_innate")
for(x in cluster1){
  print(x)
  diff_res <- diffExp(total_remove,paste0(x,"-CAS"),paste0(x,"-Ctl"))
  write.csv(diff_res,paste0("./DEG_innate/",x,"_diff_sig.csv"),row.names=T)
}

###Emapplot
library(clusterProfiler)
library(enrichplot)
library(ReactomePA)

GO_database <- 'org.Hs.eg.db'
cluster1 <- levels(total_immune$cluster_name_slim)

res <- list()
res2 <- c()
i=1
for(x in cluster1){
  print(x)
  diff_res <- read.csv(paste0("./DEG_innate/",x,"_diff_sig.csv"),row.names = 1)
  diff_res <- na.omit(diff_res)
  diff_gene <- diff_res[(diff_res$p_val<0.01)&(abs(diff_res$avg_log2FC)>0.5),]
  diff_gene <- diff_gene[!(row.names(diff_gene) %in% remove),]
  diff_gene <- diff_gene[order(diff_gene$avg_log2FC,decreasing = T),]
  diff_gene$change <- "Up"
  diff_gene$change[diff_gene$avg_log2FC<0] <- "Down"
  gene <- bitr(row.names(diff_gene[diff_gene$change=="Up",]),fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  res[[i]] <- gene$ENTREZID
  names(res)[i] <- x
  res2 <- c(res2,row.names(diff_gene))
  i=i+1
}

res2 <- unique(res2)

xx <- compareCluster(res, fun="enrichPathway")
xx2 <- pairwise_termsim(xx)
p <- emapplot(xx2,pie="Count",showCategory = 5,shadowtext=F)+
  scale_fill_manual(values = pal[1:10])+theme(legend.position = "none");p
ggsave("Figure4E_DEG_Reactome_network_down_10.pdf",p,width = 7,height = 6)

#Cytokines
cytokine <- c("IL1A","IL1B","IL1RN","IL18","IL2","IL4","IL7","IL9","IL13","IL15","IL3",
              "IL5","CSF2","IL6","IL11","CSF3","IL12A","IL12B","LIF","OSM","IL10","IL20","TXLNA",
              "IL16","IL17A","IL17B","IFNA1","IFNB1","IFNG","CD40LG","LTB","TNF","TNFSF9","TNFSF13",
              "CD70","TNFSF8","FASLG","TNFSF18","TNFSF14","TNFSF4","TNFSF13B","TNFSF10",
              "TNFSF12","TNFSF11","TGFB1","TGFB2","TGFB3","EPO","TPO","FLT3LG","KITLG",
              "CSF1","MST1")

exp <- FetchData(total_immune,vars = cytokine)
exp <- exp[,colnames(exp)[!str_detect(colnames(exp),pattern = "rna")]]
total_immune[["cytokine"]] <- CreateAssayObject(t(exp))

matrix <- AverageExpression(total_immune,assays = "cytokine",group.by = "cell_disease") %>% as.data.frame()

test <- matrix(nrow = nrow(matrix),ncol = ncol(matrix)) %>%
  as.data.frame(row.names = row.names(matrix))
colnames(test) <- colnames(matrix)
for(x in unique(total_immune$cluster_name_slim)){
  print(x)
  DefaultAssay(total_immune) <- "cytokine"
  Idents(total_immune) <- "cell_disease"
  res <- FindMarkers(total_immune,ident.1 = paste(x,"CAS",sep = "-"),ident.2 = paste(x,"Ctl",sep = "-"), min.pct = 0, logfc.threshold =0, only.pos = TRUE, test.use = "wilcox")
  cytokine <- res[res$p_val<0.05,] %>% row.names()
  x <- gsub("_",".",x)
  x <- gsub("\\+",".",x)
  test[cytokine,paste("cytokine",x,"CAS",sep = ".")] <- 1
}
test[is.na(test)] <- 0

test <- test[rowSums(test)>0,]
matrix <- matrix[row.names(test),]
test[test==1] <- "*"
test[test==0] <- ""

anno_col <- data.frame(Group = rep(levels(total_immune$disease),10),
                       Cell = rep(levels(total_immune$cluster_name_slim),each=2),
                       row.names = colnames(matrix))
anno_color <- list(Cell = c(NK_CD56bright = pal[1],
                            NK_CD56dim = pal[2],
                            NK_tr = pal[3],
                            cDC1 = pal[4],
                            cDC2 = pal[5],
                            "DC_LAMP3+" = pal[6],
                            Monocyte = pal[7],
                            Neutrophil = pal[8],
                            Macrophage = pal[9],
                            pDC = pal[10]),
                   Group = c(Ctl = "#2b90d9", CAS = "#f26d5b"))

pheatmap::pheatmap(matrix,scale = "row",cluster_cols = F,
                   annotation_col = anno_col,
                   annotation_colors = anno_color,
                   display_numbers = test,
                   fontsize_number = 12,
                   show_colnames = F,
                   gaps_col = seq(2,18,2),
                   treeheight_row = 20,
                   width = 12,
                   cellwidth = 16,cellheight = 12,
                   fontsize_row = 14,
                   filename = "Figure4F_cytokine_heatmap.pdf")


DotPlot(total_immune,features = cytokine,group.by = "cell_disease")+RotatedAxis()

######################Figure 4 NK cell##############################

total_immune <- readRDS("total_innate_immune.rds")

total_NK <- subset(total_immune,subset = cluster_name_slim %in% c("NK_CD56bright","NK_CD56dim","NK_tr"))

data.list <- SplitObject(total_NK,split.by = "sample")

data.list <- lapply(X=data.list,FUN=function(x){
  DefaultAssay(x) <- "RNA"
  x[['integrated']] <- NULL
  x <- SCTransform(x)
})

options(future.globals.maxSize= 100000000000000)

data.features <- SelectIntegrationFeatures(object.list = data.list)
data.list <- PrepSCTIntegration(object.list = data.list,assay = "SCT", anchor.features = data.features,verbose = FALSE)
data.anchors <- FindIntegrationAnchors(object.list = data.list, normalization.method = "SCT", anchor.features = data.features, verbose = FALSE,k.filter = 50)
total_NK <- IntegrateData(anchorset = data.anchors, normalization.method = "SCT",k.weight = 50, verbose = FALSE)

total_NK <- Reduct(total_NK)

Idents(total_NK) <- "cluster_name_slim"

DimPlot(total_NK,group.by = "integrated_snn_res.1",label = T)

markers <- c("PTPRC","NCAM1","FCGR3A","ITGA1","B3GAT1","KLRD1","GNLY","NKG7","IL1B","S100A8","S100A9","CSF3R","CD68","C1QA","LYVE1","XCR1","CLEC9A","CD1C","CD1E","LAMP3","FSCN1","PTGDS","IL3RA","LILRA4")
Idents(total_NK) <- "integrated_snn_res.1"
DefaultAssay(total_NK) <- "SCT"
DotPlot(total_NK, features = markers, dot.scale = 6)  + theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) + RotatedAxis()

total_NK$cluster_name_slim <- "NK_CD56bright"
total_NK$cluster_name_slim[total_NK$integrated_snn_res.1 %in% c(0,3)] <- "NK_CD56dim"
total_NK$cluster_name_slim[total_NK$integrated_snn_res.1 %in% c(4)] <- "NK_tr"
total_NK$cluster_name_slim[total_NK$integrated_snn_res.1 %in% c(7)] <- "Remove"

total_NK <- subset(total_NK,subset = cluster_name_slim != "Remove")

total_NK$cluster_name_slim <- factor(total_NK$cluster_name_slim,
                                     levels = c("NK_CD56bright","NK_CD56dim","NK_tr"))

Idents(total_NK) <- "cluster_name_slim"

pdf(file = "Figure5_dimplot_NK.pdf",width = 6,height = 4)
DimPlot_scCustom(total_NK,colors_use = pal[1:3],figure_plot = T,label = T,label.size = 5,pt.size = 0.5)
dev.off()

total_NK$disease <- factor(total_NK$disease,levels = c("Ctl","CAS"))
Idents(total_NK) <- "disease"
pdf(file = "Figure5_dimplot_NK_disease.pdf",width = 5.5,height = 4)
DimPlot_scCustom(total_NK,reduction = "umap",label = F,pt.size = 0.5,colors_use = c("#2b90d9","#f26d5b"),figure_plot = T)
dev.off()

markers <- FindAllMarkers(total_NK,logfc.threshold = 0.25,test.use = "wilcox",min.pct = 0.5,min.diff.pct = 0.2,only.pos = T,recorrect_umi = FALSE)
write.csv(markers,file = "Figure5_total_markers_NK.csv")

#cell proportion in bar plot

cell.prop <- table(total_NK$cluster_name_slim,total_NK$disease) %>% as.data.frame.matrix()
#cell.prop <- t(t(cell.prop)/colSums(cell.prop)) %>% as.data.frame()
cell.prop$cell <- row.names(cell.prop)
cell.prop <- pivot_longer(cell.prop,cols = c("Ctl","CAS"),names_to = "disease",values_to = "num")
cell.prop = ddply(cell.prop, 'cell', transform, percent=num/sum(num)*100)

cell.prop$cell <- factor(cell.prop$cell,levels = rev(c("NK_CD56bright","NK_CD56dim","NK_tr","cDC1","cDC2","DC_LAMP3+","Monocyte","Neutrophil","Macrophage","pDC")))
cell.prop$disease <- factor(cell.prop$disease,levels = c("CAS","Ctl"))

p <- ggplot(cell.prop,aes(x=percent,y=cell,fill=disease))+
  geom_bar(stat = 'identity',width = 0.8,colour='black')+
  geom_vline(aes(xintercept = 26.94), colour="black", linetype="dashed",linewidth = 1)+
  scale_fill_manual(values = c("#f26d5b","#2b90d9"))+
  theme_classic()+
  theme(axis.text = element_text(size = 12,colour = "black"),
        axis.title.x = element_text(size = 14,colour = "black"))+
  ylab("")+xlab("Cell ratio (%)");p
ggsave("Figure5_cell_proportion_bar_NK.pdf",p,width = 5,height = 3.5)



#cell proportion box plot
cell.prop <- as.data.frame.matrix((table(total_NK$cluster_name_slim, total_NK$sample)))
cell.prop <- t(cell.prop)/colSums(cell.prop)
cell.prop <- as.data.frame(cell.prop)
cell.prop$Patient <- row.names(cell.prop)
cell.prop <- pivot_longer(cell.prop,cols = setdiff(colnames(cell.prop),"Patient"),names_to = "Cell",values_to = "Freq")

meta <- total_NK@meta.data[,c("sample",'disease')] %>% distinct()
row.names(meta) <- meta$sample
cell.prop$disease <- factor(meta[cell.prop$Patient,"disease"],
                            levels = c("Ctl","CAS"))
cell.prop$Cell <- factor(cell.prop$Cell,
                         levels = levels(total_NK$cluster_name_slim))

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
  facet_wrap(~Cell,ncol = 3,strip.position="bottom",scales = "free_y");p

ggsave("Figure5_NK_Cell proportion between Ctl and CAS.pdf",p,width = 6,height = 3)

#Roe
meta.tb <- total_NK@meta.data
meta.tb$disease <- as.character(meta.tb$disease) %>% factor(.,levels=c("Ctl","CAS"))
OR.list <- do.tissueDist(cellInfo.tb=meta.tb,
                         out.prefix="./Roe",
                         meta.cluster=meta.tb$cluster_name_slim,
                         loc=meta.tb$disease,
                         pdf.width=4,pdf.height=6,verbose=1)

matrix <- OR.list[["OR.dist.tb"]]
matrix <- column_to_rownames(matrix,var = "rid")
matrix <- matrix[levels(total_NK$cluster_name_slim),]
write.csv(matrix,file = "Figure5_NK_Roe_res.csv")

label <- matrix
label[label>1] <- "+++"
label[label<1&label>0.8] <- "++"
label[label<0.8&label>0.2] <- "+"
label[label<0.2&label>0] <- "+/-"

#Correlation heatmap
library(psych)

total_SCT_mat <- as.matrix(total_NK@assays$SCT@data)
total_SCT_mat_ClusterMeanSplit <- apply(total_SCT_mat,1,function(input){aggregate(input,list(total_NK@meta.data$cluster_name_slim),mean)[,2]})
rownames(total_SCT_mat_ClusterMeanSplit) <- levels(total_NK@meta.data$cluster_name_slim)
total_SCT_mat_ClusterMeanSplit <- t(total_SCT_mat_ClusterMeanSplit)
total_SCT_mat_ClusterMeanSplit <- total_SCT_mat_ClusterMeanSplit[rownames(total_SCT_mat),]

res <- corr.test(total_SCT_mat_ClusterMeanSplit)
r <- res$r
#bk <- seq(0.7,1,by=0.01)
pheatmap::pheatmap(r,cellwidth = 16,cellheight = 16,
                   width = 5,height = 4,
                   fontsize_col = 8,fontsize_row = 8,
                   #color = colorRampPalette(c("steelblue","white","chocolate1"))(length(bk)),
                   #legend_breaks=seq(0.75,1,0.05),breaks=bk,
                   treeheight_col = 20,treeheight_row = 20,
                   filename = "Figure5_NK_correlation.pdf")


#DEGs

library(msigdbr)
library(GSVA)

total_NK$disease <- factor(total_NK$disease,levels = c("Ctl","CAS"))
total_NK$cell_disease <- paste(total_NK$cluster_name_slim,total_NK$disease,sep = "-")
total_NK$cell_disease <- factor(total_NK$cell_disease,
                                    levels = paste(rep(levels(total_NK$cluster_name_slim),each=2),rep(levels(total_NK$disease),3),sep = "-"))


Idents(total_NK) <- "cell_disease"

diffExp <- function(total,ident1,ident2){
  DefaultAssay(total) <- "SCT"
  Cell_SCT_diff <- FindMarkers(total, ident.1 = ident1, ident.2 = ident2, only.pos = F, test.use="wilcox",logfc.threshold = 0.2,min.pct = 0.4,min.diff.pct = 0.1,recorrect_umi = FALSE)
  return(Cell_SCT_diff)
}

cluster1 <- levels(total_NK$cluster_name_slim)
DefaultAssay(total_NK) <- "SCT"
data.features <- row.names(total_NK)
RB_genes = c(grep("^RPS", data.features, v=T),grep("^RPL", data.features, v=T))
MT_genes = grep("^MT-", data.features, v=T)
remove <- c(RB_genes,MT_genes)
total_remove <- subset(total_NK,features = setdiff(data.features,remove))


dir.create("DEG_NK")
for(x in cluster1){
  print(x)
  diff_res <- diffExp(total_remove,paste0(x,"-CAS"),paste0(x,"-Ctl"))
  write.csv(diff_res,paste0("./DEG_NK/",x,"_diff_sig.csv"),row.names=T)
}

###Emapplot
library(clusterProfiler)
library(enrichplot)
library(ReactomePA)

GO_database <- 'org.Hs.eg.db'
cluster1 <- levels(total_NK$cluster_name_slim)

res <- list()
res2 <- c()
i=1
for(x in cluster1){
  print(x)
  diff_res <- read.csv(paste0("./DEG_NK/",x,"_diff_sig.csv"),row.names = 1)
  diff_res <- na.omit(diff_res)
  diff_gene <- diff_res[(diff_res$p_val<0.01)&(abs(diff_res$avg_log2FC)>0.5),]
  diff_gene <- diff_gene[!(row.names(diff_gene) %in% remove),]
  diff_gene <- diff_gene[order(diff_gene$avg_log2FC,decreasing = T),]
  diff_gene$change <- "Up"
  diff_gene$change[diff_gene$avg_log2FC<0] <- "Down"
  gene <- bitr(row.names(diff_gene[diff_gene$change=="Down",]),fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  res[[i]] <- gene$ENTREZID
  names(res)[i] <- x
  res2 <- c(res2,row.names(diff_gene))
  i=i+1
}

res2 <- unique(res2)

xx <- compareCluster(res, fun="enrichPathway")
xx2 <- pairwise_termsim(xx)
p <- emapplot(xx2,pie="Count",showCategory = 10,shadowtext=F)+
  scale_fill_manual(values = pal[1:3])+theme(legend.position = "none");p
ggsave("Figure5_NK_DEG_Reactome_network_down_10.pdf",p,width = 7,height = 6)

#Cytokines
cytokine <- c("IL1A","IL1B","IL1RN","IL18","IL2","IL4","IL7","IL9","IL13","IL15","IL3",
              "IL5","CSF2","IL6","IL11","CSF3","IL12A","IL12B","LIF","OSM","IL10","IL20","TXLNA",
              "IL16","IL17A","IL17B","IFNA1","IFNB1","IFNG","CD40LG","LTB","TNF","TNFSF9","TNFSF13",
              "CD70","TNFSF8","FASLG","TNFSF18","TNFSF14","TNFSF4","TNFSF13B","TNFSF10",
              "TNFSF12","TNFSF11","TGFB1","TGFB2","TGFB3","EPO","TPO","FLT3LG","KITLG",
              "CSF1","MST1")

exp <- FetchData(total_NK,vars = cytokine)
exp <- exp[,colnames(exp)[!str_detect(colnames(exp),pattern = "rna")]]
total_NK[["cytokine"]] <- CreateAssayObject(t(exp))

matrix <- AverageExpression(total_NK,assays = "cytokine",group.by = "cell_disease") %>% as.data.frame()

test <- matrix(nrow = nrow(matrix),ncol = ncol(matrix)) %>%
  as.data.frame(row.names = row.names(matrix))
colnames(test) <- colnames(matrix)
for(x in unique(total_NK$cluster_name_slim)){
  print(x)
  DefaultAssay(total_NK) <- "cytokine"
  Idents(total_NK) <- "cell_disease"
  res <- FindMarkers(total_NK,ident.1 = paste(x,"CAS",sep = "-"),ident.2 = paste(x,"Ctl",sep = "-"), min.pct = 0, logfc.threshold =0, only.pos = TRUE, test.use = "wilcox")
  cytokine <- res[res$p_val<0.05,] %>% row.names()
  x <- gsub("_",".",x)
  x <- gsub("\\+",".",x)
  test[cytokine,paste("cytokine",x,"CAS",sep = ".")] <- 1
}
test[is.na(test)] <- 0

test <- test[rowSums(test)>0,]
matrix <- matrix[row.names(test),]
test[test==1] <- "*"
test[test==0] <- ""

anno_col <- data.frame(Group = rep(levels(total_NK$disease),3),
                       Cell = rep(levels(total_NK$cluster_name_slim),each=2),
                       row.names = colnames(matrix))
anno_color <- list(Cell = c(NK_CD56bright = pal[1],
                            NK_CD56dim = pal[2],
                            NK_tr = pal[3]),
                   Group = c(Ctl = "#2b90d9", CAS = "#f26d5b"))

pheatmap::pheatmap(matrix,scale = "row",cluster_cols = F,
                   annotation_col = anno_col,
                   annotation_colors = anno_color,
                   display_numbers = test,
                   fontsize_number = 12,
                   show_colnames = F,
                   gaps_col = c(2,4),
                   treeheight_row = 20,
                   width = 12,
                   cellwidth = 16,cellheight = 12,
                   fontsize_row = 14,
                   filename = "Figure5_NK_cytokine_heatmap.pdf")


######################Figure 4 macropahge##############################

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

ggsave("Macrophage_marker_UMAP.pdf",p,width = 5,height = 4)

markers <- c("CD68","CD163","C1QA","APOE","CXCL8","CCL3")

p <- VlnPlot(total_macro,features = markers,stack = F,
             group.by = "cluster_name_slim",cols = pal[5:8],
             pt.size = 0,ncol = 3);p

ggsave("Macrophage_marker_VlnPlot.pdf",p,width = 9,height = 6)

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
ggsave("Macrophage_cell_proportion_bar.pdf",p,width = 6,height = 3)

#cell proportion box plot
cell.prop <- as.data.frame.matrix((table(total_macro$cluster_name_slim, total_macro$sample)))
cell.prop <- t(cell.prop)/colSums(cell.prop)
cell.prop <- as.data.frame(cell.prop)
cell.prop$Patient <- row.names(cell.prop)
cell.prop <- pivot_longer(cell.prop,cols = setdiff(colnames(cell.prop),"Patient"),names_to = "Cell",values_to = "Freq")

meta <- total_macro@meta.data[,c("sample",'disease')] %>% distinct()
row.names(meta) <- meta$sample
cell.prop$disease <- factor(meta[cell.prop$Patient,"disease"],
                            levels = c("Ctl","CAS"))
cell.prop$Cell <- factor(cell.prop$Cell,
                         levels = levels(total_macro$cluster_name_slim))

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
  facet_wrap(~Cell,ncol = 2,strip.position="bottom",scales = "free_y");p

ggsave("Macrophage_Cell proportion between Ctl and CAS.pdf",p,width = 5,height = 5)

#Roe
meta.tb <- total_macro@meta.data
meta.tb$disease <- as.character(meta.tb$disease) %>% factor(.,levels=c("Ctl","CAS"))
OR.list <- do.tissueDist(cellInfo.tb=meta.tb,
                         out.prefix="./Roe",
                         meta.cluster=meta.tb$cluster_name_slim,
                         loc=meta.tb$disease,
                         pdf.width=4,pdf.height=6,verbose=1)

matrix <- OR.list[["OR.dist.tb"]]
matrix <- column_to_rownames(matrix,var = "rid")
matrix <- matrix[levels(total_macro$cluster_name_slim),]
write.csv(matrix,file = "Figure5D_Roe_res_B.csv")

label <- matrix
label[label>1] <- "+++"
label[label<1&label>0.8] <- "++"
label[label<0.8&label>0.2] <- "+"
label[label<0.2&label>0] <- "+/-"

#DEGs
total_macro$disease <- factor(total_macro$disease,levels = c("Ctl","CAS"))
total_macro$cell_disease <- paste(total_macro$cluster_name_slim,total_macro$disease,sep = "-")
total_macro$cell_disease <- factor(total_macro$cell_disease,
                               levels = paste(rep(levels(total_macro$cluster_name_slim),each=2),rep(levels(total_macro$disease),4),sep = "-"))
total_macro <- SCTransform(total_macro)

Idents(total_macro) <- "cell_disease"

diffExp <- function(total,ident1,ident2){
  DefaultAssay(total) <- "SCT"
  Cell_SCT_diff <- FindMarkers(total, ident.1 = ident1, ident.2 = ident2, only.pos = F, test.use="wilcox",logfc.threshold = 0.2,min.pct = 0.4,min.diff.pct = 0.1,recorrect_umi = FALSE)
  return(Cell_SCT_diff)
}

cluster1 <- levels(total_macro$cluster_name_slim)
DefaultAssay(total_macro) <- "SCT"
data.features <- row.names(total_macro)
RB_genes = c(grep("^RPS", data.features, v=T),grep("^RPL", data.features, v=T))
MT_genes = grep("^MT-", data.features, v=T)
IG_genes = c(grep("^IGKV", data.features, v=T), grep("^IGHV",data.features,v=T),grep("^IGLV", data.features, v=T))
remove <- c(RB_genes,MT_genes,IG_genes)
total_macro_remove <- subset(total_macro,features = setdiff(data.features,remove))


dir.create("DEG_macro")
for(x in setdiff(cluster1,"Macro_CXCL8")){
  print(x)
  diff_res <- diffExp(total_macro_remove,paste0(x,"-CAS"),paste0(x,"-Ctl"))
  write.csv(diff_res,paste0("./DEG_macro/",x,"_diff_sig.csv"),row.names=T)
}

###Emapplot
library(clusterProfiler)
library(enrichplot)
library(ReactomePA)

GO_database <- 'org.Hs.eg.db'
cluster1 <- setdiff(levels(total_macro$cluster_name_slim),"Macro_CXCL8")

res <- list()
res2 <- c()
i=1
for(x in cluster1){
  print(x)
  diff_res <- read.csv(paste0("./DEG_macro/",x,"_diff_sig.csv"),row.names = 1)
  diff_res <- na.omit(diff_res)
  diff_gene <- diff_res[(diff_res$p_val<0.01)&(abs(diff_res$avg_log2FC)>0.5),]
  diff_gene <- diff_gene[!(row.names(diff_gene) %in% remove),]
  diff_gene <- diff_gene[order(diff_gene$avg_log2FC,decreasing = T),]
  diff_gene$change <- "Up"
  diff_gene$change[diff_gene$avg_log2FC<0] <- "Down"
  gene <- bitr(row.names(diff_gene[diff_gene$change=="Down",]),fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  res[[i]] <- gene$ENTREZID
  names(res)[i] <- x
  res2 <- c(res2,row.names(diff_gene))
  i=i+1
}

res2 <- unique(res2)

xx <- compareCluster(res, fun="enrichPathway")
xx2 <- pairwise_termsim(xx)
p <- emapplot(xx2,pie="Count",showCategory = 20,shadowtext=F)+
  scale_fill_manual(values = pal[5:8])+theme(legend.position = "none");p
ggsave("Macrophage_DEG_network_up_20.pdf",p,width = 7,height = 6)


#Cytokines
cytokine <- c("IL1A","IL1B","IL1RN","IL18","IL2","IL4","IL7","IL9","IL13","IL15","IL3",
              "IL5","CSF2","IL6","IL11","CSF3","IL12A","IL12B","LIF","OSM","IL10","IL20","TXLNA",
              "IL16","IL17A","IL17B","IFNA1","IFNB1","IFNG","CD40LG","LTB","TNF","TNFSF9","TNFSF13",
              "CD70","TNFSF8","FASLG","TNFSF18","TNFSF14","TNFSF4","TNFSF13B","TNFSF10",
              "TNFSF12","TNFSF11","TGFB1","TGFB2","TGFB3","EPO","TPO","FLT3LG","KITLG",
              "CSF1","MST1")

exp <- FetchData(total_macro,vars = cytokine)
exp <- exp[,colnames(exp)[!str_detect(colnames(exp),pattern = "rna")]]
total_macro[["cytokine"]] <- CreateAssayObject(t(exp))

#disease
total_macro_sub <- subset(total_macro,subset = cluster_name_slim != "Macro_CXCL8")

matrix <- AverageExpression(total_macro_sub,assays = "cytokine",group.by = "cell_disease") %>% as.data.frame()

test <- matrix(nrow = nrow(matrix),ncol = ncol(matrix)) %>%
  as.data.frame(row.names = row.names(matrix))
colnames(test) <- colnames(matrix)
for(x in unique(total_macro_sub$cluster_name_slim)){
  print(x)
  DefaultAssay(total_macro_sub) <- "cytokine"
  Idents(total_macro_sub) <- "cell_disease"
  res <- FindMarkers(total_macro_sub,ident.1 = paste(x,"CAS",sep = "-"),ident.2 = paste(x,"Ctl",sep = "-"), min.pct = 0, logfc.threshold =0, only.pos = TRUE, test.use = "wilcox")
  cytokine <- res[res$p_val<0.05,] %>% row.names()
  x <- gsub("_",".",x)
  x <- gsub(" ",".",x)
  x <- gsub("\\+",".",x)
  x <- gsub("\\(",".",x)
  x <- gsub("\\)",".",x)
  test[cytokine,paste("cytokine",x,"CAS",sep = ".")] <- 1
}
test[is.na(test)] <- 0

test <- test[rowSums(test)>0,]
matrix <- matrix[row.names(test),]
test[test==1] <- "*"
test[test==0] <- ""

anno_col <- data.frame(Group = rep(levels(total_macro_sub$disease),3),
                       Cell = rep(levels(total_macro_sub$cluster_name_slim),each=2),
                       row.names = colnames(matrix))
anno_color <- list(Cell = c(Macro_C1QA = pal[5],
                            Macro_APOE = pal[6],
                            Macro_CCL3 = pal[8]),
                   Group = c(Ctl = "#2b90d9", CAS = "#f26d5b"))

pheatmap::pheatmap(matrix,scale = "row",cluster_cols = F,
                   annotation_col = anno_col,
                   annotation_colors = anno_color,
                   display_numbers = test,
                   fontsize_number = 12,
                   show_colnames = F,
                   gaps_col = c(2,4),
                   width = 12,
                   cellwidth = 16,cellheight = 12,
                   fontsize_row = 14,
                   treeheight_row = 20,
                   filename = "Macro_cytokine_heatmap.pdf")



######################Figure 5 B cell##############################

total <- readRDS("total_final.2.rds")

total_B <- subset(total,idents=c("B cell","Plasma cell"))

data.list <- SplitObject(total_B,split.by = "sample")

data.list <- lapply(X=data.list,FUN=function(x){
  DefaultAssay(x) <- "RNA"
  x[['integrated']] <- NULL
  x <- SCTransform(x)
})

options(future.globals.maxSize= 100000000000000)

data.features <- SelectIntegrationFeatures(object.list = data.list)
IG_genes = c(grep("^IGKV", data.features, v=T), grep("^IGHV",data.features,v=T),grep("^IGLV", data.features, v=T))
data.features <- setdiff(data.features,IG_genes)
data.list <- PrepSCTIntegration(object.list = data.list,assay = "SCT", anchor.features = data.features,verbose = FALSE)
data.anchors <- FindIntegrationAnchors(object.list = data.list, normalization.method = "SCT", anchor.features = data.features, verbose = FALSE)
total_B <- IntegrateData(anchorset = data.anchors, normalization.method = "SCT", verbose = FALSE)

total_B <- Reduct(total_B)

DimPlot(total_B, reduction = "umap", label=T,group.by = "integrated_snn_res.0.5")

markers <- c("PTPRC","MS4A1","CD19","BANK1","JCHAIN","SDC1","IGKC","CD3D","CD4","CCR7","SELL","CD8A","GZMK","CCL5","KLRD1","GNLY","NKG7","MKI67","STMN1","TUBA1B","IL1B","S100A8","S100A9","CD68","C1QA","LYVE1","CD1E","LAMP3","FSCN1","PTGDS","IL3RA","LILRA4","F3","CXCL14","CD55")
Idents(total_B) <- "integrated_snn_res.0.5"
DefaultAssay(total_B) <- "SCT"
DotPlot(total_B, features = markers, dot.scale = 6)  + theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) + RotatedAxis()

total_B <- subset(total_B,subset = integrated_snn_res.0.5 %in% c(0:5,7:11))
total_B <- Reduct(total_B)

markers <- FindMarkers(total_B,ident.1 = c(2),logfc.threshold = 0.25,test.use = "wilcox",min.pct = 0.5,min.diff.pct = 0.2,only.pos = T,recorrect_umi = FALSE)

markers <- c("IL4R","IGHM","IGHD","CD27","CD38","GCSAM","MZB1")

DotPlot(total_B,features = markers,group.by = "integrated_snn_res.0.5")

total_B$cluster_name_slim <- "B_naive"
total_B$cluster_name_slim[total_B$integrated_snn_res.0.5 %in% c(4:6,9)] <- "B_memory(unswitched)"
total_B$cluster_name_slim[total_B$integrated_snn_res.0.5 %in% c(0,1)] <- "B_memory(switched)"
total_B$cluster_name_slim[total_B$integrated_snn_res.0.5 %in% c(8)] <- "Plasma cell"

total_B$cluster_name_slim <- factor(total_B$cluster_name_slim,
                                    levels = c("B_naive","B_memory(unswitched)","B_memory(switched)","Plasma cell"))
Idents(total_B) <- "cluster_name_slim"

pdf(file = "Figure5A_dimplot_B.pdf",width = 8,height = 6)
DimPlot_scCustom(total_B,colors_use = pal[1:10],figure_plot = T,label = T,label.size = 5,pt.size = 0.5)
dev.off()

total_B$disease <- factor(total_B$disease,levels = c("Ctl","CAS"))
Idents(total_B) <- "disease"
pdf(file = "Figure5C_dimplot_disease_B.pdf",width = 7,height = 6)
DimPlot_scCustom(total_B,reduction = "umap",label = F,pt.size = 0.1,colors_use = c("#2b90d9","#f26d5b"),figure_plot = T)
dev.off()

markers <- FindAllMarkers(total_B,logfc.threshold = 0.25,test.use = "wilcox",min.pct = 0.5,min.diff.pct = 0.2,only.pos = T,recorrect_umi = FALSE)
write.csv(markers,file = "total_markers_B.csv")

markers <- c("IL4R","IGHM","IGHD","CD27","CD38","MZB1")

p <- VlnPlot(total_B,features = markers,stack = T,group.by = "cluster_name_slim",cols = pal[1:4],fill.by = "ident")+
  theme(strip.text.x = element_text(angle = 0,face = "plain",size = 16),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 16),
        legend.position = "none");p
ggsave("Figure5B_ViolinPlot.pdf",p,width = 9,height = 4)

#cell proportion in bar plot

cell.prop <- table(total_B$cluster_name_slim,total_B$disease) %>% as.data.frame.matrix()
#cell.prop <- t(t(cell.prop)/colSums(cell.prop)) %>% as.data.frame()
cell.prop$cell <- row.names(cell.prop)
cell.prop <- pivot_longer(cell.prop,cols = c("Ctl","CAS"),names_to = "disease",values_to = "num")
cell.prop = ddply(cell.prop, 'cell', transform, percent=num/sum(num)*100)

cell.prop$cell <- factor(cell.prop$cell,levels = rev(levels(total_B$cluster_name_slim)))
cell.prop$disease <- factor(cell.prop$disease,levels = c("CAS","Ctl"))

index <- table(total_T$disease)
gap <- 100*index["Ctl"]/sum(index)

p <- ggplot(cell.prop,aes(x=percent,y=cell,fill=disease))+
  geom_bar(stat = 'identity',width = 0.8,colour='black')+
  geom_vline(aes(xintercept = gap), colour="black", linetype="dashed",linewidth = 1)+
  scale_fill_manual(values = c("#f26d5b","#2b90d9"))+
  theme_classic()+
  theme(axis.text = element_text(size = 12,colour = "black"),
        axis.title.x = element_text(size = 14,colour = "black"))+
  ylab("")+xlab("Cell ratio (%)");p
ggsave("Figure5D_cell_proportion_bar_Bcell.pdf",p,width = 6,height = 3)


#cell proportion box plot
cell.prop <- as.data.frame.matrix((table(total_B$cluster_name_slim, total_B$sample)))
cell.prop <- t(cell.prop)/colSums(cell.prop)
cell.prop <- as.data.frame(cell.prop)
cell.prop$Patient <- row.names(cell.prop)
cell.prop <- pivot_longer(cell.prop,cols = setdiff(colnames(cell.prop),"Patient"),names_to = "Cell",values_to = "Freq")

meta <- total_B@meta.data[,c("sample",'disease')] %>% distinct()
row.names(meta) <- meta$sample
cell.prop$disease <- factor(meta[cell.prop$Patient,"disease"],
                            levels = c("Ctl","CAS"))
cell.prop$Cell <- factor(cell.prop$Cell,
                         levels = levels(total_B$cluster_name_slim))

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
  facet_wrap(~Cell,ncol = 2,strip.position="bottom",scales = "free_y");p

ggsave("Figure5D_Cell proportion between Ctl and CAS_innate.pdf",p,width = 5,height = 5)

#Roe
meta.tb <- total_B@meta.data
meta.tb$disease <- as.character(meta.tb$disease) %>% factor(.,levels=c("Ctl","CAS"))
OR.list <- do.tissueDist(cellInfo.tb=meta.tb,
                         out.prefix="./Roe",
                         meta.cluster=meta.tb$cluster_name_slim,
                         loc=meta.tb$disease,
                         pdf.width=4,pdf.height=6,verbose=1)

matrix <- OR.list[["OR.dist.tb"]]
matrix <- column_to_rownames(matrix,var = "rid")
matrix <- matrix[levels(total_B$cluster_name_slim),]
write.csv(matrix,file = "Figure5D_Roe_res_B.csv")

label <- matrix
label[label>1] <- "+++"
label[label<1&label>0.8] <- "++"
label[label<0.8&label>0.2] <- "+"
label[label<0.2&label>0] <- "+/-"

#Correlation heatmap
library(psych)

total_SCT_mat <- as.matrix(total_B@assays$SCT@data)
total_SCT_mat_ClusterMeanSplit <- apply(total_SCT_mat,1,function(input){aggregate(input,list(total_B@meta.data$cluster_name_slim),mean)[,2]})
rownames(total_SCT_mat_ClusterMeanSplit) <- levels(total_B@meta.data$cluster_name_slim)
total_SCT_mat_ClusterMeanSplit <- t(total_SCT_mat_ClusterMeanSplit)
total_SCT_mat_ClusterMeanSplit <- total_SCT_mat_ClusterMeanSplit[rownames(total_SCT_mat),]

res <- corr.test(total_SCT_mat_ClusterMeanSplit)
r <- res$r
#bk <- seq(0.7,1,by=0.01)
pheatmap::pheatmap(r,cellwidth = 16,cellheight = 16,
                   width = 5,height = 4,
                   fontsize_col = 8,fontsize_row = 8,
                   #color = colorRampPalette(c("steelblue","white","chocolate1"))(length(bk)),
                   #legend_breaks=seq(0.75,1,0.05),breaks=bk,
                   treeheight_col = 20,treeheight_row = 20,
                   filename = "FigureS3B_correlation_Bcell.pdf")


#DEGs
total_B$disease <- factor(total_B$disease,levels = c("Ctl","CAS"))
total_B$cell_disease <- paste(total_B$cluster_name_slim,total_B$disease,sep = "-")
total_B$cell_disease <- factor(total_B$cell_disease,
                                    levels = paste(rep(levels(total_B$cluster_name_slim),each=2),rep(levels(total_B$disease),4),sep = "-"))


Idents(total_B) <- "cell_disease"

diffExp <- function(total,ident1,ident2){
  DefaultAssay(total) <- "SCT"
  Cell_SCT_diff <- FindMarkers(total, ident.1 = ident1, ident.2 = ident2, only.pos = F, test.use="wilcox",logfc.threshold = 0.2,min.pct = 0.4,min.diff.pct = 0.1,recorrect_umi = FALSE)
  return(Cell_SCT_diff)
}

cluster1 <- levels(total_B$cluster_name_slim)
DefaultAssay(total_B) <- "SCT"
data.features <- row.names(total_B)
RB_genes = c(grep("^RPS", data.features, v=T),grep("^RPL", data.features, v=T))
MT_genes = grep("^MT-", data.features, v=T)
IG_genes = c(grep("^IGKV", data.features, v=T), grep("^IGHV",data.features,v=T),grep("^IGLV", data.features, v=T))
remove <- c(RB_genes,MT_genes,IG_genes)
total_B_remove <- subset(total_B,features = setdiff(data.features,remove))


dir.create("DEG_B")
for(x in cluster1){
  print(x)
  diff_res <- diffExp(total_B_remove,paste0(x,"-CAS"),paste0(x,"-Ctl"))
  write.csv(diff_res,paste0("./DEG_B/",x,"_diff_sig.csv"),row.names=T)
}

###Emapplot
library(clusterProfiler)
library(enrichplot)
library(ReactomePA)

GO_database <- 'org.Hs.eg.db'
cluster1 <- levels(total_B$cluster_name_slim)

res <- list()
res2 <- c()
i=1
for(x in cluster1){
  print(x)
  diff_res <- read.csv(paste0("./DEG_B/",x,"_diff_sig.csv"),row.names = 1)
  diff_res <- na.omit(diff_res)
  diff_gene <- diff_res[(diff_res$p_val<0.01)&(abs(diff_res$avg_log2FC)>0.5),]
  diff_gene <- diff_gene[!(row.names(diff_gene) %in% remove),]
  diff_gene <- diff_gene[order(diff_gene$avg_log2FC,decreasing = T),]
  diff_gene$change <- "Up"
  diff_gene$change[diff_gene$avg_log2FC<0] <- "Down"
  gene <- bitr(row.names(diff_gene[diff_gene$change=="Up",]),fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  res[[i]] <- gene$ENTREZID
  names(res)[i] <- x
  res2 <- c(res2,row.names(diff_gene))
  i=i+1
}

res2 <- unique(res2)

xx <- compareCluster(res, fun="enrichPathway")
xx2 <- pairwise_termsim(xx)
p <- emapplot(xx2,pie="Count",showCategory = 20,shadowtext=F)+
  scale_fill_manual(values = pal[1:4])+theme(legend.position = "none");p
ggsave("Figure5E_DEG_Reactome_network_up_20.pdf",p,width = 7,height = 6)

#DEGs_Smoke
total_B$condition <- factor(total_B$condition,levels = c("Ctl","CAS_NS","CAS_Smo"))
total_B$cell_condition <- paste(total_B$cluster_name_slim,total_B$condition,sep = "-")
total_B$cell_condition <- factor(total_B$cell_condition,
                               levels = paste(rep(levels(total_B$cluster_name_slim),each=3),rep(levels(total_B$condition),4),sep = "-"))


Idents(total_B) <- "cell_condition"

diffExp <- function(total,ident1,ident2){
  DefaultAssay(total) <- "SCT"
  Cell_SCT_diff <- FindMarkers(total, ident.1 = ident1, ident.2 = ident2, only.pos = F, test.use="wilcox",logfc.threshold = 0.2,min.pct = 0.4,min.diff.pct = 0.1,recorrect_umi = FALSE)
  return(Cell_SCT_diff)
}

cluster1 <- levels(total_B$cluster_name_slim)
DefaultAssay(total_B) <- "SCT"
data.features <- row.names(total_B)
RB_genes = c(grep("^RPS", data.features, v=T),grep("^RPL", data.features, v=T))
MT_genes = grep("^MT-", data.features, v=T)
IG_genes = c(grep("^IGKV", data.features, v=T), grep("^IGHV",data.features,v=T),grep("^IGLV", data.features, v=T))
remove <- c(RB_genes,MT_genes,IG_genes)
total_B_remove <- subset(total_B,features = setdiff(data.features,remove))


dir.create("DEG_B_Smoke")
for(x in cluster1){
  print(x)
  diff_res <- diffExp(total_B_remove,paste0(x,"-CAS_Smo"),paste0(x,"-CAS_NS"))
  write.csv(diff_res,paste0("./DEG_B_Smoke/",x,"_diff_sig.csv"),row.names=T)
}

###Emapplot_Smoke
library(clusterProfiler)
library(enrichplot)
library(ReactomePA)

GO_database <- 'org.Hs.eg.db'
cluster1 <- levels(total_B$cluster_name_slim)

res <- list()
res2 <- c()
i=1
for(x in cluster1){
  print(x)
  diff_res <- read.csv(paste0("./DEG_B_Smoke/",x,"_diff_sig.csv"),row.names = 1)
  diff_res <- na.omit(diff_res)
  diff_gene <- diff_res[(diff_res$p_val<0.01)&(abs(diff_res$avg_log2FC)>0.5),]
  diff_gene <- diff_gene[!(row.names(diff_gene) %in% remove),]
  diff_gene <- diff_gene[order(diff_gene$avg_log2FC,decreasing = T),]
  diff_gene$change <- "Up"
  diff_gene$change[diff_gene$avg_log2FC<0] <- "Down"
  gene <- bitr(row.names(diff_gene[diff_gene$change=="Up",]),fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  res[[i]] <- gene$ENTREZID
  names(res)[i] <- x
  res2 <- c(res2,row.names(diff_gene))
  i=i+1
}

res2 <- unique(res2)

xx <- compareCluster(res, fun="enrichKEGG")
xx2 <- pairwise_termsim(xx)
p <- emapplot(xx2,pie="Count",showCategory = 20,shadowtext=F)+
  scale_fill_manual(values = pal[1:4])+theme(legend.position = "none");p
ggsave("FigureS3_DEG_Bcell_network_up_20.pdf",p,width = 7,height = 6)



#Cytokines
cytokine <- c("IL1A","IL1B","IL1RN","IL18","IL2","IL4","IL7","IL9","IL13","IL15","IL3",
              "IL5","CSF2","IL6","IL11","CSF3","IL12A","IL12B","LIF","OSM","IL10","IL20","TXLNA",
              "IL16","IL17A","IL17B","IFNA1","IFNB1","IFNG","CD40LG","LTB","TNF","TNFSF9","TNFSF13",
              "CD70","TNFSF8","FASLG","TNFSF18","TNFSF14","TNFSF4","TNFSF13B","TNFSF10",
              "TNFSF12","TNFSF11","TGFB1","TGFB2","TGFB3","EPO","TPO","FLT3LG","KITLG",
              "CSF1","MST1")

exp <- FetchData(total_B,vars = cytokine)
exp <- exp[,colnames(exp)[!str_detect(colnames(exp),pattern = "rna")]]
total_B[["cytokine"]] <- CreateAssayObject(t(exp))

#disease
matrix <- AverageExpression(total_B,assays = "cytokine",group.by = "cell_disease") %>% as.data.frame()

test <- matrix(nrow = nrow(matrix),ncol = ncol(matrix)) %>%
  as.data.frame(row.names = row.names(matrix))
colnames(test) <- colnames(matrix)
for(x in unique(total_B$cluster_name_slim)){
  print(x)
  DefaultAssay(total_B) <- "cytokine"
  Idents(total_B) <- "cell_disease"
  res <- FindMarkers(total_B,ident.1 = paste(x,"CAS",sep = "-"),ident.2 = paste(x,"Ctl",sep = "-"), min.pct = 0, logfc.threshold =0, only.pos = TRUE, test.use = "wilcox")
  cytokine <- res[res$p_val<0.05,] %>% row.names()
  x <- gsub("_",".",x)
  x <- gsub(" ",".",x)
  x <- gsub("\\+",".",x)
  x <- gsub("\\(",".",x)
  x <- gsub("\\)",".",x)
  test[cytokine,paste("cytokine",x,"CAS",sep = ".")] <- 1
}
test[is.na(test)] <- 0

test <- test[rowSums(test)>0,]
matrix <- matrix[row.names(test),]
test[test==1] <- "*"
test[test==0] <- ""

anno_col <- data.frame(Group = rep(levels(total_B$disease),4),
                       Cell = rep(levels(total_B$cluster_name_slim),each=2),
                       row.names = colnames(matrix))
anno_color <- list(Cell = c(B_naive = pal[1],
                            "B_memory(unswitched)" = pal[2],
                            "B_memory(switched)" = pal[3],
                            "Plasma cell" = pal[4]),
                   Group = c(Ctl = "#2b90d9", CAS = "#f26d5b"))

pheatmap::pheatmap(matrix,scale = "row",cluster_cols = F,
                   annotation_col = anno_col,
                   annotation_colors = anno_color,
                   display_numbers = test,
                   fontsize_number = 12,
                   show_colnames = F,
                   gaps_col = seq(2,6,2),
                   width = 12,
                   cellwidth = 16,cellheight = 12,
                   fontsize_row = 14,
                   treeheight_row = 20,
                   filename = "Figure5F_cytokine_heatmap_B.pdf")

#condition
total_B_sub <- subset(total_B,subset = condition != "Ctl")
total_B_sub$condition <- droplevels(total_B_sub$condition)
matrix <- AverageExpression(total_B_sub,assays = "cytokine",group.by = "cell_condition") %>% as.data.frame()

test <- matrix(nrow = nrow(matrix),ncol = ncol(matrix)) %>%
  as.data.frame(row.names = row.names(matrix))
colnames(test) <- colnames(matrix)
for(x in unique(total_B_sub$cluster_name_slim)){
  print(x)
  DefaultAssay(total_B_sub) <- "cytokine"
  Idents(total_B_sub) <- "cell_condition"
  res <- FindMarkers(total_B_sub,ident.1 = paste(x,"CAS_Smo",sep = "-"),ident.2 = paste(x,"CAS_NS",sep = "-"), min.pct = 0, logfc.threshold =0, only.pos = TRUE, test.use = "wilcox")
  cytokine <- res[res$p_val<0.05,] %>% row.names()
  x <- gsub("_",".",x)
  x <- gsub(" ",".",x)
  x <- gsub("\\+",".",x)
  x <- gsub("\\(",".",x)
  x <- gsub("\\)",".",x)
  test[cytokine,paste("cytokine",x,"CAS.Smo",sep = ".")] <- 1
}
test[is.na(test)] <- 0

test <- test[rowSums(test)>0,]
matrix <- matrix[row.names(test),]
test[test==1] <- "*"
test[test==0] <- ""

anno_col <- data.frame(Group = rep(levels(total_B_sub$condition),4),
                       Cell = rep(levels(total_B_sub$cluster_name_slim),each=2),
                       row.names = colnames(matrix))
anno_color <- list(Cell = c(B_naive = pal[1],
                            "B_memory(unswitched)" = pal[2],
                            "B_memory(switched)" = pal[3],
                            "Plasma cell" = pal[4]),
                   Group = c(CAS_NS = "#2b90d9", CAS_Smo = "#f26d5b"))

pheatmap::pheatmap(matrix,scale = "row",cluster_cols = F,
                   annotation_col = anno_col,
                   annotation_colors = anno_color,
                   display_numbers = test,
                   fontsize_number = 12,
                   show_colnames = F,
                   gaps_col = seq(2,6,2),
                   width = 12,
                   cellwidth = 16,cellheight = 12,
                   fontsize_row = 14,
                   treeheight_row = 20,
                   filename = "FigureS3_cytokine_heatmap_B.pdf")

######################Figure 5 T cell##############################

total <- readRDS("total_final.2.rds")

total_T <- subset(total,idents=c("CD4T cell","CD8T cell"))

data.list <- SplitObject(total_T,split.by = "sample")

data.list <- lapply(X=data.list,FUN=function(x){
  DefaultAssay(x) <- "RNA"
  x[['integrated']] <- NULL
  x <- SCTransform(x)
})

options(future.globals.maxSize= 100000000000000)

data.features <- SelectIntegrationFeatures(object.list = data.list)
nms <- data.features
tr_genes = c(grep("^TRAV", nms, v=T),grep("^TRBV",nms,v=T),grep("^TRDV", nms, v=T),grep("^TRGV", nms, v=T),grep("^TRAJ", nms, v=T),grep("^TRBJ", nms, v=T))
data.features <- setdiff(data.features,tr_genes)
data.list <- PrepSCTIntegration(object.list = data.list,assay = "SCT", anchor.features = data.features,verbose = FALSE)
data.anchors <- FindIntegrationAnchors(object.list = data.list, normalization.method = "SCT", anchor.features = data.features, verbose = FALSE)
total_T <- IntegrateData(anchorset = data.anchors, normalization.method = "SCT", verbose = FALSE)
saveRDS(total_T,file = "total_T.rds")

total_T <- Reduct(total_T)

DimPlot(total_T, reduction = "umap", label=T,group.by = "integrated_snn_res.0.5")

markers <- c("PTPRC","MS4A1","CD19","BANK1","JCHAIN","SDC1","IGKC","CD3D","CD4","CCR7","SELL","CD8A","GZMK","CCL5","KLRD1","GNLY","NKG7","MKI67","STMN1","TUBA1B","IL1B","S100A8","S100A9","CD68","C1QA","LYVE1","CD1E","LAMP3","FSCN1","PTGDS","IL3RA","LILRA4","F3","CXCL14","CD55")
Idents(total_T) <- "integrated_snn_res.0.5"
DefaultAssay(total_T) <- "SCT"
DotPlot(total_T, features = markers, dot.scale = 6)  + theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) + RotatedAxis()

total_T <- subset(total_T,subset = integrated_snn_res.0.5 %in% c(0:14,16))
total_T <- Reduct(total_T)

marker <- c("CD3D","CD8A","CD8B","CD4","GATA3","RORC","TRDC","FOXP3","IL2RA","CCR7","LEF1","SELL","TCF7","CD27","CD28","S1PR1","IL7R","PRF1","GZMA","CCL5","GPR183","KLRG1","CX3CR1","FCGR3A","FGFBP2","GZMH","TBX21","EOMES","S1PR5","GZMK","CXCR4","CXCR3","CD44","CD6","XCL1","XCL2","MYADM","CAPG","RORA","NR4A1","NR4A2","CD69","ITGAE","ZNF683","CD160","KIR2DL4","TMIGD2","KLRC1","IKZF2","ENTPD1","HAVCR2","CXCL13","PDCD1","LAYN","TOX","IFNG","TNF","GZMB","MKI67","STMN1","SLC4A10","KLRB1","ZBTB16","NCR3","RORC","RORA")
DotPlot(total_T, features = unique(marker), dot.scale = 6)  + theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) + RotatedAxis()

markers <- FindAllMarkers(total_T,logfc.threshold = 0.25,test.use = "wilcox",min.pct = 0.5,min.diff.pct = 0.2,only.pos = T,recorrect_umi = FALSE)

total_T$cluster_name_slim <- "CD4Tmemory"
total_T$cluster_name_slim[total_T$integrated_snn_res.0.5 %in% c(12)] <- "CD8Tmemory"
total_T$cluster_name_slim[total_T$integrated_snn_res.0.5 %in% c(6,7)] <- "CD8Teffector"
total_T$cluster_name_slim[total_T$integrated_snn_res.0.5 %in% c(9)] <- "CD4Treg"
total_T$cluster_name_slim[total_T$integrated_snn_res.0.5 %in% c(13)] <- "CD4Texhausted"
total_T$cluster_name_slim[total_T$integrated_snn_res.0.5 %in% c(14)] <- "CD4Tstress"

total_T$cluster_name_slim <- factor(total_T$cluster_name_slim,
                                    levels = c("CD4Tmemory",'CD4Texhausted',"CD4Tstress","CD4Treg","CD8Tmemory","CD8Teffector"))

Idents(total_T) <- "cluster_name_slim"

pdf(file = "Figure5G_dimplot_T.pdf",width = 8,height = 6)
DimPlot_scCustom(total_T,colors_use = hue_pal()(6),figure_plot = T,label = T,label.size = 5,pt.size = 0.5)
dev.off()

total_T$disease <- factor(total_T$disease,levels = c("Ctl","CAS"))
Idents(total_T) <- "disease"
pdf(file = "Figure5I_dimplot_disease_T.pdf",width = 7,height = 6)
DimPlot_scCustom(total_T,reduction = "umap",label = F,pt.size = 0.1,colors_use = c("#2b90d9","#f26d5b"),figure_plot = T)
dev.off()

markers <- FindAllMarkers(total_T,logfc.threshold = 0.25,test.use = "wilcox",min.pct = 0.5,min.diff.pct = 0.2,only.pos = T,recorrect_umi = FALSE)
write.csv(markers,file = "total_markers_T.csv")

markers <- c("CD3D",'CD4',"CD8A","IL7R","TCF7","PDCD1","HSPA1A","FOXP3","GZMA")
p <- VlnPlot(total_T,features = markers,pt.size = 0,ncol = 3)
ggsave("Figure5H_VlnPlot_T.pdf",p,width = 8,height = 8)

saveRDS(total_T,file = "total_T.rds")

#cell proportion in bar plot

cell.prop <- table(total_T$cluster_name_slim,total_T$disease) %>% as.data.frame.matrix()
#cell.prop <- t(t(cell.prop)/colSums(cell.prop)) %>% as.data.frame()
cell.prop$cell <- row.names(cell.prop)
cell.prop <- pivot_longer(cell.prop,cols = c("Ctl","CAS"),names_to = "disease",values_to = "num")
cell.prop = ddply(cell.prop, 'cell', transform, percent=num/sum(num)*100)

cell.prop$cell <- factor(cell.prop$cell,levels = rev(levels(total_T$cluster_name_slim)))
cell.prop$disease <- factor(cell.prop$disease,levels = c("CAS","Ctl"))

index <- table(total_T$disease)
gap <- 100*index["Ctl"]/sum(index)

p <- ggplot(cell.prop,aes(x=percent,y=cell,fill=disease))+
  geom_bar(stat = 'identity',width = 0.8,colour='black')+
  geom_vline(aes(xintercept = gap), colour="black", linetype="dashed",linewidth = 1)+
  scale_fill_manual(values = c("#f26d5b","#2b90d9"))+
  theme_classic()+
  theme(axis.text = element_text(size = 12,colour = "black"),
        axis.title.x = element_text(size = 14,colour = "black"))+
  ylab("")+xlab("Cell ratio (%)");p
ggsave("Figure5J_cell_proportion_bar_Tcell.pdf",p,width = 5,height = 3.5)


#cell proportion box plot
cell.prop <- as.data.frame.matrix((table(total_T$cluster_name_slim, total_T$sample)))
cell.prop <- t(cell.prop)/colSums(cell.prop)
cell.prop <- as.data.frame(cell.prop)
cell.prop$Patient <- row.names(cell.prop)
cell.prop <- pivot_longer(cell.prop,cols = setdiff(colnames(cell.prop),"Patient"),names_to = "Cell",values_to = "Freq")

meta <- total_T@meta.data[,c("sample",'disease')] %>% distinct()
row.names(meta) <- meta$sample
cell.prop$disease <- factor(meta[cell.prop$Patient,"disease"],
                            levels = c("Ctl","CAS"))
cell.prop$Cell <- factor(cell.prop$Cell,
                         levels = levels(total_T$cluster_name_slim))

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
  facet_wrap(~Cell,ncol = 3,strip.position="bottom",scales = "free_y");p

ggsave("Figure5J_Cell proportion between Ctl and CAS_T.pdf",p,width = 7,height = 5)

#Roe
meta.tb <- total_T@meta.data
meta.tb$disease <- as.character(meta.tb$disease) %>% factor(.,levels=c("Ctl","CAS"))
OR.list <- do.tissueDist(cellInfo.tb=meta.tb,
                         out.prefix="./Roe",
                         meta.cluster=meta.tb$cluster_name_slim,
                         loc=meta.tb$disease,
                         pdf.width=4,pdf.height=6,verbose=1)

matrix <- OR.list[["OR.dist.tb"]]
matrix <- column_to_rownames(matrix,var = "rid")
matrix <- matrix[levels(total_T$cluster_name_slim),]
write.csv(matrix,file = "Figure5J_Roe_res_T.csv")

label <- matrix
label[label>1] <- "+++"
label[label<1&label>0.8] <- "++"
label[label<0.8&label>0.2] <- "+"
label[label<0.2&label>0] <- "+/-"

#stacked bar plot
total_T <- readRDS("total_T.rds")

total_T$sample <- factor(total_T$sample,
                         levels = c("Ctl1","Ctl2","Ctl3","CAS1","CAS2","CAS3","CAS4","CAS5","CAS6","CAS7"))

cell.prop <- as.data.frame(prop.table(table(total_T$cluster_name_slim, total_T$sample)))
colnames(cell.prop) <- c("cluster","patient","proportion")
cell.prop$cluster <- factor(cell.prop$cluster,levels = levels(total_T$cluster_name_slim))
cell.prop$patient <- factor(cell.prop$patient,levels = rev(levels(total_T$sample)))

p <- ggplot(cell.prop,aes(proportion,patient,fill=cluster))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  scale_fill_manual(values=hue_pal()(6))+
  theme_bw()+
  theme_classic()+
  theme(axis.ticks.length=unit(0.2,'cm'),
        axis.text.x = element_text(size=16,colour = "black"),
        axis.text.y = element_text(size=16,colour = "black"),
        axis.title.x = element_text(size=16,colour = "black"),
        axis.title.y = element_blank(),
        legend.text = element_text(size=14,colour = "black"))+
  guides(fill=guide_legend(title=NULL))+
  xlab("Cell proportion");p
ggsave("Cell_proportion_T_sample.pdf",p,width = 7,height = 5)

cell.prop <- as.data.frame.matrix(table(total_T$cluster_name_slim, total_T$sample))
cell.prop <- t(t(cell.prop)/colSums(cell.prop)) %>% as.data.frame()


#Correlation heatmap
library(psych)

total_SCT_mat <- as.matrix(total_T@assays$SCT@data)
total_SCT_mat_ClusterMeanSplit <- apply(total_SCT_mat,1,function(input){aggregate(input,list(total_T@meta.data$cluster_name_slim),mean)[,2]})
rownames(total_SCT_mat_ClusterMeanSplit) <- levels(total_T@meta.data$cluster_name_slim)
total_SCT_mat_ClusterMeanSplit <- t(total_SCT_mat_ClusterMeanSplit)
total_SCT_mat_ClusterMeanSplit <- total_SCT_mat_ClusterMeanSplit[rownames(total_SCT_mat),]

res <- corr.test(total_SCT_mat_ClusterMeanSplit)
r <- res$r
#bk <- seq(0.7,1,by=0.01)
pheatmap::pheatmap(r,cellwidth = 16,cellheight = 16,
                   width = 5,height = 4,
                   fontsize_col = 8,fontsize_row = 8,
                   #color = colorRampPalette(c("steelblue","white","chocolate1"))(length(bk)),
                   #legend_breaks=seq(0.75,1,0.05),breaks=bk,
                   treeheight_col = 20,treeheight_row = 20,
                   filename = "Figure5L_correlation_Tcell.pdf")

#DEGs
total_T$disease <- factor(total_T$disease,levels = c("Ctl","CAS"))
total_T$cell_disease <- paste(total_T$cluster_name_slim,total_T$disease,sep = "-")
total_T$cell_disease <- factor(total_T$cell_disease,
                               levels = paste(rep(levels(total_T$cluster_name_slim),each=2),rep(levels(total_T$disease),6),sep = "-"))


Idents(total_T) <- "cell_disease"

diffExp <- function(total,ident1,ident2){
  DefaultAssay(total) <- "SCT"
  Cell_SCT_diff <- FindMarkers(total, ident.1 = ident1, ident.2 = ident2, only.pos = F, test.use="wilcox",logfc.threshold = 0.2,min.pct = 0.4,min.diff.pct = 0.1,recorrect_umi = FALSE)
  return(Cell_SCT_diff)
}

cluster1 <- levels(total_T$cluster_name_slim)
DefaultAssay(total_T) <- "SCT"
data.features <- row.names(total_T)
RB_genes = c(grep("^RPS", data.features, v=T),grep("^RPL", data.features, v=T))
MT_genes = grep("^MT-", data.features, v=T)
tr_genes = c(grep("^TRAV", data.features, v=T),grep("^TRBV",data.features,v=T),grep("^TRDV", data.features, v=T),grep("^TRGV", data.features, v=T),grep("^TRAJ", data.features, v=T),grep("^TRBJ", data.features, v=T))
remove <- c(RB_genes,MT_genes,tr_genes)
total_T_remove <- subset(total_T,features = setdiff(data.features,remove))


dir.create("DEG_T")
for(x in cluster1){
  print(x)
  diff_res <- diffExp(total_T_remove,paste0(x,"-CAS"),paste0(x,"-Ctl"))
  write.csv(diff_res,paste0("./DEG_T/",x,"_diff_sig.csv"),row.names=T)
}

###Emapplot
library(clusterProfiler)
library(enrichplot)
library(ReactomePA)

GO_database <- 'org.Hs.eg.db'
cluster1 <- levels(total_T$cluster_name_slim)

res <- list()
res2 <- c()
i=1
for(x in cluster1){
  print(x)
  diff_res <- read.csv(paste0("./DEG_T/",x,"_diff_sig.csv"),row.names = 1)
  diff_res <- na.omit(diff_res)
  diff_gene <- diff_res[(diff_res$p_val<0.01)&(abs(diff_res$avg_log2FC)>0.5),]
  diff_gene <- diff_gene[!(row.names(diff_gene) %in% remove),]
  diff_gene <- diff_gene[order(diff_gene$avg_log2FC,decreasing = T),]
  diff_gene$change <- "Up"
  diff_gene$change[diff_gene$avg_log2FC<0] <- "Down"
  gene <- bitr(row.names(diff_gene[diff_gene$change=="Up",]),fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  res[[i]] <- gene$ENTREZID
  names(res)[i] <- x
  res2 <- c(res2,row.names(diff_gene))
  i=i+1
}

res2 <- unique(res2)

xx <- compareCluster(res, fun="enrichPathway")
xx2 <- pairwise_termsim(xx)
p <- emapplot(xx2,pie="Count",showCategory = 10,shadowtext=F)+
  scale_fill_manual(values = hue_pal()(6))+theme(legend.position = "none");p
ggsave("Figure5M_DEG_Reactome_network_up_10.pdf",p,width = 7,height = 6)

#DEGs-Smoke
total_T$condition <- factor(total_T$condition,levels = c("Ctl","CAS_NS","CAS_Smo"))
total_T$cell_condition <- paste(total_T$cluster_name_slim,total_T$condition,sep = "-")
total_T$cell_condition <- factor(total_T$cell_condition,
                               levels = paste(rep(levels(total_T$cluster_name_slim),each=3),rep(levels(total_T$condition),6),sep = "-"))


Idents(total_T) <- "cell_condition"

diffExp <- function(total,ident1,ident2){
  DefaultAssay(total) <- "SCT"
  Cell_SCT_diff <- FindMarkers(total, ident.1 = ident1, ident.2 = ident2, only.pos = F, test.use="wilcox",logfc.threshold = 0.2,min.pct = 0.4,min.diff.pct = 0.1,recorrect_umi = FALSE)
  return(Cell_SCT_diff)
}

cluster1 <- levels(total_T$cluster_name_slim)
DefaultAssay(total_T) <- "SCT"
data.features <- row.names(total_T)
RB_genes = c(grep("^RPS", data.features, v=T),grep("^RPL", data.features, v=T))
MT_genes = grep("^MT-", data.features, v=T)
tr_genes = c(grep("^TRAV", data.features, v=T),grep("^TRBV",data.features,v=T),grep("^TRDV", data.features, v=T),grep("^TRGV", data.features, v=T),grep("^TRAJ", data.features, v=T),grep("^TRBJ", data.features, v=T))
remove <- c(RB_genes,MT_genes,tr_genes)
total_T_remove <- subset(total_T,features = setdiff(data.features,remove))


dir.create("DEG_T_Smoke")
for(x in cluster1){
  print(x)
  diff_res <- diffExp(total_T_remove,paste0(x,"-CAS_Smo"),paste0(x,"-CAS_NS"))
  write.csv(diff_res,paste0("./DEG_T_Smoke/",x,"_diff_sig.csv"),row.names=T)
}

###Emapplot-Smoke
library(clusterProfiler)
library(enrichplot)
library(ReactomePA)

GO_database <- 'org.Hs.eg.db'
cluster1 <- levels(total_T$cluster_name_slim)

res <- list()
res2 <- c()
i=1
for(x in cluster1){
  print(x)
  diff_res <- read.csv(paste0("./DEG_T_Smoke/",x,"_diff_sig.csv"),row.names = 1)
  diff_res <- na.omit(diff_res)
  diff_gene <- diff_res[(diff_res$p_val<0.01)&(abs(diff_res$avg_log2FC)>0.5),]
  diff_gene <- diff_gene[!(row.names(diff_gene) %in% remove),]
  diff_gene <- diff_gene[order(diff_gene$avg_log2FC,decreasing = T),]
  diff_gene$change <- "Up"
  diff_gene$change[diff_gene$avg_log2FC<0] <- "Down"
  gene <- bitr(row.names(diff_gene[diff_gene$change=="Down",]),fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
  res[[i]] <- gene$ENTREZID
  names(res)[i] <- x
  res2 <- c(res2,row.names(diff_gene))
  i=i+1
}

res2 <- unique(res2)

xx <- compareCluster(res, fun="enrichPathway")
xx2 <- pairwise_termsim(xx)
p <- emapplot(xx2,pie="Count",showCategory = 10,shadowtext=F)+
  scale_fill_manual(values = hue_pal()(6))+theme(legend.position = "none");p
ggsave("FigureS3_DEG_Tcell_network_down_10.pdf",p,width = 7,height = 6)


#Cytokines
cytokine <- c("IL1A","IL1B","IL1RN","IL18","IL2","IL4","IL7","IL9","IL13","IL15","IL3",
              "IL5","CSF2","IL6","IL11","CSF3","IL12A","IL12B","LIF","OSM","IL10","IL20","TXLNA",
              "IL16","IL17A","IL17B","IFNA1","IFNB1","IFNG","CD40LG","LTB","TNF","TNFSF9","TNFSF13",
              "CD70","TNFSF8","FASLG","TNFSF18","TNFSF14","TNFSF4","TNFSF13B","TNFSF10",
              "TNFSF12","TNFSF11","TGFB1","TGFB2","TGFB3","EPO","TPO","FLT3LG","KITLG",
              "CSF1","MST1")

exp <- FetchData(total_T,vars = cytokine)
exp <- exp[,colnames(exp)[!str_detect(colnames(exp),pattern = "rna")]]
total_T[["cytokine"]] <- CreateAssayObject(t(exp))

#disease
matrix <- AverageExpression(total_T,assays = "cytokine",group.by = "cell_disease") %>% as.data.frame()

test <- matrix(nrow = nrow(matrix),ncol = ncol(matrix)) %>%
  as.data.frame(row.names = row.names(matrix))
colnames(test) <- colnames(matrix)
for(x in unique(total_T$cluster_name_slim)){
  print(x)
  DefaultAssay(total_T) <- "cytokine"
  Idents(total_T) <- "cell_disease"
  res <- FindMarkers(total_T,ident.1 = paste(x,"CAS",sep = "-"),ident.2 = paste(x,"Ctl",sep = "-"), min.pct = 0, logfc.threshold =0, only.pos = TRUE, test.use = "wilcox")
  cytokine <- res[res$p_val<0.05,] %>% row.names()
  x <- gsub("_",".",x)
  x <- gsub(" ",".",x)
  x <- gsub("\\+",".",x)
  x <- gsub("\\(",".",x)
  x <- gsub("\\)",".",x)
  test[cytokine,paste("cytokine",x,"CAS",sep = ".")] <- 1
}
test[is.na(test)] <- 0

test <- test[rowSums(test)>0,]
matrix <- matrix[row.names(test),]
test[test==1] <- "*"
test[test==0] <- ""

anno_col <- data.frame(Group = rep(levels(total_T$disease),6),
                       Cell = rep(levels(total_T$cluster_name_slim),each=2),
                       row.names = colnames(matrix))
anno_color <- list(Cell = c(CD4Tmemory = hue_pal()(6)[1],
                            "CD4Texhausted" = hue_pal()(6)[2],
                            "CD4Tstress" = hue_pal()(6)[3],
                            "CD4Treg" = hue_pal()(6)[4],
                            "CD8Tmemory" = hue_pal()(6)[5],
                            "CD8Teffector" = hue_pal()(6)[6]),
                   
                   Group = c(Ctl = "#2b90d9", CAS = "#f26d5b"))

pheatmap::pheatmap(matrix,scale = "row",cluster_cols = F,
                   annotation_col = anno_col,
                   annotation_colors = anno_color,
                   display_numbers = test,
                   fontsize_number = 12,
                   show_colnames = F,
                   gaps_col = seq(2,10,2),
                   width = 12,
                   cellwidth = 16,cellheight = 12,
                   fontsize_row = 14,
                   treeheight_row = 20,
                   filename = "Figure5H_cytokine_heatmap_T.pdf")


#condition
total_T_sub <- subset(total_T,subset = condition != "Ctl")
total_T_sub$condition <- droplevels(total_T_sub$condition)
matrix <- AverageExpression(total_T_sub,assays = "cytokine",group.by = "cell_condition") %>% as.data.frame()

test <- matrix(nrow = nrow(matrix),ncol = ncol(matrix)) %>%
  as.data.frame(row.names = row.names(matrix))
colnames(test) <- colnames(matrix)
for(x in unique(total_T_sub$cluster_name_slim)){
  print(x)
  DefaultAssay(total_T_sub) <- "cytokine"
  Idents(total_T_sub) <- "cell_condition"
  res <- FindMarkers(total_T_sub,ident.1 = paste(x,"CAS_Smo",sep = "-"),ident.2 = paste(x,"CAS_NS",sep = "-"), min.pct = 0, logfc.threshold =0, only.pos = TRUE, test.use = "wilcox")
  cytokine <- res[res$p_val<0.05,] %>% row.names()
  x <- gsub("_",".",x)
  x <- gsub(" ",".",x)
  x <- gsub("\\+",".",x)
  x <- gsub("\\(",".",x)
  x <- gsub("\\)",".",x)
  test[cytokine,paste("cytokine",x,"CAS.Smo",sep = ".")] <- 1
}
test[is.na(test)] <- 0

test <- test[rowSums(test)>0,]
matrix <- matrix[row.names(test),]
test[test==1] <- "*"
test[test==0] <- ""

anno_col <- data.frame(Group = rep(levels(total_T_sub$condition),6),
                       Cell = rep(levels(total_T_sub$cluster_name_slim),each=2),
                       row.names = colnames(matrix))
anno_color <- list(Cell = c(CD4Tmemory = hue_pal()(6)[1],
                            "CD4Texhausted" = hue_pal()(6)[2],
                            "CD4Tstress" = hue_pal()(6)[3],
                            "CD4Treg" = hue_pal()(6)[4],
                            "CD8Tmemory" = hue_pal()(6)[5],
                            "CD8Teffector" = hue_pal()(6)[6]),
                   
                   Group = c(CAS_NS = "#2b90d9", CAS_Smo = "#f26d5b"))

pheatmap::pheatmap(matrix,scale = "row",cluster_cols = F,
                   annotation_col = anno_col,
                   annotation_colors = anno_color,
                   display_numbers = test,
                   fontsize_number = 12,
                   show_colnames = F,
                   gaps_col = seq(2,10,2),
                   width = 12,
                   cellwidth = 16,cellheight = 12,
                   fontsize_row = 14,
                   treeheight_row = 20,
                   filename = "FigureS3_cytokine_heatmap_T.pdf")

######################Figure 6 CellChat##############################
library(CellChat)

CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")

#total
total <- readRDS("total_final.2.rds")

data.input <- as.matrix(total@assays$SCT@data)
meta <- total@meta.data

for(x in c("Ctl","CAS")){
  cell.use <- rownames(meta)[meta$disease == x]
  data.input_2 <- data.input[, cell.use]
  meta_2 = meta[cell.use, ]
  
  cellchat <- createCellChat(object = data.input_2, meta = meta_2, group.by = "cluster_name_slim")
  groupSize <- as.numeric(table(cellchat@idents))
  cellchat@DB <- CellChatDB.use
  
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)
  cellchat <- computeCommunProb(cellchat, raw.use = T,type = "truncatedMean",trim = 0.2,population.size = F)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  df.net <- subsetCommunication(cellchat)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  assign(paste0("cellchat_",x),cellchat)
}


object.list <- list(Ctl = cellchat_Ctl, CAS = cellchat_CAS)
saveRDS(object.list, file="cellchat.object.list.rds")

cellchat <- mergeCellChat(object.list, add.names = names(object.list))

weight.max <- getMaxWeight(object.list, attribute = c("idents","weight"))
pdf("Figure6A_InteractionWeights_ADSC.pdf",width = 10,height = 5)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$weight,idents.use = c("ADSC"), weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Interaction weights - ", names(object.list)[i]))
}
dev.off()

weight.max <- getMaxWeight(object.list, attribute = c("idents","weight"))
pdf("FigureS4_InteractionWeights.pdf",width = 10,height = 5)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$weight, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Interaction weights - ", names(object.list)[i]))
}
dev.off()

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(unlist(num.link)), max(unlist(num.link))) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]], slot.name = "netP")
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)+xlim(0,3.5)+ylim(0,2)
}
pdf("Figure6B_OutIncomeStrength.pdf",width = 10,height = 5)
patchwork::wrap_plots(plots = gg)
dev.off()

pdf("Figure6C_InteractionStrengthDiff.pdf",width = 5,height = 4)
netVisual_heatmap(cellchat, measure = "weight")
dev.off()

#ADSC and innate immune cell
total_immune <- readRDS("total_innate_immune.rds")
total_ADSC <- readRDS("total_ADSC.rds")
total <- merge(total_ADSC,total_immune)

data.input <- as.matrix(total@assays$SCT@data)
meta <- total@meta.data

for(x in c("Ctl","CAS")){
  cell.use <- rownames(meta)[meta$disease == x]
  data.input_2 <- data.input[, cell.use]
  meta_2 = meta[cell.use, ]
  
  cellchat <- createCellChat(object = data.input_2, meta = meta_2, group.by = "cluster_name_slim")
  groupSize <- as.numeric(table(cellchat@idents))
  cellchat@DB <- CellChatDB.use
  
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)
  cellchat <- computeCommunProb(cellchat, raw.use = T,type = "truncatedMean",trim = 0.2,population.size = F)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  df.net <- subsetCommunication(cellchat)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  assign(paste0("cellchat_",x),cellchat)
}


object.list <- list(Ctl = cellchat_Ctl, CAS = cellchat_CAS)
saveRDS(object.list, file="cellchat.object.list.innate.rds")

object.list <- readRDS("cellchat.object.list.innate.rds")
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

weight.max <- getMaxWeight(object.list, attribute = c("idents","weight"))
pdf("Figure6D_InteractionWeights_ADSC_innate.pdf",width = 10,height = 5)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$weight,idents.use = c("ADSC_CXCL14","ADSC_CD55"), weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Interaction weights - ", names(object.list)[i]))
}
dev.off()

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(unlist(num.link)), max(unlist(num.link))) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]], slot.name = "netP")
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)+xlim(0,6)+ylim(0,3)
}
pdf("Figure6E_OutIncomeStrength_innate.pdf",width = 10,height = 5)
patchwork::wrap_plots(plots = gg)
dev.off()

pdf("Figure6F_InteractionStrengthDiff_innate.pdf",width = 5,height = 4)
netVisual_heatmap(cellchat, measure = "weight")
dev.off()

p1 <- netVisual_bubble(cellchat,sources.use = c(1,2),targets.use = c(3:12),comparison = c(1, 2),max.dataset = 2, title.name = "Increased signaling in CAS",remove.isolate = T,color.text = c("#00BFC4","#F8766D"))
p2 <- netVisual_bubble(cellchat,sources.use = c(1,2),targets.use = c(3:12),comparison = c(1, 2),max.dataset = 1, title.name = "Decreased signaling in CAS",remove.isolate = T,color.text = c("#00BFC4","#F8766D"))
pdf("Figure6J_Bubble_ADSC_to_Myeloid.pdf",width = 15,height = 7)
p1+p2
dev.off()

p1 <- netVisual_bubble(cellchat,sources.use = c(3:12),targets.use = c(1,2),comparison = c(1, 2),max.dataset = 2, title.name = "Increased signaling in CAS",remove.isolate = T,color.text = c("#00BFC4","#F8766D"))
p2 <- netVisual_bubble(cellchat,sources.use = c(3:12),targets.use = c(1,2),comparison = c(1, 2),max.dataset = 1, title.name = "Decreased signaling in CAS",remove.isolate = T,color.text = c("#00BFC4","#F8766D"))
pdf("Figure6J_Bubble_Myeloid_to_ADSC.pdf",width = 15,height = 7)
p1+p2
dev.off()


df.net <- subsetCommunication(object.list[[2]])

pathways.show <- c("PTN") 
LR.show <- c("PTN_NCL")
pdf("Figure6K_Pathway_ChordPlot_innate_PTN.pdf",width = 5,height = 5)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_individual(object.list[[i]], signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
}
dev.off()

pathways.show <- c("ANNEXIN") 
LR.show <- c("ANXA1_FPR1")
pdf("Figure6K_Pathway_ChordPlot_innate_ANXA1.pdf",width = 5,height = 5)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_individual(object.list[[i]], signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
}
dev.off()

pathways.show <- c("COMPLEMENT") 
LR.show <- c("C3_ITGAX_ITGB2")
pdf("Figure6K_Pathway_ChordPlot_innate_C3.pdf",width = 5,height = 5)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_individual(object.list[[i]], signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
}
dev.off()

pathways.show <- c("MIF") 
LR.show <- c("MIF_CD74_CXCR4")
pdf("Figure6K_Pathway_ChordPlot_innate_MIF.pdf",width = 5,height = 5)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_individual(object.list[[i]], signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
}
dev.off()

pathways.show <- c("VEGF") 
LR.show <- c("VEGFA_VEGFR1")
pdf("Figure6L_Pathway_ChordPlot_innate_VEGFA.pdf",width = 5,height = 5)
netVisual_individual(object.list[[2]], signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
dev.off()

pathways.show <- c("VEGF") 
LR.show <- c("VEGFB_VEGFR1")
pdf("Figure6L_Pathway_ChordPlot_innate_VEGFB.pdf",width = 5,height = 5)
netVisual_individual(object.list[[2]], signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
dev.off()

list <- c("LGALS9_CD44","EREG_EGFR","AREG_EGFR","MIF_ACKR3","NAMPT_ITGA5_ITGB1","TNFSF12_TNFRSF12A",
          "GAS6_AXL","CXCL12_ACKR3","CXCL12_CXCR4","MIF_CD74_CXCR4","C3_ITGAX_ITGB2","ANXA1_FPR1",
          "PTN_NCL","ANGPTL1_PIRB","CSF1_CSF1R","FGF2_FGFR1","RARRES2_CMKLR1","SEMA3C_PLXND1")

df.net.sub <- df.net[df.net$interaction_name %in% list,c("interaction_name","pathway_name")] %>% distinct()

for(x in 1:nrow(df.net.sub)){
  pathways.show <- df.net.sub$pathway_name[x] 
  LR.show <- df.net.sub$interaction_name[x]
  pdf(paste0("Figure6_Pathway_ChordPlot_innate_",LR.show,".pdf"),width = 5,height = 5)
  par(mfrow = c(1,2), xpd=TRUE)
  for (i in 1:length(object.list)) {
    netVisual_individual(object.list[[i]], signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
  }
  dev.off()
}



#ADSC and T cell
total_T <- readRDS("total_T.rds")
total_ADSC <- readRDS("total_ADSC.rds")

total <- merge(total_ADSC,total_T)

data.input <- as.matrix(total@assays$SCT@data)
meta <- total@meta.data

for(x in c("Ctl","CAS")){
  cell.use <- rownames(meta)[meta$disease == x]
  data.input_2 <- data.input[, cell.use]
  meta_2 = meta[cell.use, ]
  
  cellchat <- createCellChat(object = data.input_2, meta = meta_2, group.by = "cluster_name_slim")
  groupSize <- as.numeric(table(cellchat@idents))
  cellchat@DB <- CellChatDB.use
  
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)
  cellchat <- computeCommunProb(cellchat, raw.use = T,type = "truncatedMean",trim = 0.2,population.size = F)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  df.net <- subsetCommunication(cellchat)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  assign(paste0("cellchat_",x),cellchat)
}


object.list <- list(Ctl = cellchat_Ctl, CAS = cellchat_CAS)
saveRDS(object.list, file="cellchat.object.list.T.rds")

object.list <- readRDS("cellchat.object.list.T.rds")
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

weight.max <- getMaxWeight(object.list, attribute = c("idents","weight"))
pdf("Figure6G_InteractionWeights_ADSC_T.pdf",width = 10,height = 5)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$weight,idents.use = c("ADSC_CXCL14","ADSC_CD55"), weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Interaction weights - ", names(object.list)[i]))
}
dev.off()

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(unlist(num.link)), max(unlist(num.link))) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]], slot.name = "netP")
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)+xlim(0,3)+ylim(0,2)
}
pdf("Figure6H_OutIncomeStrength_T.pdf",width = 10,height = 5)
patchwork::wrap_plots(plots = gg)
dev.off()

pdf("Figure6I_InteractionStrengthDiff_T.pdf",width = 5,height = 4)
netVisual_heatmap(cellchat, measure = "weight")
dev.off()

p1 <- netVisual_bubble(cellchat,sources.use = c(1,2),targets.use = c(3:12),comparison = c(1, 2),max.dataset = 2, title.name = "Increased signaling in CAS",remove.isolate = T,color.text = c("#00BFC4","#F8766D"))
p2 <- netVisual_bubble(cellchat,sources.use = c(1,2),targets.use = c(3:12),comparison = c(1, 2),max.dataset = 1, title.name = "Decreased signaling in CAS",remove.isolate = T,color.text = c("#00BFC4","#F8766D"))
pdf("Figure6M_Bubble_ADSC_to_T.pdf",width = 15,height = 3.5)
p1+p2
dev.off()

p1 <- netVisual_bubble(cellchat,sources.use = c(3:12),targets.use = c(1,2),comparison = c(1, 2),max.dataset = 2, title.name = "Increased signaling in CAS",remove.isolate = T,color.text = c("#00BFC4","#F8766D"))
p2 <- netVisual_bubble(cellchat,sources.use = c(3:12),targets.use = c(1,2),comparison = c(1, 2),max.dataset = 1, title.name = "Decreased signaling in CAS",remove.isolate = T,color.text = c("#00BFC4","#F8766D"))
pdf("Figure6M_Bubble_T_to_ADSC.pdf",width = 15,height = 6)
p1+p2
dev.off()

df.net <- subsetCommunication(object.list[[2]])

pathways.show <- c("CXCL") 
LR.show <- c("CXCL12_CXCR4")
pdf("Figure6N_Pathway_CirclePlot_T_CXCL.pdf",width = 5,height = 5)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_individual(object.list[[i]], signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
}
dev.off()

pathways.show <- c("PTN") 
LR.show <- c("PTN_NCL")
pdf("Figure6N_Pathway_CirclePlot_T_PTN.pdf",width = 5,height = 5)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_individual(object.list[[i]], signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
}
dev.off()

pathways.show <- c("MIF") 
LR.show <- c("MIF_CD74_CXCR4")
pdf("Figure6K_Pathway_CirclePlot_T_MIF.pdf",width = 5,height = 5)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_individual(object.list[[i]], signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
}
dev.off()

list <- c("MIF_ACKR3","MIF_CD74_CD44","IFNG_IFNGR1_IFNGR2","CD40LG_ITGA5_ITGB1",
          "AREG_EGFR","CD70_CD27","PTN_NCL","BTLA_TNFRSF14")

df.net.sub <- df.net[df.net$interaction_name %in% list,c("interaction_name","pathway_name")] %>% distinct()

for(x in c(2,3,4,6,7,8)){
  pathways.show <- df.net.sub$pathway_name[x] 
  LR.show <- df.net.sub$interaction_name[x]
  pdf(paste0("Figure6_Pathway_ChordPlot_T_",LR.show,".pdf"),width = 5,height = 5)
  par(mfrow = c(1,2), xpd=TRUE)
  for (i in 1:length(object.list)) {
    netVisual_individual(object.list[[i]], signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
  }
  dev.off()
}

x=1
pathways.show <- df.net.sub$pathway_name[x] 
LR.show <- df.net.sub$interaction_name[x]
pdf(paste0("Figure6_Pathway_ChordPlot_T_",LR.show,"_onlyCAS.pdf"),width = 5,height = 5)
i=2
netVisual_individual(object.list[[i]], signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
dev.off()

x=5
pathways.show <- df.net.sub$pathway_name[x] 
LR.show <- df.net.sub$interaction_name[x]
pdf(paste0("Figure6_Pathway_ChordPlot_T_",LR.show,"_onlyCAS.pdf"),width = 5,height = 5)
i=2
netVisual_individual(object.list[[i]], signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
dev.off()
######################Public ADSC data############################

raw <- Read10X(data.dir="./PublicData/Science_2019/data/")
data <- CreateSeuratObject(raw,project="Science",min.cells=10)
data <- RenameCells(object=data,add.cell.id="Science")
data$sample <- "Science"
data <- PercentageFeatureSet(data, pattern = "^MT-", col.name = "percent.mt")
data <- subset(data, subset = nFeature_RNA > 400 & nFeature_RNA < 8000 & percent.mt < 15)
data <- SCTransform(data)

data <- RunPCA(data, verbose = FALSE)
ElbowPlot(data, ndims = 50)
data <- RunUMAP(data, dims = 1:30)
data <- FindNeighbors(data, dims = 1:30)
for(i in c(0.5,0.8,1,1.5)){
  data <- FindClusters(data, resolution = i)
}

Idents(data) <- "SCT_snn_res.0.5"
DimPlot(data,label = T)

markers <- c("PTPRC","MS4A1","CD19","BANK1","JCHAIN","SDC1","IGKC","CD3D","CD4","CCR7","SELL","CD8A","GZMK","CCL5","KLRD1","GNLY","NKG7","MKI67","STMN1","TUBA1B","IL1B","S100A8","S100A9","CD68","C1QA","LYVE1","CD1E","LAMP3","FSCN1","PTGDS","IL3RA","LILRA4","F3","CXCL14","CD55","HBB","HBA1","VWF")
DotPlot(data, features = markers, dot.scale = 6)  + theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) + RotatedAxis()

data <- subset(data, subset = SCT_snn_res.0.5 != "11")

VlnPlot(data,features = c("DPP4","CD55","ICAM1",'CXCL14',"TAGLN"))

data$cluster_name_slim <- "ADSC_CXCL14"
data$cluster_name_slim[data$SCT_snn_res.0.5 == "2"] <- "ADSC_CD55"
data$cluster_name_slim[data$SCT_snn_res.0.5 == "6"] <- "Endothelial"
data$cluster_name_slim[data$SCT_snn_res.0.5 == "8"] <- "SmoothMuscle"

Idents(data) <- "cluster_name_slim"

saveRDS(data,file = "./PublicData/Science_2019/data.rds")

public_ADSC_2 <- subset(data, subset = cluster_name_slim %in% c("ADSC_CXCL14","ADSC_CD55"))
public_ADSC_2$disease <- "Adipose"

#Nature 2020 BAT and WAT

total <- readRDS("PublicData/total_BAT_WAT.rds")

DimPlot(total,label = T)

public_ADSC <- subset(total,subset = cluster_name_slim == "ADSC")
public_ADSC$disease <- public_ADSC$group

public_ADSC <- FindVariableFeatures(public_ADSC,selection.method = "vst", nfeatures = 5000)
public_ADSC <- ScaleData(public_ADSC)
public_ADSC <- RunPCA(public_ADSC, verbose = FALSE)

public_ADSC@meta.data$sample <- as.factor(public_ADSC@meta.data$sample)
public_ADSC <- RunHarmony(public_ADSC,group.by.vars="sample",plot_convergence = TRUE)
ElbowPlot(public_ADSC,ndims = 50)
public_ADSC <- RunUMAP(public_ADSC,reduction = "harmony", dims = 1:30)
public_ADSC <- FindNeighbors(public_ADSC,reduction = "harmony", dims = 1:30)

DimPlot(public_ADSC,split.by = "group")

FeaturePlot(public_ADSC,features = c("CXCL14","CD55"),order = T)

#整合效果不好

####################Integrate public and in-house data (ADSC)###########################

total_ADSC <- readRDS("total_ADSC.rds")

DimPlot(total_ADSC)

total_ADSC[["integrated"]] <- NULL

total_ADSC$cohort <- "InHouse"

data.list1 <- SplitObject(total_ADSC,split.by = "sample")

#public data

#ATVB_2024

data_AVTB2024 <- qread("../PublicSCData/ATVB_2024/total_ADSC.qs")
data_AVTB2024$cohort <- "AVTB2024"
data.list2 <- SplitObject(data_AVTB2024,split.by = "sample")

#Nature_2020

data_Nature2020 <- qread("../PublicSCData/Nature_2020/total_ADSC.qs")
data_Nature2020$cohort <- paste0("Nature2020_",data_Nature2020$group)
data.list3 <- SplitObject(data_Nature2020,split.by = "sample")

#NM_2020

data_NM2020 <- qread("../PublicSCData/NatureMeta_2020/total_ADSC.qs")
data_NM2020$SampleName <- droplevels(data_NM2020$SampleName)
data_NM2020$sample <- data_NM2020$SampleName
data_NM2020$cohort <- "NM2020"
data.list4 <- SplitObject(data_NM2020,split.by = "sample")

#NM_2021

data_NM2021 <- qread("../PublicSCData/NatureMeta_2021/total_ADSC.qs")
data_NM2021$cohort <- "NM2021"
data.list5 <- SplitObject(data_NM2021,split.by = "sample")

#Science_2019

data_Science2019 <- qread("../PublicSCData/Science_2019/total_ADSC.qs")
data_Science2019$cohort <- "Science2019"
data_Science2019$sample <- "Science2019"

data.list <- c(data.list1,data.list2,data.list3,data.list4,data.list5,data_Science2019)

keep_elements <- sapply(data.list, function(x) {
  # 检查元素是否有列数属性（如矩阵、数据框）
  cols <- tryCatch(ncol(x), error = function(e) NULL)
  # 如果列数存在且大于等于40，则保留（返回TRUE），否则不保留（返回FALSE）
  ifelse(!is.null(cols) && cols >= 40, TRUE, FALSE)
})

data.list <- data.list[keep_elements]

data.list <- lapply(data.list,FUN = function(x){DefaultAssay(x) <- "RNA";x <- SCTransform(x)})

options(future.globals.maxSize= 100000000000000)

data.features <- SelectIntegrationFeatures(object.list = data.list)

#getAnchor.features=function(obj){
#   obj <- GetResidual(object = obj,
#                      features = data.features,
#                      verbose = FALSE)
# 
#   scale.data <- GetAssayData(object = obj, assay = "SCT",
#                              slot = "scale.data")
# 
#   row.names(scale.data)
# }

#c=lapply(data.list,getAnchor.features)
#data.features=Reduce(intersect, c(c,list(data.features=data.features)))

data.list <- PrepSCTIntegration(object.list = data.list,assay = "SCT", anchor.features = data.features,verbose = FALSE)
data.anchors <- FindIntegrationAnchors(object.list = data.list, normalization.method = "SCT", anchor.features = data.features,k.filter = 40, verbose = FALSE)
total_ADSC <- IntegrateData(anchorset = data.anchors, normalization.method = "SCT",k.weight = 40, verbose = FALSE)

total_ADSC <- Reduct(total_ADSC)

Idents(total_ADSC) <- "integrated_snn_res.0.5"
DimPlot(total_ADSC,label = T)
DimPlot(total_ADSC,label = T,group.by = "cluster_name_slim")

DefaultAssay(total_ADSC) <- "SCT"
VlnPlot(total_ADSC,features = c("CXCL14","CD55"))

total_ADSC$cluster_name_slim <- "ADSC_CXCL14"
total_ADSC$cluster_name_slim[total_ADSC$integrated_snn_res.0.5 %in% c(1,2,6,7)] <- "ADSC_CD55"

total_ADSC$cluster_name_slim <- factor(total_ADSC$cluster_name_slim,
                                       levels = c("ADSC_CXCL14","ADSC_CD55"))

Idents(total_ADSC) <- "cluster_name_slim"

total_ADSC$cohort[total_ADSC$cohort == "NM2020"] <- "NM2020_SAT"
total_ADSC$cohort[total_ADSC$Tissue == "VAT"] <- "NM2020_VAT"

total_ADSC$cohort <- factor(total_ADSC$cohort,levels = c("InHouse","Nature2020_BAT","NM2021","AVTB2024",
                                                         "Science2019","NM2020_SAT","Nature2020_WAT",
                                                         "NM2020_VAT"))

pdf(file = "Dimplot_ADSC_mergePublic.pdf",width = 20,height = 4)
DimPlot(total_ADSC,reduction = "umap",label = F,pt.size = 0.15,cols = pal[c(1,3)],split.by = "cohort")+theme(legend.position = "bottom")
dev.off()

pdf(file = "DensityPlot_ADSC_mergePublic.pdf",width = 10,height = 5)
Plot_Density_Custom(seurat_object = total_ADSC, features = c("CXCL14","CD55"))
dev.off()

total_ADSC$group <- "PVAT"
total_ADSC$group[total_ADSC$cohort == "InHouse"] <- "InHouse"
total_ADSC$group[total_ADSC$cohort %in% c("Science2019","NM2020_SAT","Nature2020_WAT")] <- "SAT"
total_ADSC$group[total_ADSC$cohort %in% c("NM2020_VAT")] <- "VAT"

i=1
diff_gene <- list()
for(x in c("ADSC_CXCL14","ADSC_CD55")){
  total_ADSC_CD55 <- subset(total_ADSC,subset = cluster_name_slim == x)
  
  DefaultAssay(total_ADSC_CD55) <- "RNA"
  total_ADSC_CD55 <- NormalizeData(total_ADSC_CD55)
  
  data.features <- row.names(total_ADSC_CD55)
  RB_genes = c(grep("^RPS", data.features, v=T),grep("^RPL", data.features, v=T))
  MT_genes = grep("^MT-", data.features, v=T)
  remove <- c(RB_genes,MT_genes)
  total_ADSC_CD55_remove <- subset(total_ADSC_CD55,features = setdiff(data.features,remove))
  
  for(y in c("PVAT","SAT","VAT")){
    Idents(total_ADSC_CD55_remove) <- "group"
    marker <- FindMarkers(total_ADSC_CD55_remove, ident.1 = "InHouse", ident.2 = y,logfc.threshold = 0, min.pct = 0.3, min.diff.pct = 0.01,only.pos = F,recorrect_umi=FALSE)
    
    gene_up <- filter(marker,p_val_adj < 0.05,avg_log2FC > 1) %>% row.names()
    gene_down <- filter(marker,p_val_adj < 0.05,avg_log2FC < (-1)) %>% row.names()
    
    diff_gene[[i]] <- gene_up
    diff_gene[[i+1]] <- gene_down
    names(diff_gene)[i] <- paste(x,y,"up",sep = "_")
    names(diff_gene)[i+1] <- paste(x,y,"down",sep = "_")
    i=i+2
    
    exp <- AverageExpression(total_ADSC_CD55,assays = "RNA",group.by = "group") %>% as.data.frame()
    exp <- log2(exp+1)
    
    colnames(exp) <- c("Inhouse","PVAT","SAT","VAT")
    
    exp <- exp[,c("Inhouse",y)]
    colnames(exp)[2] <- "Ctrl"
    
    exp$Sig <- "NS"
    exp[gene_up,"Sig"] <- "Up"
    exp[gene_down,"Sig"] <- "Down"
    
    p <- ggplot(exp, aes(x = Inhouse, y = Ctrl, color = Sig)) +
      geom_point(aes(size=Sig),size=1) +
      scale_color_manual(values = c("Up" = alpha("red",0.8), "Down" = alpha("blue",0.8), "NS" = alpha("lightgrey",0.8))) + 
      labs(
        x = "Expression in in house data",
        y = paste0("Expression in ",y),
        color = "Significance"
      ) +
      theme_classic()+
      theme(
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.grid =element_blank(),
        axis.line = element_line(size=0),
        axis.text = element_text(color = "black",size = 12),
        axis.title = element_text(color = "black",size = 14),
        legend.position = "none",
        legend.key.size = unit(0.5, "cm"),        # 图例项大小
        legend.text = element_text(size = 12),    # 图例文字大小
        legend.title = element_text(size = 12),   # 图例标题大小
        legend.spacing = unit(0.1, "cm")
      ) +
      guides(color = guide_legend(override.aes = list(size = 3)));p
    
    ggsave(paste0("ExpVolcanoPlot_",x,"_",y,".pdf"),p,width = 4,height = 4)
    
  }
}

qsave(diff_gene,file = "diff_gene_mergePublic.qs")

##############################WebGseatlas#######################################

library(WebGestaltR)
library(Hmisc)

for(i in 1:length(diff_gene)){
  gene <- diff_gene[[i]]
  name <- names(diff_gene)[i]
  
  if (i %% 2 == 0) {color <- "blue"}else{color <- "red"}
  
  GOBP <- WebGestaltR(enrichMethod = "ORA", organism = "hsapiens",
                      enrichDatabase = "geneontology_Biological_Process_noRedundant",
                      enrichDatabaseType = "genesymbol",
                      interestGene = gene, interestGeneType = "genesymbol",
                      referenceSet = "genome",
                      fdrThr = 0.2,
                      isOutput = F)
  GOBP <- GOBP[GOBP$FDR<0.05&GOBP$overlap>3,]
  
  gobp_plot <- GOBP
  gobp_plot$FDR <- (-log10(gobp_plot$pValue))
  gobp_plot <- gobp_plot[order(gobp_plot$enrichmentRatio,decreasing = T),]
  gobp_plot <- gobp_plot[1:10,] %>% na.omit()
  gobp_plot$description <- capitalize(gobp_plot$description)
  gobp_plot$description <- factor(gobp_plot$description,levels = rev(gobp_plot$description))
  
  p <- ggplot(gobp_plot,aes(x=enrichmentRatio,y=description))+
    geom_bar(stat="identity",fill=alpha(color,0.5))+
    geom_text(aes(x=0.1,y=description,label=description),hjust=0)+
    #scale_fill_gradientn(colours = rev(brewer.pal(9,"Reds")))+
    theme_classic()+
    scale_x_continuous(expand = c(0,0))+
    theme(legend.position = "right",
          legend.text = element_text(color = "black",size = 14),
          legend.title = element_blank(),
          axis.text = element_text(color = "black",size = 14),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(color = "black",size = 16),
          axis.title.y = element_blank())+
    xlab("Enrichment Ratio");p
  
  ggsave(paste0("GeneEnrichment_",name,".pdf"),p,width = 5,height = 4)
  
}


##############################Cell Proportion of ADSC#######################################

cell.prop <- table(total_ADSC$cluster_name_slim,total_ADSC$cohort) %>% as.data.frame.matrix()
#cell.prop <- t(t(cell.prop)/colSums(cell.prop)) %>% as.data.frame()
cell.prop$cell <- row.names(cell.prop)
cell.prop <- pivot_longer(cell.prop,cols = setdiff(colnames(cell.prop),"cell"),names_to = "cohort",values_to = "num")
cell.prop = ddply(cell.prop, 'cohort', transform, percent=num/sum(num)*100)

cell.prop$cell <- factor(cell.prop$cell,levels = rev(levels(total_ADSC$cluster_name_slim)))
cell.prop$cohort <- factor(cell.prop$cohort,levels = levels(total_ADSC$cohort))

p <- ggplot(cell.prop,aes(x=cohort,y=percent,fill=cell))+
  geom_bar(stat = 'identity',width = 0.8,colour='black')+
  scale_fill_manual(values = pal[c(3,1)])+
  scale_y_continuous(expand = c(0,0))+
  theme_classic()+
  theme(axis.text = element_text(size = 14,colour = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 16,colour = "black"))+
  ylab("Cell ratio (%)")+xlab("");p
ggsave("Cell_proportion_bar_ADSC_mergePublic.pdf",p,width = 6,height = 3.5)


##############################Pseudo-Bulk#######################################
library(SingleCellExperiment)
library(scuttle)
library(Matrix.utils)
library(IOBR)
library(DESeq2)
library(sva)
library(UpSetR)

#pseudo-bulk

pseudoBulk <- function(total){
  count <- total@assays$RNA@counts %>% as.data.frame()
  pesudo_metadata <- total@meta.data
  pesudo_metadata$sample_id <- factor(total$sample)
  
  sce_bulk <- SingleCellExperiment(assay = list(counts = count),
                                   colData = pesudo_metadata)
  
  group <- colData(sce_bulk)[,c("sample_id")]
  pb <- aggregate.Matrix(t(counts(sce_bulk)), 
                         groupings = group, fun = "sum")
  
  pb_t <- t(pb) %>% as.matrix()
  
  return(pb_t)
}

anno <- read.table('annotation.txt')
anno <- filter(anno,V1 == "protein_coding")
proCod <- anno$V2

total_ADSC_CD55 <- subset(total_ADSC,subset = cluster_name_slim == "ADSC_CD55")

meta <- total_ADSC_CD55@meta.data[,c("sample","cohort","group")] %>% distinct()
row.names(meta) <- meta$sample

pb2 <- pseudoBulk(total_ADSC_CD55)

qsave(pb2,file = "PseudoBulk_ADSC_CD55.qs")

for(x in c("PVAT","SAT","VAT")){
  
  meta_sub <- filter(meta,group %in% c("InHouse",x))
  pb2_sub <- pb2[,row.names(meta_sub)]
  
  meta_sub$group <- factor(meta_sub$group,levels = c("InHouse",x))
  
  pb2_sub <- pb2_sub[row.names(pb2_sub) %in% proCod,]
  
  data.features <- row.names(pb2_sub)
  RB_genes = c(grep("^RPS", data.features, v=T),grep("^RPL", data.features, v=T))
  MT_genes = grep("^MT-", data.features, v=T)
  remove <- c(RB_genes,MT_genes)
  pb2_sub <- pb2_sub[!(row.names(pb2_sub) %in% remove),]
  
  dds <- DESeqDataSetFromMatrix(countData=pb2_sub,colData=meta_sub,design=~group)
  dds <- DESeq(dds)
  
  DEG <- results(dds) %>% as.data.frame() %>% na.omit()
  
  write.csv(DEG,file = paste0("DEG_ADSC_CD55_",x,".csv"),row.names = T,col.names = T)
  
}


#Volcano plot

library(ggrepel)

for(x in c("PVAT","SAT","VAT")){
  deg <- read.csv(paste0("DEG_ADSC_CD55_",x,".csv"),row.names = 1,header = T)
  
  data_vol <- subset(deg,select=c(6,2))
  colnames(data_vol) <- c("pvalue","log2FoldChange")
  data_vol <- na.omit(data_vol)
  data_vol$sig[((data_vol$log2FoldChange<2)&(data_vol$log2FoldChange>-2))|(data_vol$pvalue>0.05)] <- "no"
  data_vol$sig[(data_vol$log2FoldChange>=2)&(data_vol$pvalue<=0.05)] <- "up"
  data_vol$sig[(data_vol$log2FoldChange<=-2)&(data_vol$pvalue<=0.05)] <- "down"
  data_vol$label <- row.names(data_vol)
  data_vol$label[data_vol$sig=="no"] <- ""
  
  data_vol <- data_vol[order(data_vol$pvalue,decreasing = F),]
  
  top5 <- data_vol[data_vol$log2FoldChange>0,] %>% .[1:5,] %>% row.names()
  bottom5 <- data_vol[data_vol$log2FoldChange<0,] %>% .[1:5,] %>% row.names()
  
  select <- c(top5,bottom5)
  
  data_vol$label[!row.names(data_vol) %in% select] <- ""
  
  cols <- c("steelblue","grey","firebrick")
  
  names(cols) <- c("down","no","up")
  
  p <- ggplot(data_vol,aes(log2FoldChange,-1*log10(pvalue),
                           color = sig))+geom_point()+
    geom_text_repel(aes(log2FoldChange,-1*log10(pvalue),label=label),size=3)+
    labs(x="log2(FoldChange)",y="-log10(FDR)") + 
    scale_color_manual(values = cols)+
    geom_hline(yintercept=-log10(0.05),linetype=4)+
    geom_vline(xintercept=c(-2,2),linetype=4)+
    theme_bw()+
    ylim(0,100)+
    theme(panel.grid.major =element_blank(),
          panel.grid.minor =element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color="black"),
          strip.background.x = element_rect(colour = "white"),
          axis.text=element_text(size=20),
          axis.title=element_text(size=20),
          legend.position = "none")
  
  ggsave(filename = paste0("./Volcano_ADSC_CD55_",x,".pdf"),p,device = "pdf",width = 6,height = 6)
}



##############################DEG of ADSC_CD55#######################################

total_CD55 <- subset(total_ADSC,subset = cluster_name_slim == "ADSC_CD55")

Idents(total_CD55) <- "disease"
marker <- FindMarkers(total_CD55, ident.1 = "CAS", ident.2 = "Adipose",logfc.threshold = 0, min.pct = 0.3, min.diff.pct = 0.01,only.pos = F,recorrect_umi=FALSE)
write.csv(marker,file = "CD55_CASvsAdiposeTissue.csv")

marker <- FindMarkers(total_CD55, ident.1 = "CAS", ident.2 = "BAT",logfc.threshold = 0, min.pct = 0.3, min.diff.pct = 0.01,only.pos = F,recorrect_umi=FALSE)
write.csv(marker,file = "CD55_CASvsDeepNeck.csv")

marker <- read.csv("CD55_CASvsDeepNeck.csv",header = T,row.names = 1)

##############################GSEA#######################################
library(clusterProfiler)
library(msigdbr)
library(gggsea)
library(ggsci)
library(enrichplot)

m_df_gobp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
m_df_gobp <- m_df_gobp %>% dplyr::select(gs_name,gene_symbol)

marker <- marker[order(marker$avg_log2FC,decreasing = T),]
genelist <- marker$avg_log2FC
names(genelist) <- rownames(marker)

gsea.re1 <- GSEA(genelist,  #待富集的基因列表
                TERM2GENE = m_df_gobp,  #基因集
                pvalueCutoff = 1,  #指定 p 值阈值（可指定 1 以输出全部）
                pAdjustMethod = 'BH')

g1 <- as.data.frame(gsea.re1)


col_gsea1<-pal_simpsons()(16)

num1=1
p <- gseaplot2(gsea.re1,geneSetID = "GOBP_VASCULATURE_DEVELOPMENT",
          title = "GOBP_VASCULATURE_DEVELOPMENT",#标题
          color = col_gsea1[1:num1],#颜色
          base_size = 14,#基础大小
          rel_heights = c(1, 0.2, 0.4),#小图相对高度
          subplots = 1:3,#展示小图
          pvalue_table = FALSE,#p值表格
          ES_geom = "line"#line or dot
)
ggsave("GSEA_DeepNeck.pdf",p,width = 6,height = 6)

gene <- g1["GOBP_VASCULATURE_DEVELOPMENT","core_enrichment"] %>% strsplit(.,"\\/")
gene <- gene[[1]]

gene <- marker[gene,]

gene <- intersect(row.names(gene_adipose),row.names(gene_deepNeck))

total_CD55_sub <- subset(total_CD55, subset = disease %in% c("BAT","Ctl","CAS"))

total_CD55_sub$disease <- factor(total_CD55_sub$disease,
                                 levels = c("BAT","Ctl","CAS"))

Idents(total_CD55_sub) <- "disease"

p <- DoHeatmap(total_CD55_sub,assay = "SCT", slot="data",
          features = gene,size = 3,label = T) + 
  scale_fill_gradientn(colors = c("navy","white","firebrick"))+
  theme(legend.title = element_text(size = 15,family = "Arial"),
        legend.text = element_text(size = 12,family = "Arial"),
        axis.text = element_text(family = "Arial"),
        strip.text = element_text(family = "Arial"));p
ggsave("Heatmap_DeepNeck.png",p,dpi = "retina",width = 6,height = 10)

######################Merge pla data######################

x <- "Patient 1"
count1 <- Read10X(data.dir=paste0(x,"/C230517003-pla"),gene.column=1) 
count1 <- CreateSeuratObject(counts = as.sparse(count1),project=x)
count2 <- Read10X(data.dir=paste0(x,"/C230517004-pla"),gene.column=1)
count2 <- CreateSeuratObject(counts = as.sparse(count2),project=x)
count3 <- Read10X(data.dir=paste0(x,"/C230517011-pla"),gene.column=1)
count3 <- CreateSeuratObject(counts = as.sparse(count3),project=x)
count4 <- Read10X(data.dir=paste0(x,"/C230517012-pla"),gene.column=1)
count4 <- CreateSeuratObject(counts = as.sparse(count4),project=x)

count <- merge(count1,list(count2,count3,count4))
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
data <- RenameCells(object=data,add.cell.id="Patient1")
data$sample <- "Patient1"
data <- PercentageFeatureSet(data, pattern = "^MT-", col.name = "percent.mt")
data <- subset(data, subset = nFeature_RNA > 400 & nFeature_RNA < 8000 & percent.mt < 15)
data <- SCTransform(data)

saveRDS(data,file = "Patient1.rds")

x <- "Patient 2"
count1 <- Read10X(data.dir=paste0(x,"/C231227001-pla"),gene.column=1) 
count1 <- CreateSeuratObject(counts = as.sparse(count1),project=x)
count2 <- Read10X(data.dir=paste0(x,"/C231227002-pla"),gene.column=1)
count2 <- CreateSeuratObject(counts = as.sparse(count2),project=x)
count3 <- Read10X(data.dir=paste0(x,"/C231227007-pla"),gene.column=1)
count3 <- CreateSeuratObject(counts = as.sparse(count3),project=x)
count4 <- Read10X(data.dir=paste0(x,"/C231227008-pla"),gene.column=1)
count4 <- CreateSeuratObject(counts = as.sparse(count4),project=x)

count <- merge(count1,list(count2,count3,count4))
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
data <- RenameCells(object=data,add.cell.id="Patient2")
data$sample <- "Patient2"
data <- PercentageFeatureSet(data, pattern = "^MT-", col.name = "percent.mt")
data <- subset(data, subset = nFeature_RNA > 400 & nFeature_RNA < 8000 & percent.mt < 15)
data <- SCTransform(data)

saveRDS(data,file = "Patient2.rds")

x <- "Patient 3"
count1 <- Read10X(data.dir=paste0(x,"/C231122001-pla"),gene.column=1) 
count1 <- CreateSeuratObject(counts = as.sparse(count1),project=x)
count2 <- Read10X(data.dir=paste0(x,"/C231122002-pla"),gene.column=1)
count2 <- CreateSeuratObject(counts = as.sparse(count2),project=x)


count <- merge(count1,count2)
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
data <- RenameCells(object=data,add.cell.id="Patient3")
data$sample <- "Patient3"
data <- PercentageFeatureSet(data, pattern = "^MT-", col.name = "percent.mt")
data <- subset(data, subset = nFeature_RNA > 400 & nFeature_RNA < 8000 & percent.mt < 15)
data <- SCTransform(data)

saveRDS(data,file = "Patient3.rds")

data.list <- list()
data.list[[1]] <- readRDS("Patient1.rds")
data.list[[2]] <- readRDS("Patient2.rds")
data.list[[3]] <- readRDS("Patient3.rds")

total2 <- readRDS("total_final.2.rds")

options(future.globals.maxSize= 100000000000000)

data.features <- SelectIntegrationFeatures(object.list = data.list)
data.list <- PrepSCTIntegration(object.list = data.list,assay = "SCT", anchor.features = data.features,verbose = FALSE)
data.anchors <- FindIntegrationAnchors(object.list = data.list, normalization.method = "SCT", anchor.features = data.features, verbose = FALSE)
total <- IntegrateData(anchorset = data.anchors, normalization.method = "SCT", verbose = FALSE)

total <- Reduct(total)

DimPlot(total, reduction = "umap", label=T,group.by = "integrated_snn_res.1")

DefaultAssay(total) <- "SCT"
Idents(total) <- "integrated_snn_res.1"
markers <- c("PTPRC","MS4A1","CD19","BANK1","JCHAIN","SDC1","IGKC","CD3D","CD4","CCR7","SELL","CD8A","GZMK","CCL5","KLRD1","GNLY","NKG7","MKI67","STMN1","TUBA1B","IL1B","S100A8","S100A9","CD68","C1QA","LYVE1","CD1E","LAMP3","FSCN1","PTGDS","IL3RA","LILRA4","F3","ITGB1","CD34","COL1A1","CXCL14","CD55","HBB","HBA1","VWF","PECAM1")
DotPlot(total, features = markers, dot.scale = 6)  + theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) + RotatedAxis()

total <- PrepSCTFindMarkers(total)
clu6.markers <- FindMarkers(total,ident.1 = "Pericyte",min.pct = 0.5, min.diff.pct = 0.2, logfc.threshold =0.25, only.pos = TRUE, test.use = "wilcox")

total$cluster_name_slim <- "CD4T cell"
total$cluster_name_slim[total$integrated_snn_res.0.5 %in% c(12)] <- "B cell"
total$cluster_name_slim[total$integrated_snn_res.0.5 %in% c(18)] <- "Plasma cell"
total$cluster_name_slim[total$integrated_snn_res.0.5 %in% c(2,3)] <- "CD8T cell"
total$cluster_name_slim[total$integrated_snn_res.1 %in% c(12)] <- "NK"
total$cluster_name_slim[total$integrated_snn_res.0.5 %in% c(15)] <- "Proliferating Lymphocyte"
total$cluster_name_slim[total$integrated_snn_res.0.5 %in% c(5)] <- "Monocyte/Neutrophil"
total$cluster_name_slim[total$integrated_snn_res.0.5 %in% c(4,7)] <- "Macrophage"
total$cluster_name_slim[total$integrated_snn_res.0.5 %in% c(9)] <- "DC"
total$cluster_name_slim[total$integrated_snn_res.0.5 %in% c(17)] <- "pDC"
total$cluster_name_slim[total$integrated_snn_res.0.5 %in% c(0,11)] <- "Fibroblast"
total$cluster_name_slim[total$integrated_snn_res.0.5 %in% c(6)] <- "Endothelial cell"
total$cluster_name_slim[total$integrated_snn_res.0.5 %in% c(10)] <- "Mast cell"
total$cluster_name_slim[total$integrated_snn_res.0.5 %in% c(8,13,14,16)] <- "Doublets"

total <- subset(total,subset = cluster_name_slim != "Doublets")

total$cluster_name_slim <- factor(total$cluster_name_slim,
                                  levels = c("B cell","Plasma cell","CD4T cell","CD8T cell","NK","Proliferating Lymphocyte","Monocyte/Neutrophil","Macrophage","DC","pDC","Mast cell","Fibroblast","Endothelial cell"))

DefaultAssay(total) <- "integrated"
total <- RunUMAP(total, dims = 1:30)

DimPlot(total,group.by = "cluster_name_slim",label=T)

qsave(total,file = "total_pla.qs")

Idents(total) <- "cluster_name_slim"

pal = DiscretePalette_scCustomize(num_colors = 14, palette = "ditto_seq")[-c(8)]

pdf(file = "Dimplot_mergePLA.pdf",width = 8,height = 6)
DimPlot_scCustom(total,colors_use = pal[1:13],figure_plot = T,label = T,label.size = 5,pt.size = 0.5)
dev.off()

markers <- c("PTPRC","MS4A1","CD19","BANK1","JCHAIN","SDC1","IGKC","CD3D","CD4","CCR7","SELL","CD8A","GZMK","CCL5","KLRD1","GNLY","NKG7","MKI67","STMN1","TUBA1B","IL1B","S100A8","S100A9","CD68","C1QA","LYVE1","CD1E","LAMP3","FSCN1","PTGDS","IL3RA","LILRA4","TPSB2","TPSAB1","CPA3","ITGB1","COL1A1","FN1","CD34","PECAM1","VWF")
Idents(total) <- "cluster_name_slim"
DefaultAssay(total) <- "SCT"
pdf(file = "DotPlot_mergePLA.pdf",width = 14,height = 5)
DotPlot(total, features = markers, dot.scale = 6)  + theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) + RotatedAxis()
dev.off()

####################Integrate public and in-house data (T cell)###########################

total_T <- readRDS("total_T.rds")

DimPlot(total_T)

DefaultAssay(total_T) <- "RNA"
total_T[["integrated"]] <- NULL

total_T$cohort <- "InHouse"

data.list1 <- SplitObject(total_T,split.by = "sample")

#public data

#ATVB_2024

data_AVTB2024 <- qread("../PublicSCData/ATVB_2024/total_T.qs")
data_AVTB2024$cohort <- "AVTB2024"
data.list2 <- SplitObject(data_AVTB2024,split.by = "sample")

#Nature_2020

data_Nature2020 <- qread("../PublicSCData/Nature_2020/total_T.qs")
data_Nature2020$cohort <- paste0("Nature2020_",data_Nature2020$group)
data.list3 <- SplitObject(data_Nature2020,split.by = "sample")

#NM_2020

data_NM2020 <- qread("../PublicSCData/NatureMeta_2020/total_T.qs")
data_NM2020$SampleName <- droplevels(data_NM2020$SampleName)
data_NM2020$sample <- data_NM2020$SampleName
data_NM2020$cohort <- paste0("NM2020_",data_NM2020$Tissue)
data.list4 <- SplitObject(data_NM2020,split.by = "sample")

#NM_2021

data_NM2021 <- qread("../PublicSCData/NatureMeta_2021/total_T.qs")
data_NM2021$cohort <- "NM2021"
data.list5 <- SplitObject(data_NM2021,split.by = "sample")


data.list <- c(data.list1,data.list2,data.list3,data.list4,data.list5)

keep_elements <- sapply(data.list, function(x) {
  # 检查元素是否有列数属性（如矩阵、数据框）
  cols <- tryCatch(ncol(x), error = function(e) NULL)
  # 如果列数存在且大于等于40，则保留（返回TRUE），否则不保留（返回FALSE）
  ifelse(!is.null(cols) && cols >= 40, TRUE, FALSE)
})

data.list <- data.list[keep_elements]

rm(list = setdiff(ls(),"data.list"))

data.list <- lapply(data.list,FUN = function(x){DefaultAssay(x) <- "RNA";x <- SCTransform(x)})

options(future.globals.maxSize= 100000000000000)

data.features <- SelectIntegrationFeatures(object.list = data.list)
data.list <- PrepSCTIntegration(object.list = data.list,assay = "SCT", anchor.features = data.features,verbose = FALSE)
data.anchors <- FindIntegrationAnchors(object.list = data.list, normalization.method = "SCT", anchor.features = data.features,k.filter = 40, verbose = FALSE)
total_T <- IntegrateData(anchorset = data.anchors, normalization.method = "SCT",k.weight = 40, verbose = FALSE)

total_T <- Reduct(total_T)

Idents(total_T) <- "integrated_snn_res.0.5"
DimPlot(total_T,label = T)
DimPlot(total_T,label = T,group.by = "cluster_name_slim")

qsave(total_T,file = "total_T_mergePublic.qs")

#harmony

DefaultAssay(total_T) <- "RNA"

total_T <- NormalizeData(total_T)
total_T <- FindVariableFeatures(total_T,selection.method = "vst", nfeatures = 1000)
total_T <- ScaleData(total_T)
total_T <- RunPCA(total_T, verbose = FALSE)

total_T@meta.data$sample <- as.factor(total_T@meta.data$sample)
total_T <- RunHarmony(total_T,group.by.vars="sample",plot_convergence = TRUE)
ElbowPlot(total_T,ndims = 50)
total_T <- RunUMAP(total_T,reduction = "harmony", dims = 1:30)
total_T <- FindNeighbors(total_T,reduction = "harmony", dims = 1:30)
total_T <- FindClusters(total_T, resolution = 0.5)

DimPlot(total_T, reduction = "umap", label=T,group.by = "RNA_snn_res.0.5")
DimPlot(total_T,label = T,group.by = "cluster_name_slim")

#FastMNN
library(SeuratWrappers)

for(i in 1:length(data.list)){
  DefaultAssay(data.list[[i]]) <- "RNA"
  data.list[[i]][["SCT"]] <- NULL
  data.list[[i]] <- NormalizeData(data.list[[i]])
  data.list[[i]] <- FindVariableFeatures(data.list[[i]])
}

total_T <- RunFastMNN(object.list = data.list)
total_T <- RunUMAP(total_T, reduction = "mnn", dims = 1:30)
total_T <- FindNeighbors(total_T, reduction = "mnn", dims = 1:30) 
total_T <- FindClusters(total_T,resolution = 0.5)

DimPlot(total_T, reduction = "umap", label=T,group.by = "RNA_snn_res.0.5")
DimPlot(total_T,label = T,group.by = "cluster_name_slim")

clu13.marker <- FindMarkers(total_T, ident.1 = "13", logfc.threshold = 0.2, min.pct = 0.5, min.diff.pct = 0.2,only.pos = T)

Plot_Density_Custom(total_T,features = c("PDCD1","HSPA1A"))
Plot_Density_Custom(total_T,features = c("NCAM1","FCGR3A"))

total_T$cluster_name_slim <- "CD4/CD8Tmemory"
total_T$cluster_name_slim[total_T$RNA_snn_res.0.5 %in% c(5,9)] <- "CD4Treg"
total_T$cluster_name_slim[total_T$RNA_snn_res.0.5 == 4] <- "CD8Teffector"
total_T$cluster_name_slim[total_T$RNA_snn_res.0.5 %in% c(3,6)] <- "Remove"
total_T$cluster_name_slim[total_T$RNA_snn_res.0.5 %in% c(13)] <- "NK"

total_NK <- subset(total_T,subset = cluster_name_slim == "NK")

total_T <- subset(total_T,subset = cluster_name_slim %in% c("CD4/CD8Tmemory","CD4Treg","CD8Teffector"))

total_T$cluster_name_slim <- factor(total_T$cluster_name_slim, levels = c("CD4/CD8Tmemory","CD4Treg","CD8Teffector"))

total_T <- subset(total_T,subset = cohort %in% c("InHouse","Nature2020_BAT","AVTB2024","NM2020_SAT"))

total_T$cohort <- factor(total_T$cohort,levels = c("InHouse","Nature2020_BAT","AVTB2024",
                                                  "NM2020_SAT"))

Idents(total_T) <- "cluster_name_slim"
pdf(file = "Dimplot_T_mergePublic.pdf",width = 15,height = 5)
DimPlot(total_T,reduction = "umap",label = F,pt.size = 0.15,split.by = "cohort")+theme(legend.position = "bottom")
dev.off()

#cell proportion

cell.prop <- table(total_T$cluster_name_slim,total_T$cohort) %>% as.data.frame.matrix()
#cell.prop <- t(t(cell.prop)/colSums(cell.prop)) %>% as.data.frame()
cell.prop$cell <- row.names(cell.prop)
cell.prop <- pivot_longer(cell.prop,cols = setdiff(colnames(cell.prop),"cell"),names_to = "cohort",values_to = "num")
cell.prop = ddply(cell.prop, 'cohort', transform, percent=num/sum(num)*100)

cell.prop$cell <- factor(cell.prop$cell,levels = rev(levels(total_T$cluster_name_slim)))
cell.prop$cohort <- factor(cell.prop$cohort,levels = levels(total_T$cohort))

p <- ggplot(cell.prop,aes(x=cohort,y=percent,fill=cell))+
  geom_bar(stat = 'identity',width = 0.8,colour='black')+
  scale_fill_manual(values = rev(hue_pal()(3)))+
  scale_y_continuous(expand = c(0,0))+
  theme_classic()+
  theme(axis.text = element_text(size = 14,colour = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 16,colour = "black"))+
  ylab("Cell ratio (%)")+xlab("");p
ggsave("Cell_proportion_bar_T_mergePublic.pdf",p,width = 5,height = 3.5)

qsave(total_T,file = "total_T_mergePublic.qs")

qsave(total_NK,file = "total_NK_Public.qs")

####################Integrate public and in-house data (Myeloid cell)###########################

total_Myeloid <- readRDS("total_myeloid.rds")

Idents(total_Myeloid) <- "cluster_name_slim"

DimPlot(total_Myeloid)

DefaultAssay(total_Myeloid) <- "RNA"
total_Myeloid[["integrated"]] <- NULL

total_Myeloid$cohort <- "InHouse"

data.list1 <- SplitObject(total_Myeloid,split.by = "sample")

#public data

#ATVB_2024

data_AVTB2024 <- qread("./PublicSCData/ATVB_2024/total_Mye.qs")
data_AVTB2024$cohort <- "AVTB2024"
data.list2 <- SplitObject(data_AVTB2024,split.by = "sample")

#Nature_2020

data_Nature2020 <- qread("./PublicSCData/Nature_2020/total_Mye.qs")
data_Nature2020$cohort <- paste0("Nature2020_",data_Nature2020$group)
data.list3 <- SplitObject(data_Nature2020,split.by = "sample")

#NM_2020

data_NM2020 <- qread("./PublicSCData/NatureMeta_2020/total_Mye.qs")
data_NM2020$SampleName <- droplevels(data_NM2020$SampleName)
data_NM2020$sample <- data_NM2020$SampleName
data_NM2020$cohort <- paste0("NM2020_",data_NM2020$Tissue)
data.list4 <- SplitObject(data_NM2020,split.by = "sample")

#NM_2021

data_NM2021 <- qread("./PublicSCData/NatureMeta_2021/total_Mye.qs")
data_NM2021$cohort <- "NM2021"
data.list5 <- SplitObject(data_NM2021,split.by = "sample")


data.list <- c(data.list1,data.list2,data.list3,data.list4,data.list5)

keep_elements <- sapply(data.list, function(x) {
  # 检查元素是否有列数属性（如矩阵、数据框）
  cols <- tryCatch(ncol(x), error = function(e) NULL)
  # 如果列数存在且大于等于40，则保留（返回TRUE），否则不保留（返回FALSE）
  ifelse(!is.null(cols) && cols >= 40, TRUE, FALSE)
})

data.list <- data.list[keep_elements]

rm(list = setdiff(ls(),"data.list"))

data.list <- lapply(data.list,FUN = function(x){DefaultAssay(x) <- "RNA";x <- SCTransform(x)})

options(future.globals.maxSize= 100000000000000)

data.features <- SelectIntegrationFeatures(object.list = data.list)
data.list <- PrepSCTIntegration(object.list = data.list,assay = "SCT", anchor.features = data.features,verbose = FALSE)
data.anchors <- FindIntegrationAnchors(object.list = data.list, normalization.method = "SCT", anchor.features = data.features,k.filter = 40, verbose = FALSE)
total_Myeloid <- IntegrateData(anchorset = data.anchors, normalization.method = "SCT",k.weight = 40, verbose = FALSE)

total_Myeloid <- Reduct(total_Myeloid)

Idents(total_Myeloid) <- "integrated_snn_res.0.5"
DimPlot(total_Myeloid,label = T)
DimPlot(total_Myeloid,label = T,group.by = "cluster_name_slim")

markers <- c("NCAM1",'FCGR3A',"CD69","NKG7","XCR1","CLEC9A","CD1C","CLEC10A","LAMP3","CCL19","CDKN1C",'LILRB2',"FCN1","S100A8","S100A9","CSF3R","CD163","C1QA","MRC1","IL3RA","LILRA4","TCF4")
Idents(total_Myeloid) <- "integrated_snn_res.0.5"
DefaultAssay(total_Myeloid) <- "SCT"
DotPlot(total_Myeloid, features = markers)  + theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) + RotatedAxis()

#FastMNN
library(SeuratWrappers)

count <- total_Myeloid@assays$RNA@counts
meta <- total_Myeloid@meta.data
total_Myeloid <- CreateSeuratObject(counts = count,meta.data = meta)

data.list <- SplitObject(total_Myeloid,split.by = "sample")

for(i in 1:length(data.list)){
  DefaultAssay(data.list[[i]]) <- "RNA"
  data.list[[i]] <- NormalizeData(data.list[[i]])
  data.list[[i]] <- FindVariableFeatures(data.list[[i]])
}

total_Myeloid <- RunFastMNN(object.list = data.list)
total_Myeloid <- RunUMAP(total_Myeloid, reduction = "mnn", dims = 1:30)
total_Myeloid <- FindNeighbors(total_Myeloid, reduction = "mnn", dims = 1:30) 
total_Myeloid <- FindClusters(total_Myeloid,resolution = 0.5)

DimPlot(total_Myeloid, reduction = "umap", label=T,group.by = "RNA_snn_res.0.5")
DimPlot(total_Myeloid,label = T,group.by = "cluster_name_slim")

markers <- c("PTPRC","NCAM1","FCGR3A","ITGA1","B3GAT1","KLRD1","GNLY","NKG7","IL1B","S100A8","S100A9","CSF3R","CD68","C1QA","LYVE1","XCR1","CLEC9A","CD1C","CD1E","LAMP3","FSCN1","PTGDS","IL3RA","LILRA4")
Idents(total_Myeloid) <- "RNA_snn_res.0.5"
DotPlot(total_Myeloid, features = markers, dot.scale = 6)  + theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) + RotatedAxis()

clu12.marker <- FindMarkers(total_Myeloid, ident.1 = "12", logfc.threshold = 0.2, min.pct = 0.5, min.diff.pct = 0.2,only.pos = T)

total_Myeloid$cluster_name_slim <- "Macrophage"
total_Myeloid$cluster_name_slim[total_Myeloid$RNA_snn_res.0.5 %in% c(4)] <- "Neutrophil/Monocyte"
total_Myeloid$cluster_name_slim[total_Myeloid$RNA_snn_res.0.5 == 8] <- "cDC2"
total_Myeloid$cluster_name_slim[total_Myeloid$RNA_snn_res.0.5 %in% c(10)] <- "pDC"
total_Myeloid$cluster_name_slim[total_Myeloid$RNA_snn_res.0.5 %in% c(11)] <- "cDC1"
total_Myeloid$cluster_name_slim[total_Myeloid$RNA_snn_res.0.5 %in% c(13)] <- "DC_LAMP3+"
total_Myeloid$cluster_name_slim[total_Myeloid$RNA_snn_res.0.5 %in% c(7,12)] <- "Remove"

total_Myeloid <- subset(total_Myeloid,subset = cluster_name_slim != "Remove")

total_Myeloid$cluster_name_slim <- factor(total_Myeloid$cluster_name_slim,
                                          levels = c("Macrophage","Neutrophil/Monocyte","cDC1","cDC2","DC_LAMP3+","pDC"))


total_Myeloid$cohort <- factor(total_Myeloid$cohort,levels = c("InHouse","Nature2020_BAT","NM2021","AVTB2024",
                                                         "NM2020_SAT",
                                                         "NM2020_VAT"))

Idents(total_Myeloid) <- "cluster_name_slim"
pdf(file = "Dimplot_Myeloid_mergePublic.pdf",width = 16,height = 5)
DimPlot(total_Myeloid,reduction = "umap",label = F,pt.size = 0.15,split.by = "cohort")+theme(legend.position = "bottom")
dev.off()

#cell proportion

cell.prop <- table(total_Myeloid$cluster_name_slim,total_Myeloid$cohort) %>% as.data.frame.matrix()
#cell.prop <- t(t(cell.prop)/colSums(cell.prop)) %>% as.data.frame()
cell.prop$cell <- row.names(cell.prop)
cell.prop <- pivot_longer(cell.prop,cols = setdiff(colnames(cell.prop),"cell"),names_to = "cohort",values_to = "num")
cell.prop = ddply(cell.prop, 'cohort', transform, percent=num/sum(num)*100)

cell.prop$cell <- factor(cell.prop$cell,levels = rev(levels(total_Myeloid$cluster_name_slim)))
cell.prop$cohort <- factor(cell.prop$cohort,levels = levels(total_Myeloid$cohort))

p <- ggplot(cell.prop,aes(x=cohort,y=percent,fill=cell))+
  geom_bar(stat = 'identity',width = 0.8,colour='black')+
  scale_fill_manual(values = rev(hue_pal()(6)))+
  scale_y_continuous(expand = c(0,0))+
  theme_classic()+
  theme(axis.text = element_text(size = 14,colour = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 16,colour = "black"))+
  ylab("Cell ratio (%)")+xlab("");p
ggsave("Cell_proportion_bar_Myeloid_mergePublic.pdf",p,width = 6,height = 3.5)

qsave(total_Myeloid,file = "total_Myeloid_mergePublic.qs")

####################Integrate public and in-house data (B cell)###########################

total_B <- readRDS("total_B.rds")

Idents(total_B) <- "cluster_name_slim"

DimPlot(total_B)

DefaultAssay(total_B) <- "RNA"
total_B[["integrated"]] <- NULL

total_B$cohort <- "InHouse"

data.list1 <- SplitObject(total_B,split.by = "sample")

#public data

#ATVB_2024

data_AVTB2024 <- qread("./PublicSCData/ATVB_2024/total_B.qs")
data_AVTB2024$cohort <- "AVTB2024"
data.list2 <- SplitObject(data_AVTB2024,split.by = "sample")

#Nature_2020

data_Nature2020 <- qread("./PublicSCData/Nature_2020/total_B.qs")
data_Nature2020$cohort <- paste0("Nature2020_",data_Nature2020$group)
data.list3 <- SplitObject(data_Nature2020,split.by = "sample")

#NM_2020

data_NM2020 <- qread("./PublicSCData/NatureMeta_2020/total_B.qs")
data_NM2020$SampleName <- droplevels(data_NM2020$SampleName)
data_NM2020$sample <- data_NM2020$SampleName
data_NM2020$cohort <- paste0("NM2020_",data_NM2020$Tissue)
data.list4 <- SplitObject(data_NM2020,split.by = "sample")

data.list <- c(data.list1,data.list2,data.list3,data.list4)

keep_elements <- sapply(data.list, function(x) {
  # 检查元素是否有列数属性（如矩阵、数据框）
  cols <- tryCatch(ncol(x), error = function(e) NULL)
  # 如果列数存在且大于等于40，则保留（返回TRUE），否则不保留（返回FALSE）
  ifelse(!is.null(cols) && cols >= 40, TRUE, FALSE)
})

data.list <- data.list[keep_elements]

rm(list = setdiff(ls(),"data.list"))

data.list <- lapply(data.list,FUN = function(x){DefaultAssay(x) <- "RNA";x <- SCTransform(x)})

options(future.globals.maxSize= 100000000000000)

data.features <- SelectIntegrationFeatures(object.list = data.list)
data.list <- PrepSCTIntegration(object.list = data.list,assay = "SCT", anchor.features = data.features,verbose = FALSE)
data.anchors <- FindIntegrationAnchors(object.list = data.list, normalization.method = "SCT", anchor.features = data.features,k.filter = 40, verbose = FALSE)
total_B <- IntegrateData(anchorset = data.anchors, normalization.method = "SCT",k.weight = 40, verbose = FALSE)

total_B <- Reduct(total_B)

Idents(total_B) <- "integrated_snn_res.0.5"
DimPlot(total_B,label = T)
DimPlot(total_B,label = T,group.by = "cluster_name_slim")

qsave(total_B,file = "total_B_mergePublic.qs")

#FastMNN
library(SeuratWrappers)

count <- total_B@assays$RNA@counts
meta <- total_B@meta.data
total_B <- CreateSeuratObject(counts = count,meta.data = meta)

data.list <- SplitObject(total_B,split.by = "sample")

for(i in 1:length(data.list)){
  DefaultAssay(data.list[[i]]) <- "RNA"
  data.list[[i]] <- NormalizeData(data.list[[i]])
  data.list[[i]] <- FindVariableFeatures(data.list[[i]])
}

total_B <- RunFastMNN(object.list = data.list)
total_B <- RunUMAP(total_B, reduction = "mnn", dims = 1:30)
total_B <- FindNeighbors(total_B, reduction = "mnn", dims = 1:30) 
total_B <- FindClusters(total_B,resolution = 0.5)

DimPlot(total_B, reduction = "umap", label=T,group.by = "RNA_snn_res.0.5")
DimPlot(total_B,label = T,group.by = "cluster_name_slim")

total_B$cluster_name_slim <- "B_naive"
total_B$cluster_name_slim[total_B$RNA_snn_res.0.5 %in% c(3,4,5,8)] <- "B_memory(unswitched)"
total_B$cluster_name_slim[total_B$RNA_snn_res.0.5 %in% c(1,2,7)] <- "B_memory(switched)"
total_B$cluster_name_slim[total_B$RNA_snn_res.0.5 %in% c(9,10,11,12)] <- "Plasma cell"

total_B$cluster_name_slim <- factor(total_B$cluster_name_slim,
                                    levels = c("B_naive","B_memory(unswitched)","B_memory(switched)","Plasma cell"))
Idents(total_B) <- "cluster_name_slim"

total_B$cohort <- factor(total_B$cohort,levels = c("InHouse","Nature2020_BAT","AVTB2024",
                                                               "NM2020_SAT"))

Idents(total_B) <- "cluster_name_slim"
pdf(file = "Dimplot_B_mergePublic.pdf",width = 15,height = 5)
DimPlot(total_B,reduction = "umap",label = F,pt.size = 0.15,split.by = "cohort",cols = pal[1:10])+theme(legend.position = "bottom")
dev.off()

#cell proportion

cell.prop <- table(total_B$cluster_name_slim,total_B$cohort) %>% as.data.frame.matrix()
#cell.prop <- t(t(cell.prop)/colSums(cell.prop)) %>% as.data.frame()
cell.prop$cell <- row.names(cell.prop)
cell.prop <- pivot_longer(cell.prop,cols = setdiff(colnames(cell.prop),"cell"),names_to = "cohort",values_to = "num")
cell.prop = ddply(cell.prop, 'cohort', transform, percent=num/sum(num)*100)

cell.prop$cell <- factor(cell.prop$cell,levels = rev(levels(total_B$cluster_name_slim)))
cell.prop$cohort <- factor(cell.prop$cohort,levels = levels(total_B$cohort))

p <- ggplot(cell.prop,aes(x=cohort,y=percent,fill=cell))+
  geom_bar(stat = 'identity',width = 0.8,colour='black')+
  scale_fill_manual(values = rev(pal[1:4]))+
  scale_y_continuous(expand = c(0,0))+
  theme_classic()+
  theme(axis.text = element_text(size = 14,colour = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 16,colour = "black"))+
  ylab("Cell ratio (%)")+xlab("");p
ggsave("Cell_proportion_bar_B_mergePublic.pdf",p,width = 6,height = 3.5)

qsave(total_B,file = "total_B_mergePublic.qs")

####################Integrate public and in-house data (NK cell)###########################

total_NK <- readRDS("total_NK.rds")

Idents(total_NK) <- "cluster_name_slim"

DimPlot(total_NK)

DefaultAssay(total_NK) <- "RNA"
total_NK[["integrated"]] <- NULL

total_NK$cohort <- "InHouse"

data.list1 <- SplitObject(total_NK,split.by = "sample")

data_public <- qread("total_NK_Public.qs")

data.list2 <- SplitObject(data_public,split.by = "sample")

data.list <- c(data.list1,data.list2)

keep_elements <- sapply(data.list, function(x) {
  # 检查元素是否有列数属性（如矩阵、数据框）
  cols <- tryCatch(ncol(x), error = function(e) NULL)
  # 如果列数存在且大于等于40，则保留（返回TRUE），否则不保留（返回FALSE）
  ifelse(!is.null(cols) && cols >= 40, TRUE, FALSE)
})

data.list <- data.list[keep_elements]

for(i in 1:length(data.list)){
  DefaultAssay(data.list[[i]]) <- "RNA"
  data.list[[i]] <- NormalizeData(data.list[[i]])
  data.list[[i]] <- FindVariableFeatures(data.list[[i]])
}

total_NK <- RunFastMNN(object.list = data.list)
total_NK <- RunUMAP(total_NK, reduction = "mnn", dims = 1:30)
total_NK <- FindNeighbors(total_NK, reduction = "mnn", dims = 1:30) 
total_NK <- FindClusters(total_NK,resolution = 0.5)

DimPlot(total_NK, reduction = "umap", label=T,group.by = "RNA_snn_res.0.5")
DimPlot(total_NK,label = T,group.by = "cluster_name_slim")

total_NK$cluster_name_slim <- "NK_CD56bright"
total_NK$cluster_name_slim[total_NK$RNA_snn_res.0.5 %in% c(0)] <- "NK_CD56dim"
total_NK$cluster_name_slim[total_NK$RNA_snn_res.0.5 %in% c(4)] <- "NK_tr"

total_NK$cluster_name_slim <- factor(total_NK$cluster_name_slim,
                                     levels = c("NK_CD56bright","NK_CD56dim","NK_tr"))

total_NK$cohort <- factor(total_NK$cohort,levels = c("InHouse","AVTB2024"))

Idents(total_NK) <- "cluster_name_slim"
pdf(file = "Dimplot_NK_mergePublic.pdf",width = 10,height = 5)
DimPlot(total_NK,reduction = "umap",label = F,pt.size = 0.15,split.by = "cohort",cols = pal[1:3])+theme(legend.position = "bottom")
dev.off()

#cell proportion

cell.prop <- table(total_NK$cluster_name_slim,total_NK$cohort) %>% as.data.frame.matrix()
#cell.prop <- t(t(cell.prop)/colSums(cell.prop)) %>% as.data.frame()
cell.prop$cell <- row.names(cell.prop)
cell.prop <- pivot_longer(cell.prop,cols = setdiff(colnames(cell.prop),"cell"),names_to = "cohort",values_to = "num")
cell.prop = ddply(cell.prop, 'cohort', transform, percent=num/sum(num)*100)

cell.prop$cell <- factor(cell.prop$cell,levels = rev(levels(total_NK$cluster_name_slim)))
cell.prop$cohort <- factor(cell.prop$cohort,levels = levels(total_NK$cohort))

p <- ggplot(cell.prop,aes(x=cohort,y=percent,fill=cell))+
  geom_bar(stat = 'identity',width = 0.8,colour='black')+
  scale_fill_manual(values = rev(pal[1:3]))+
  scale_y_continuous(expand = c(0,0))+
  theme_classic()+
  theme(axis.text = element_text(size = 14,colour = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 16,colour = "black"))+
  ylab("Cell ratio (%)")+xlab("");p
ggsave("Cell_proportion_bar_NK_mergePublic.pdf",p,width = 4,height = 3.5)

qsave(total_NK,file = "total_NK_mergePublic.qs")

