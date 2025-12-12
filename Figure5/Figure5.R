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

######################Annotated all innate immune cell##############################

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

#Check if there are any cluster that expressed feature of other cell lineages

markers <- c("PTPRC","MS4A1","CD19","BANK1","JCHAIN","SDC1","IGKC","CD3D","CD4","CCR7","SELL","CD8A","GZMK","CCL5","KLRD1","GNLY","NKG7","MKI67","STMN1","TUBA1B","IL1B","S100A8","S100A9","CD68","C1QA","LYVE1","CD1E","LAMP3","FSCN1","PTGDS","IL3RA","LILRA4","F3","CXCL14","CD55")
Idents(total_immune) <- "integrated_snn_res.1.5"
DefaultAssay(total_immune) <- "SCT"
DotPlot(total_immune, features = markers, dot.scale = 6)  + theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) + RotatedAxis()

markers <- FindMarkers(total_immune,ident.1 = c(4,10),logfc.threshold = 0.25,test.use = "wilcox",min.pct = 0.5,min.diff.pct = 0.2,only.pos = T,recorrect_umi = FALSE)

#extract pure innate immune cell

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

######################## NK cell##############################

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

#Check if there are any cluster that expressed feature of other cell lineages

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

pdf(file = "Figure5A_dimplot_NK.pdf",width = 6,height = 4)
DimPlot_scCustom(total_NK,colors_use = pal[1:3],figure_plot = T,label = T,label.size = 5,pt.size = 0.5)
dev.off()

total_NK$disease <- factor(total_NK$disease,levels = c("Ctl","CAS"))
Idents(total_NK) <- "disease"
pdf(file = "Figure5B_dimplot_NK_disease.pdf",width = 5.5,height = 4)
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
ggsave("Figure5C_cell_proportion_bar_NK.pdf",p,width = 5,height = 3.5)


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
                   filename = "Figure5D_NK_correlation.pdf")


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
ggsave("Figure5E_NK_DEG_Reactome_network_down_10.pdf",p,width = 7,height = 6)


#Cytokines
cytokine <- c("IL1A","IL1B","IL1RN","IL18","IL2","IL4","IL7","IL9","IL13","IL15","IL3",
              "IL5","CSF2","IL6","IL11","CSF3","IL12A","IL12B","LIF","OSM","IL10","IL20","TXLNA",
              "IL16","IL17A","IL17B","IFNA1","IFNB1","IFNG","CD40LG","LTB","TNF","TNFSF9","TNFSF13",
              "CD70","TNFSF8","FASLG","TNFSF18","TNFSF14","TNFSF4","TNFSF13B","TNFSF10",
              "TNFSF12","TNFSF11","TGFB1","TGFB2","TGFB3","EPO","TPO","FLT3LG","KITLG",
              "CSF1","MST1")

DefaultAssay(total_NK) <- "SCT"
exp <- FetchData(total_NK,vars = cytokine)
exp <- exp[,colnames(exp)[!str_detect(colnames(exp),pattern = "rna")]]
total_NK[["cytokine"]] <- CreateAssayObject(t(exp))

total_NK$cell_sample <- paste(total_NK$cluster_name_slim,total_NK$sample,sep = "-")
matrix <- AverageExpression(total_NK,assays = "cytokine",group.by = "cell_sample") %>% as.data.frame()

meta <- total_NK@meta.data[,c("sample","disease","condition")] %>% distinct()
meta <- remove_rownames(meta) %>% column_to_rownames(.,var = "sample")

anno <- data.frame(item = colnames(matrix) %>% gsub("cytokine.","",.),row.names = colnames(matrix))
#anno$item <- gsub("..CAS",".CAS",anno$item) %>% gsub("..Ctl",".Ctl",.)
anno$item <- gsub(".CAS","_CAS",anno$item) %>% gsub(".Ctl","_Ctl",.)
anno$sample <- strsplit(anno$item,split = "_") %>% lapply(.,function(x){x[2]}) %>% unlist()
anno$cluster <- strsplit(anno$item,split = "_") %>% lapply(.,function(x){x[1]}) %>% unlist()
anno$disease <- meta[anno$sample,"disease"]
anno$condition <- meta[anno$sample,"condition"]

#disease

matrix_merge <- AverageExpression(total_NK,assays = "cytokine",group.by = "cell_disease") %>% as.data.frame()

test <- matrix(nrow = nrow(matrix_merge),ncol = ncol(matrix_merge)) %>%
  as.data.frame(row.names = row.names(matrix_merge))
colnames(test) <- colnames(matrix_merge)

for(x in unique(total_NK$cluster_name_slim)){
  print(x)
  
  x <- gsub("_",".",x)
  x <- gsub(" ",".",x)
  x <- gsub("\\+",".",x)
  x <- gsub("\\(",".",x)
  x <- gsub("\\)",".",x)
  
  sample1 <- anno[anno$disease=="CAS","sample"] %>% unique()
  index1 <- paste("cytokine",x,sample1,sep = ".")
  
  sample2 <- anno[anno$disease=="Ctl","sample"] %>% unique()
  index2 <- paste("cytokine",x,sample2,sep = ".")
  
  cytokine <- c()
  
  for(y in row.names(matrix)){
    data1 <- matrix[y,index1] %>% as.numeric()
    data2 <- matrix[y,index2] %>% as.numeric()
    t <- t.test(data1,data2)
    if(!is.na(t$p.value)&t$p.value<0.05){cytokine <- c(cytokine,y)}
  }
  
  
  test[cytokine,paste("cytokine",x,"CAS",sep = ".")] <- 1
}
test[is.na(test)] <- 0

test <- test[rowSums(test)>0,]
matrix_merge <- matrix_merge[row.names(test),]
test[test==1] <- "*"
test[test==0] <- ""

anno_col <- data.frame(Group = rep(levels(total_NK$disease),3),
                       Cell = rep(levels(total_NK$cluster_name_slim),each=2),
                       row.names = colnames(matrix_merge))
anno_color <- list(Cell = c(NK_CD56bright = pal[1],
                            NK_CD56dim = pal[2],
                            NK_tr = pal[3]),
                   Group = c(Ctl = "#2b90d9", CAS = "#f26d5b"))

pheatmap::pheatmap(matrix_merge,scale = "row",cluster_cols = F,
                   annotation_col = anno_col,
                   annotation_colors = anno_color,
                   display_numbers = test,
                   fontsize_number = 12,
                   show_colnames = F,
                   gaps_col = c(2,4),
                   width = 12,height = 6,
                   cellwidth = 16,cellheight = 12,
                   fontsize_row = 14,
                   treeheight_row = 20,
                   filename = "Figure5F_NK_cytokine_heatmap.pdf")

anno$cluster[anno$cluster=="NK.CD56bright"] <- "NK_CD56bright"
anno$cluster[anno$cluster=="NK.CD56dim"] <- "NK_CD56dim"
anno$cluster[anno$cluster=="NK.tr"] <- "NK_tr"

anno$cell_disease <- paste(anno$cluster,anno$disease,sep = "-")
anno$cell_disease <- factor(anno$cell_disease,levels = levels(total_NK$cell_disease))
anno <- anno[order(anno$cell_disease,decreasing = F),]

matrix_per <- matrix[row.names(matrix_merge),row.names(anno)]

anno_col <- data.frame(Group = anno$disease,
                       Cell = anno$cluster,
                       row.names = colnames(matrix_per))

bk <- seq(-2,2,by=0.01)
pheatmap::pheatmap(matrix_per,scale = "row",cluster_cols = F,
                   color = colorRampPalette(c("navy","white","firebrick"))(length(bk)),
                   legend_breaks=seq(-2,2,1),breaks=bk,
                   annotation_col = anno_col,
                   annotation_colors = anno_color,
                   show_colnames = F,
                   gaps_col = seq(10,30,10),
                   border_color = "white",
                   width = 12,height = 6,
                   cellwidth = 7,cellheight = 12,
                   fontsize_row = 12,
                   treeheight_row = 20,
                   filename = "Figure5F_NK_cytokine_heatmap_perSample.pdf")

