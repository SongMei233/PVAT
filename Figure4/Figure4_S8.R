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


######################T cell##############################

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

#Check if there are any cluster that expressed feature of other cell lineages

markers <- c("PTPRC","MS4A1","CD19","BANK1","JCHAIN","SDC1","IGKC","CD3D","CD4","CCR7","SELL","CD8A","GZMK","CCL5","KLRD1","GNLY","NKG7","MKI67","STMN1","TUBA1B","IL1B","S100A8","S100A9","CD68","C1QA","LYVE1","CD1E","LAMP3","FSCN1","PTGDS","IL3RA","LILRA4","F3","CXCL14","CD55")
Idents(total_T) <- "integrated_snn_res.0.5"
DefaultAssay(total_T) <- "SCT"
DotPlot(total_T, features = markers, dot.scale = 6)  + theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) + RotatedAxis()

#extract pure T cell

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

pdf(file = "Figure4A_dimplot_T.pdf",width = 8,height = 6)
DimPlot_scCustom(total_T,colors_use = hue_pal()(6),figure_plot = T,label = T,label.size = 5,pt.size = 0.5)
dev.off()

markers <- c("CD3D",'CD4',"CD8A","IL7R","TCF7","PDCD1","HSPA1A","FOXP3","GZMA")
p <- VlnPlot(total_T,features = markers,pt.size = 0,ncol = 3)
ggsave("Figure4B_VlnPlot_T.pdf",p,width = 8,height = 8)

total_T$disease <- factor(total_T$disease,levels = c("Ctl","CAS"))
Idents(total_T) <- "disease"
pdf(file = "Figure4C_dimplot_disease_T.pdf",width = 7,height = 6)
DimPlot_scCustom(total_T,reduction = "umap",label = F,pt.size = 0.1,colors_use = c("#2b90d9","#f26d5b"),figure_plot = T)
dev.off()

markers <- FindAllMarkers(total_T,logfc.threshold = 0.25,test.use = "wilcox",min.pct = 0.5,min.diff.pct = 0.2,only.pos = T,recorrect_umi = FALSE)
write.csv(markers,file = "total_markers_T.csv")

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
ggsave("Figure4D_cell_proportion_bar_Tcell.pdf",p,width = 5,height = 3.5)


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

ggsave("FigureS8A_Cell proportion between Ctl and CAS_T.pdf",p,width = 7,height = 5)

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
                   filename = "Figure4F_correlation_Tcell.pdf")

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
ggsave("Figure4G_DEG_Reactome_network_up_10.pdf",p,width = 7,height = 6)

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
ggsave("FigureS8B_DEG_Tcell_network_down_10.pdf",p,width = 7,height = 6)


#Cytokines
cytokine <- c("IL1A","IL1B","IL1RN","IL18","IL2","IL4","IL7","IL9","IL13","IL15","IL3",
              "IL5","CSF2","IL6","IL11","CSF3","IL12A","IL12B","LIF","OSM","IL10","IL20","TXLNA",
              "IL16","IL17A","IL17B","IFNA1","IFNB1","IFNG","CD40LG","LTB","TNF","TNFSF9","TNFSF13",
              "CD70","TNFSF8","FASLG","TNFSF18","TNFSF14","TNFSF4","TNFSF13B","TNFSF10",
              "TNFSF12","TNFSF11","TGFB1","TGFB2","TGFB3","EPO","TPO","FLT3LG","KITLG",
              "CSF1","MST1")

DefaultAssay(total_T) <- "SCT"
exp <- FetchData(total_T,vars = cytokine)
exp <- exp[,colnames(exp)[!str_detect(colnames(exp),pattern = "rna")]]
total_T[["cytokine"]] <- CreateAssayObject(t(exp))

total_T$cell_sample <- paste(total_T$cluster_name_slim,total_T$sample,sep = "-")
matrix <- AverageExpression(total_T,assays = "cytokine",group.by = "cell_sample") %>% as.data.frame()

meta <- total_T@meta.data[,c("sample","disease","condition")] %>% distinct()
meta <- remove_rownames(meta) %>% column_to_rownames(.,var = "sample")

anno <- data.frame(item = colnames(matrix) %>% gsub("cytokine.","",.),row.names = colnames(matrix))
anno$sample <- strsplit(anno$item,split = "\\.") %>% lapply(.,function(x){x[2]}) %>% unlist()
anno$cluster <- strsplit(anno$item,split = "\\.") %>% lapply(.,function(x){x[1]}) %>% unlist()
anno$disease <- meta[anno$sample,"disease"]
anno$condition <- meta[anno$sample,"condition"]

#disease

matrix_merge <- AverageExpression(total_T,assays = "cytokine",group.by = "cell_disease") %>% as.data.frame()

test <- matrix(nrow = nrow(matrix_merge),ncol = ncol(matrix_merge)) %>%
  as.data.frame(row.names = row.names(matrix_merge))
colnames(test) <- colnames(matrix_merge)

for(x in unique(total_T$cluster_name_slim)){
  print(x)
  
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
  
  x <- gsub("_",".",x)
  x <- gsub(" ",".",x)
  x <- gsub("\\+",".",x)
  x <- gsub("\\(",".",x)
  x <- gsub("\\)",".",x)
  test[cytokine,paste("cytokine",x,"CAS",sep = ".")] <- 1
}
test[is.na(test)] <- 0

test <- test[rowSums(test)>0,]
matrix_merge <- matrix_merge[row.names(test),]
test[test==1] <- "*"
test[test==0] <- ""

anno_col <- data.frame(Group = rep(levels(total_T$disease),6),
                       Cell = rep(levels(total_T$cluster_name_slim),each=2),
                       row.names = colnames(matrix_merge))
anno_color <- list(Cell = c(CD4Tmemory = hue_pal()(6)[1],
                            "CD4Texhausted" = hue_pal()(6)[2],
                            "CD4Tstress" = hue_pal()(6)[3],
                            "CD4Treg" = hue_pal()(6)[4],
                            "CD8Tmemory" = hue_pal()(6)[5],
                            "CD8Teffector" = hue_pal()(6)[6]),
                   
                   Group = c(Ctl = "#2b90d9", CAS = "#f26d5b"))

pheatmap::pheatmap(matrix_merge,scale = "row",cluster_cols = F,
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
                   filename = "Figure4H_cytokine_heatmap_T.pdf")

anno$cell_disease <- paste(anno$cluster,anno$disease,sep = "-")
anno$cell_disease <- factor(anno$cell_disease,levels = levels(total_T$cell_disease))
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
                   gaps_col = seq(10,50,10),
                   border_color = "white",
                   width = 12,
                   cellwidth = 7,cellheight = 12,
                   fontsize_row = 12,
                   treeheight_row = 20,
                   filename = "Figure4H_cytokine_heatmap_T_perSample.pdf")

#condition
total_T_sub <- subset(total_T,subset = condition != "Ctl")
total_T_sub$condition <- droplevels(total_T_sub$condition)
matrix_merge <- AverageExpression(total_T_sub,assays = "cytokine",group.by = "cell_condition") %>% as.data.frame()

test <- matrix(nrow = nrow(matrix_merge),ncol = ncol(matrix_merge)) %>%
  as.data.frame(row.names = row.names(matrix_merge))
colnames(test) <- colnames(matrix_merge)

for(x in unique(total_T$cluster_name_slim)){
  print(x)
  
  sample1 <- anno[anno$condition=="CAS_NS","sample"] %>% unique()
  index1 <- paste("cytokine",x,sample1,sep = ".")
  
  sample2 <- anno[anno$condition=="CAS_Smo","sample"] %>% unique()
  index2 <- paste("cytokine",x,sample2,sep = ".")
  
  cytokine <- c()
  
  for(y in row.names(matrix)){
    data1 <- matrix[y,index1] %>% as.numeric()
    data2 <- matrix[y,index2] %>% as.numeric()
    t <- t.test(data1,data2)
    if(!is.na(t$p.value)&t$p.value<0.05){cytokine <- c(cytokine,y)}
  }
  
  x <- gsub("_",".",x)
  x <- gsub(" ",".",x)
  x <- gsub("\\+",".",x)
  x <- gsub("\\(",".",x)
  x <- gsub("\\)",".",x)
  test[cytokine,paste("cytokine",x,"CAS.Smo",sep = ".")] <- 1
}
test[is.na(test)] <- 0

test <- test[rowSums(test)>0,]
matrix_merge <- matrix_merge[row.names(test),]
test[test==1] <- "*"
test[test==0] <- ""
anno_col <- data.frame(Group = rep(levels(total_T_sub$condition),6),
                       Cell = rep(levels(total_T_sub$cluster_name_slim),each=2),
                       row.names = colnames(matrix_merge))
anno_color <- list(Cell = c(CD4Tmemory = hue_pal()(6)[1],
                            "CD4Texhausted" = hue_pal()(6)[2],
                            "CD4Tstress" = hue_pal()(6)[3],
                            "CD4Treg" = hue_pal()(6)[4],
                            "CD8Tmemory" = hue_pal()(6)[5],
                            "CD8Teffector" = hue_pal()(6)[6]),
                   
                   Group = c(CAS_NS = "#2b90d9", CAS_Smo = "#f26d5b"))

pheatmap::pheatmap(matrix_merge,scale = "row",cluster_cols = F,
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
                   filename = "FigureS8C_cytokine_heatmap_T.pdf")

anno$cell_condition <- paste(anno$cluster,anno$condition,sep = "-")
anno$cell_condition <- factor(anno$cell_condition,levels = levels(total_T$cell_condition))
anno <- anno[order(anno$cell_condition,decreasing = F),]
anno <- filter(anno,condition != "Ctl")

matrix_per <- matrix[row.names(matrix_merge),row.names(anno)]

anno_col <- data.frame(Group = anno$condition,
                       Cell = anno$cluster,
                       row.names = colnames(matrix_per))

bk <- seq(-2,2,by=0.01)
pheatmap::pheatmap(matrix_per,scale = "row",cluster_cols = F,
                   color = colorRampPalette(c("navy","white","firebrick"))(length(bk)),
                   legend_breaks=seq(-2,2,1),breaks=bk,
                   annotation_col = anno_col,
                   annotation_colors = anno_color,
                   show_colnames = F,
                   gaps_col = seq(7,35,7),
                   border_color = "white",
                   width = 12,height = 6,
                   cellwidth = 7,cellheight = 12,
                   fontsize_row = 12,
                   treeheight_row = 20,
                   filename = "FigureS8C_cytokine_heatmap_T_perSample.pdf")


######################B cell##############################

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

#Check if there are any cluster that expressed feature of other cell lineages

markers <- c("PTPRC","MS4A1","CD19","BANK1","JCHAIN","SDC1","IGKC","CD3D","CD4","CCR7","SELL","CD8A","GZMK","CCL5","KLRD1","GNLY","NKG7","MKI67","STMN1","TUBA1B","IL1B","S100A8","S100A9","CD68","C1QA","LYVE1","CD1E","LAMP3","FSCN1","PTGDS","IL3RA","LILRA4","F3","CXCL14","CD55")
Idents(total_B) <- "integrated_snn_res.0.5"
DefaultAssay(total_B) <- "SCT"
DotPlot(total_B, features = markers, dot.scale = 6)  + theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) + RotatedAxis()

#extract pure B cells

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

pdf(file = "Figure4I_dimplot_B.pdf",width = 8,height = 6)
DimPlot_scCustom(total_B,colors_use = pal[1:10],figure_plot = T,label = T,label.size = 5,pt.size = 0.5)
dev.off()

markers <- c("IL4R","IGHM","IGHD","CD27","CD38","MZB1")

p <- VlnPlot(total_B,features = markers,stack = T,group.by = "cluster_name_slim",cols = pal[1:4],fill.by = "ident")+
  theme(strip.text.x = element_text(angle = 0,face = "plain",size = 16),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 16),
        legend.position = "none");p
ggsave("Figure4J_ViolinPlot.pdf",p,width = 9,height = 4)


total_B$disease <- factor(total_B$disease,levels = c("Ctl","CAS"))
Idents(total_B) <- "disease"
pdf(file = "Figure4K_dimplot_disease_B.pdf",width = 7,height = 6)
DimPlot_scCustom(total_B,reduction = "umap",label = F,pt.size = 0.1,colors_use = c("#2b90d9","#f26d5b"),figure_plot = T)
dev.off()

markers <- FindAllMarkers(total_B,logfc.threshold = 0.25,test.use = "wilcox",min.pct = 0.5,min.diff.pct = 0.2,only.pos = T,recorrect_umi = FALSE)
write.csv(markers,file = "total_markers_B.csv")

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
ggsave("Figure4L_cell_proportion_bar_Bcell.pdf",p,width = 6,height = 3)


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

ggsave("FigureS8E_Cell proportion between Ctl and CAS_innate.pdf",p,width = 5,height = 5)

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
                   filename = "FigureS8F_correlation_Bcell.pdf")


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
ggsave("Figure4M_DEG_Reactome_network_up_20.pdf",p,width = 7,height = 6)

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
ggsave("FigureS8G_DEG_Bcell_network_up_20.pdf",p,width = 7,height = 6)



#Cytokines
cytokine <- c("IL1A","IL1B","IL1RN","IL18","IL2","IL4","IL7","IL9","IL13","IL15","IL3",
              "IL5","CSF2","IL6","IL11","CSF3","IL12A","IL12B","LIF","OSM","IL10","IL20","TXLNA",
              "IL16","IL17A","IL17B","IFNA1","IFNB1","IFNG","CD40LG","LTB","TNF","TNFSF9","TNFSF13",
              "CD70","TNFSF8","FASLG","TNFSF18","TNFSF14","TNFSF4","TNFSF13B","TNFSF10",
              "TNFSF12","TNFSF11","TGFB1","TGFB2","TGFB3","EPO","TPO","FLT3LG","KITLG",
              "CSF1","MST1")

DefaultAssay(total_B) <- "SCT"
exp <- FetchData(total_B,vars = cytokine)
exp <- exp[,colnames(exp)[!str_detect(colnames(exp),pattern = "rna")]]
total_B[["cytokine"]] <- CreateAssayObject(t(exp))

total_B$cell_sample <- paste(total_B$cluster_name_slim,total_B$sample,sep = "-")
matrix <- AverageExpression(total_B,assays = "cytokine",group.by = "cell_sample") %>% as.data.frame()

meta <- total_B@meta.data[,c("sample","disease","condition")] %>% distinct()
meta <- remove_rownames(meta) %>% column_to_rownames(.,var = "sample")

anno <- data.frame(item = colnames(matrix) %>% gsub("cytokine.","",.),row.names = colnames(matrix))
anno$item <- gsub("..CAS",".CAS",anno$item) %>% gsub("..Ctl",".Ctl",.)
anno$item <- gsub(".CAS","_CAS",anno$item) %>% gsub(".Ctl","_Ctl",.)
anno$sample <- strsplit(anno$item,split = "_") %>% lapply(.,function(x){x[2]}) %>% unlist()
anno$cluster <- strsplit(anno$item,split = "_") %>% lapply(.,function(x){x[1]}) %>% unlist()
anno$disease <- meta[anno$sample,"disease"]
anno$condition <- meta[anno$sample,"condition"]

#disease

matrix_merge <- AverageExpression(total_B,assays = "cytokine",group.by = "cell_disease") %>% as.data.frame()

test <- matrix(nrow = nrow(matrix_merge),ncol = ncol(matrix_merge)) %>%
  as.data.frame(row.names = row.names(matrix_merge))
colnames(test) <- colnames(matrix_merge)

for(x in unique(total_B$cluster_name_slim)){
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

anno_col <- data.frame(Group = rep(levels(total_B$disease),4),
                       Cell = rep(levels(total_B$cluster_name_slim),each=2),
                       row.names = colnames(matrix_merge))
anno_color <- list(Cell = c(B_naive = pal[1],
                            "B_memory(unswitched)" = pal[2],
                            "B_memory(switched)" = pal[3],
                            "Plasma cell" = pal[4]),
                   Group = c(Ctl = "#2b90d9", CAS = "#f26d5b"))

pheatmap::pheatmap(matrix_merge,scale = "row",cluster_cols = F,
                   annotation_col = anno_col,
                   annotation_colors = anno_color,
                   display_numbers = test,
                   fontsize_number = 12,
                   show_colnames = F,
                   gaps_col = seq(2,6,2),
                   width = 12,height = 6,
                   cellwidth = 16,cellheight = 12,
                   fontsize_row = 14,
                   treeheight_row = 20,
                   filename = "Figure4N_cytokine_heatmap_B.pdf")

anno$cluster[anno$cluster=="B.memory.switched"] <- "B_memory(switched)"
anno$cluster[anno$cluster=="B.memory.unswitched"] <- "B_memory(unswitched)"
anno$cluster[anno$cluster=="B.naiv"] <- "B_naive"
anno$cluster[anno$cluster=="Plasma.cel"] <- "Plasma cell"

anno$cell_disease <- paste(anno$cluster,anno$disease,sep = "-")
anno$cell_disease <- factor(anno$cell_disease,levels = levels(total_B$cell_disease))
anno <- anno[order(anno$cell_disease,decreasing = F),]

matrix_per <- matrix[row.names(matrix_merge),row.names(anno)]

anno_col <- data.frame(Group = anno$disease,
                       Cell = anno$cluster,
                       row.names = colnames(matrix_per))
anno_color <- list(Cell = c(B_naive = pal[1],
                            "B_memory(unswitched)" = pal[2],
                            "B_memory(switched)" = pal[3],
                            "Plasma cell" = pal[4]),
                   Group = c(Ctl = "#2b90d9", CAS = "#f26d5b"))

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
                   filename = "Figure4N_cytokine_heatmap_B_perSample.pdf")

#condition
total_B_sub <- subset(total_B,subset = condition != "Ctl")
total_B_sub$condition <- droplevels(total_B_sub$condition)
matrix_merge <- AverageExpression(total_B_sub,assays = "cytokine",group.by = "cell_condition") %>% as.data.frame()

test <- matrix(nrow = nrow(matrix_merge),ncol = ncol(matrix_merge)) %>%
  as.data.frame(row.names = row.names(matrix_merge))
colnames(test) <- colnames(matrix_merge)

for(x in unique(total_B$cluster_name_slim)){
  print(x)
  
  x <- gsub("_",".",x)
  x <- gsub(" ",".",x)
  x <- gsub("\\+",".",x)
  x <- gsub("\\(",".",x)
  x <- gsub("\\)",".",x)
  
  sample1 <- anno[anno$condition=="CAS_NS","sample"] %>% unique()
  index1 <- paste("cytokine",x,sample1,sep = ".")
  
  sample2 <- anno[anno$condition=="CAS_Smo","sample"] %>% unique()
  index2 <- paste("cytokine",x,sample2,sep = ".")
  
  cytokine <- c()
  
  for(y in row.names(matrix)){
    data1 <- matrix[y,index1] %>% as.numeric()
    data2 <- matrix[y,index2] %>% as.numeric()
    t <- t.test(data1,data2)
    if(!is.na(t$p.value)&t$p.value<0.05){cytokine <- c(cytokine,y)}
  }
  

  test[cytokine,paste("cytokine",x,"CAS.Smo",sep = ".")] <- 1
}
test[is.na(test)] <- 0

test <- test[rowSums(test)>0,]
matrix_merge <- matrix_merge[row.names(test),]
test[test==1] <- "*"
test[test==0] <- ""

anno_col <- data.frame(Group = rep(levels(total_B_sub$condition),4),
                       Cell = rep(levels(total_B_sub$cluster_name_slim),each=2),
                       row.names = colnames(matrix_merge))
anno_color <- list(Cell = c(B_naive = pal[1],
                            "B_memory(unswitched)" = pal[2],
                            "B_memory(switched)" = pal[3],
                            "Plasma cell" = pal[4]),
                   Group = c(CAS_NS = "#2b90d9", CAS_Smo = "#f26d5b"))

pheatmap::pheatmap(matrix_merge,scale = "row",cluster_cols = F,
                   annotation_col = anno_col,
                   annotation_colors = anno_color,
                   display_numbers = test,
                   fontsize_number = 12,
                   show_colnames = F,
                   gaps_col = seq(2,6,2),
                   width = 12,height=6,
                   cellwidth = 16,cellheight = 12,
                   fontsize_row = 14,
                   treeheight_row = 20,
                   filename = "FigureS8H_cytokine_heatmap_B.pdf")

anno$cell_condition <- paste(anno$cluster,anno$condition,sep = "-")
anno$cell_condition <- factor(anno$cell_condition,levels = levels(total_B$cell_condition))
anno <- anno[order(anno$cell_condition,decreasing = F),]
anno <- filter(anno,condition != "Ctl")

matrix_per <- matrix[row.names(matrix_merge),row.names(anno)]

anno_col <- data.frame(Group = anno$condition,
                       Cell = anno$cluster,
                       row.names = colnames(matrix_per))

bk <- seq(-2,2,by=0.01)
pheatmap::pheatmap(matrix_per,scale = "row",cluster_cols = F,
                   color = colorRampPalette(c("navy","white","firebrick"))(length(bk)),
                   legend_breaks=seq(-2,2,1),breaks=bk,
                   annotation_col = anno_col,
                   annotation_colors = anno_color,
                   show_colnames = F,
                   gaps_col = seq(7,21,7),
                   border_color = "white",
                   width = 12,height = 6,
                   cellwidth = 7,cellheight = 12,
                   fontsize_row = 12,
                   treeheight_row = 20,
                   filename = "FigureS8H_cytokine_heatmap_B_perSample.pdf")
