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

total_myeloid <- subset(total_immune,subset = cluster_name_slim %in% c("cDC1","cDC2","DC_LAMP3+","Monocyte","Neutrophil","Macrophage","pDC"))

total_myeloid <- Reduct(total_myeloid)

pdf(file = "Figure6A_dimplot_innate_immune.pdf",width = 8,height = 6)
DimPlot_scCustom(total_myeloid,colors_use = pal[1:10],figure_plot = T,label = T,label.size = 5,pt.size = 0.5)
dev.off()

Idents(total_myeloid) <- "disease"
pdf(file = "Figure6B_dimplot_disease.pdf",width = 7,height = 6)
DimPlot_scCustom(total_myeloid,reduction = "umap",label = F,pt.size = 0.1,colors_use = c("#2b90d9","#f26d5b"),figure_plot = T)
dev.off()

pdf(file = "FigureS9D_VlnPlot_myeloid.pdf",width = 12,height = 6)
VlnPlot(total_myeloid,features = c("HLA-DRA","XCR1","CD1C","LAMP3","LILRB2","CSF3R","CD163","IL3RA"),pt.size = 0,cols = pal[1:10])
dev.off()

#cell proportion in bar plot

cell.prop <- table(total_myeloid$cluster_name_slim,total_myeloid$disease) %>% as.data.frame.matrix()
#cell.prop <- t(t(cell.prop)/colSums(cell.prop)) %>% as.data.frame()
cell.prop$cell <- row.names(cell.prop)
cell.prop <- pivot_longer(cell.prop,cols = c("Ctl","CAS"),names_to = "disease",values_to = "num")
cell.prop = ddply(cell.prop, 'cell', transform, percent=num/sum(num)*100)

cell.prop$cell <- factor(cell.prop$cell,levels = rev(c("cDC1","cDC2","DC_LAMP3+","Monocyte","Neutrophil","Macrophage","pDC")))
cell.prop$disease <- factor(cell.prop$disease,levels = c("CAS","Ctl"))

p <- ggplot(cell.prop,aes(x=percent,y=cell,fill=disease))+
  geom_bar(stat = 'identity',width = 0.8,colour='black')+
  geom_vline(aes(xintercept = 28.84), colour="black", linetype="dashed",linewidth = 1)+
  scale_fill_manual(values = c("#f26d5b","#2b90d9"))+
  theme_classic()+
  theme(axis.text = element_text(size = 12,colour = "black"),
        axis.title.x = element_text(size = 14,colour = "black"))+
  ylab("")+xlab("Cell ratio (%)");p
ggsave("Figure6C_cell_proportion_bar.pdf",p,width = 5,height = 3.5)



#cell proportion box plot
cell.prop <- as.data.frame.matrix((table(total_myeloid$cluster_name_slim, total_myeloid$sample)))
cell.prop <- t(cell.prop)/colSums(cell.prop)
cell.prop <- as.data.frame(cell.prop)
cell.prop$Patient <- row.names(cell.prop)
cell.prop <- pivot_longer(cell.prop,cols = setdiff(colnames(cell.prop),"Patient"),names_to = "Cell",values_to = "Freq")

meta <- total_myeloid@meta.data[,c("sample",'disease')] %>% distinct()
row.names(meta) <- meta$sample
cell.prop$disease <- factor(meta[cell.prop$Patient,"disease"],
                            levels = c("Ctl","CAS"))
cell.prop$Cell <- factor(cell.prop$Cell,
                         levels = levels(total_myeloid$cluster_name_slim))

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

ggsave("FigureS9E_Cell proportion between Ctl and CAS_myeloid.pdf",p,width = 8,height = 6)

#Roe
meta.tb <- total_myeloid@meta.data
meta.tb$disease <- as.character(meta.tb$disease) %>% factor(.,levels=c("Ctl","CAS"))
OR.list <- do.tissueDist(cellInfo.tb=meta.tb,
                         out.prefix="./Roe",
                         meta.cluster=meta.tb$cluster_name_slim,
                         loc=meta.tb$disease,
                         pdf.width=4,pdf.height=6,verbose=1)

matrix <- OR.list[["OR.dist.tb"]]
matrix <- column_to_rownames(matrix,var = "rid")
matrix <- matrix[levels(total_myeloid$cluster_name_slim),]
write.csv(matrix,file = "Figure4_Roe_res_innate.csv")

label <- matrix
label[label>1] <- "+++"
label[label<1&label>0.8] <- "++"
label[label<0.8&label>0.2] <- "+"
label[label<0.2&label>0] <- "+/-"

#Correlation heatmap
library(psych)

total_SCT_mat <- as.matrix(total_myeloid@assays$SCT@data)
total_SCT_mat_ClusterMeanSplit <- apply(total_SCT_mat,1,function(input){aggregate(input,list(total_myeloid@meta.data$cluster_name_slim),mean)[,2]})
rownames(total_SCT_mat_ClusterMeanSplit) <- levels(total_myeloid@meta.data$cluster_name_slim)
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
                   filename = "Figure6D_correlation.pdf")

#DEGs
Idents(total_myeloid) <- "cell_disease"

diffExp <- function(total,ident1,ident2){
  DefaultAssay(total) <- "SCT"
  Cell_SCT_diff <- FindMarkers(total, ident.1 = ident1, ident.2 = ident2, only.pos = F, test.use="wilcox",logfc.threshold = 0.2,min.pct = 0.4,min.diff.pct = 0.1,recorrect_umi = FALSE)
  return(Cell_SCT_diff)
}

cluster1 <- levels(total_myeloid$cluster_name_slim)
DefaultAssay(total_myeloid) <- "SCT"
data.features <- row.names(total_myeloid)
RB_genes = c(grep("^RPS", data.features, v=T),grep("^RPL", data.features, v=T))
MT_genes = grep("^MT-", data.features, v=T)
remove <- c(RB_genes,MT_genes)
total_remove <- subset(total_myeloid,features = setdiff(data.features,remove))


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
cluster1 <- levels(total_myeloid$cluster_name_slim)

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
ggsave("Figure6E_DEG_Reactome_network_down_10.pdf",p,width = 7,height = 6)

#Cytokines
cytokine <- c("IL1A","IL1B","IL1RN","IL18","IL2","IL4","IL7","IL9","IL13","IL15","IL3",
              "IL5","CSF2","IL6","IL11","CSF3","IL12A","IL12B","LIF","OSM","IL10","IL20","TXLNA",
              "IL16","IL17A","IL17B","IFNA1","IFNB1","IFNG","CD40LG","LTB","TNF","TNFSF9","TNFSF13",
              "CD70","TNFSF8","FASLG","TNFSF18","TNFSF14","TNFSF4","TNFSF13B","TNFSF10",
              "TNFSF12","TNFSF11","TGFB1","TGFB2","TGFB3","EPO","TPO","FLT3LG","KITLG",
              "CSF1","MST1")

DefaultAssay(total_myeloid) <- "SCT"
exp <- FetchData(total_myeloid,vars = cytokine)
exp <- exp[,colnames(exp)[!str_detect(colnames(exp),pattern = "rna")]]
total_myeloid[["cytokine"]] <- CreateAssayObject(t(exp))

total_myeloid <- subset(total_myeloid,subset = cluster_name_slim != "Monocyte")

total_myeloid$cell_sample <- paste(total_myeloid$cluster_name_slim,total_myeloid$sample,sep = "-")
matrix <- AverageExpression(total_myeloid,assays = "cytokine",group.by = "cell_sample") %>% as.data.frame()

meta <- total_myeloid@meta.data[,c("sample","disease","condition")] %>% distinct()
meta <- remove_rownames(meta) %>% column_to_rownames(.,var = "sample")

anno <- data.frame(item = colnames(matrix) %>% gsub("cytokine.","",.),row.names = colnames(matrix))
anno$item <- gsub("\\..CAS",".CAS",anno$item) %>% gsub("\\..Ctl",".Ctl",.)
anno$item <- gsub("\\.CAS","_CAS",anno$item) %>% gsub("\\.Ctl","_Ctl",.)
anno$sample <- strsplit(anno$item,split = "_") %>% lapply(.,function(x){x[2]}) %>% unlist()
anno$cluster <- strsplit(anno$item,split = "_") %>% lapply(.,function(x){x[1]}) %>% unlist()
anno$disease <- meta[anno$sample,"disease"]
anno$condition <- meta[anno$sample,"condition"]

#disease

matrix_merge <- AverageExpression(total_myeloid,assays = "cytokine",group.by = "cell_disease") %>% as.data.frame()

test <- matrix(nrow = nrow(matrix_merge),ncol = ncol(matrix_merge)) %>%
  as.data.frame(row.names = row.names(matrix_merge))
colnames(test) <- colnames(matrix_merge)

for(x in unique(total_myeloid$cluster_name_slim)){
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


total_myeloid$cluster_name_slim <- droplevels(total_myeloid$cluster_name_slim)
anno_col <- data.frame(Group = rep(levels(total_myeloid$disease),6),
                       Cell = rep(levels(total_myeloid$cluster_name_slim),each=2),
                       row.names = colnames(matrix_merge))
anno_color <- list(Cell = c(cDC1 = pal[4],
                            cDC2 = pal[5],
                            "DC_LAMP3+" = pal[6],
                            Monocyte = pal[7],
                            Neutrophil = pal[8],
                            Macrophage = pal[9],
                            pDC = pal[10]),
                   Group = c(Ctl = "#2b90d9", CAS = "#f26d5b"))

pheatmap::pheatmap(matrix_merge,scale = "row",cluster_cols = F,
                   annotation_col = anno_col,
                   annotation_colors = anno_color,
                   display_numbers = test,
                   fontsize_number = 12,
                   show_colnames = F,
                   gaps_col = seq(2,12,2),
                   width = 12,height = 6,
                   cellwidth = 16,cellheight = 12,
                   fontsize_row = 14,
                   treeheight_row = 20,
                   filename = "Figure6F_cytokine_heatmap.pdf")

anno$cluster[anno$cluster=="DC.LAMP3"] <- "DC_LAMP3+"

anno$cell_disease <- paste(anno$cluster,anno$disease,sep = "-")
anno$cell_disease <- factor(anno$cell_disease,levels = levels(total_myeloid$cell_disease))
anno <- anno[order(anno$cell_disease,decreasing = F),]

matrix_per <- matrix[row.names(matrix_merge),row.names(anno)]

anno_col <- data.frame(Group = anno$disease,
                       Cell = anno$cluster,
                       row.names = colnames(matrix_per))

anno_color <- list(Cell = c(cDC1 = pal[4],
                            cDC2 = pal[5],
                            "DC_LAMP3+" = pal[6],
                            Monocyte = pal[7],
                            Neutrophil = pal[8],
                            Macrophage = pal[9],
                            pDC = pal[10]),
                   Group = c(Ctl = "#2b90d9", CAS = "#f26d5b"))

bk <- seq(-2,2,by=0.01)
pheatmap::pheatmap(matrix_per,scale = "row",cluster_cols = F,
                   color = colorRampPalette(c("navy","white","firebrick"))(length(bk)),
                   legend_breaks=seq(-2,2,1),breaks=bk,
                   annotation_col = anno_col,
                   annotation_colors = anno_color,
                   show_colnames = F,
                   gaps_col = seq(10,50,10),
                   border_color = "white",
                   width = 12,height = 6,
                   cellwidth = 7,cellheight = 12,
                   fontsize_row = 12,
                   treeheight_row = 20,
                   filename = "Figure6F_cytokine_heatmap_perSample.pdf")
