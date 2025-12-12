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

markers <- c("NCAM1",'FCGR3A',"CCR7","GZMK","CCL5","IL7R")
p <- VlnPlot(total_T,features = markers,pt.size = 0,ncol = 3)
ggsave("FigureS9A_VlnPlot_NK.pdf",p,width = 8,height = 8)

#############################cell proportion box plot of NK cell################
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

ggsave("FigureS9B_NK_Cell proportion between Ctl and CAS.pdf",p,width = 6,height = 3)

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


#################cell proportion box plot of Macrophage##################
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

ggsave("FigureS9G_Macrophage_Cell proportion between Ctl and CAS.pdf",p,width = 5,height = 5)

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
ggsave("FigureS9H_Macrophage_DEG_network_up_20.pdf",p,width = 7,height = 6)
