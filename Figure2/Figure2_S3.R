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

######################Figure 2 ADSC##############################
total <- readRDS("total_final.2.rds")

total_ADSC <- subset(total,subset = cluster_name_slim == "ADSC")


#re-integrate
data.list <- SplitObject(total_ADSC,split.by = "sample")

#remove samples with ADSCs less than 50

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

#Check if there are any cluster that expressed feature of other cell lineages

markers <- c("PTPRC","MS4A1","CD19","BANK1","JCHAIN","SDC1","IGKC","CD3D","CD4","CCR7","SELL","CD8A","GZMK","CCL5","KLRD1","GNLY","NKG7","MKI67","STMN1","TUBA1B","IL1B","S100A8","S100A9","CD68","C1QA","LYVE1","CD1E","LAMP3","FSCN1","PTGDS","IL3RA","LILRA4","F3","CXCL14","CD55")
Idents(total_ADSC) <- "integrated_snn_res.1.5"
DefaultAssay(total_ADSC) <- "SCT"
DotPlot(total_ADSC, features = markers, dot.scale = 6)  + theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) + RotatedAxis()

#extract pure ADSC

total_ADSC <- subset(total_ADSC,subset = integrated_snn_res.1.5 %in% c(0:8,10,12,16))
total_ADSC <- Reduct(total_ADSC)

DimPlot(total_ADSC, reduction = "umap", label=T,group.by = "integrated_snn_res.1.5")

total_ADSC$cluster_name_slim <- "ADSC_CXCL14"
total_ADSC$cluster_name_slim[total_ADSC$integrated_snn_res.1.5 %in% c(1,5,6,8,10)] <- "ADSC_CD55"
total_ADSC$cluster_name_slim <- factor(total_ADSC$cluster_name_slim,
                                       levels = c("ADSC_CXCL14","ADSC_CD55"))

DimPlot(total_ADSC, reduction = "umap", label=T,group.by = "cluster_name_slim")

Idents(total_ADSC) <- "cluster_name_slim"

######################Figure 2##############################

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
ggsave("FigureS3B-D_feature_ADSC.pdf",p,width = 12,height = 4)

p <- FeaturePlot(total_ADSC,features = c("PDGFRA","PDGFRB"),order = T,cols = c("lightgrey","firebrick"),ncol = 2);p
ggsave("FeatureS3E_ADSC_PDGFR.pdf",p,width = 8,height = 4)

p <- FeaturePlot(total_ADSC,features = c("GPC3"),order = T,cols = c("lightgrey","firebrick"),ncol = 1);p
ggsave("FeatureS3H_ADSC_GPC3.pdf",p,width = 4,height = 4)

saveRDS(total_ADSC,file = "total_ADSC.rds")

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


######################Figure 2E-F ADSC's function##############################

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

ggsave("Figure2E_GeneEnrichment_ADSC_CD55.pdf",p,width = 10,height = 4)

#the enrichment method was the same for ADSC_CXCL14 in Figure 2D

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
ggsave(paste0("Figure2F_",features[2],".pdf"),width = 4,height = 5)

######################Figure S3J monocle##############################

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

pdf("FigureS3J_Pseudotime_cell_name.pdf",width = 5,height = 5)
plot_cell_trajectory(cds, color_by = "cluster_name_slim",show_branch_points = F,cell_size = 1)+
  scale_color_manual(values = pal[c(1,3)])
dev.off()


saveRDS(cds,file = "monocle_ADSC.rds")

######################Figure S3K cytoTrace##############################

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

######################Figure S3L highly-expressed genes####################

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
ggsave("FigureS3L_Genes_dotplot.pdf",width = 7,height = 4)

######################Figure S3N cell proportion box plot####################

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

ggsave("FigureS3N_Cell proportion between Ctl and CAS_ADSC.pdf",p,width = 5,height = 3)
