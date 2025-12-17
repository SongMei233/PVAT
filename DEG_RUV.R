library(ruvIIInb)
library(SingleCellExperiment)
library(scater)
library(scran)
library(scuttle)
library(edgeR)
library(SingleR)
library(celldex)
library(hrbrthemes)
library(tidyverse)
library(ggplot2)
library(uwot)
library(scMerge)
library(Seurat)
library(randomcoloR)
library(dittoSeq)
library(pheatmap)
library(gridExtra)
library(igraph)
library(DelayedArray)
library(qs)
library(IOBR)

total_ADSC <- qread("total_ADSC_mergePublic_V2.qs")

sce <- SingleCellExperiment(assays=list(counts=GetAssayData(total_ADSC, assay = "RNA", slot = "counts")),
                            colData=total_ADSC@meta.data)

sce <- addPerCellQCMetrics(x = sce,subsets=list(Mito=grep("MT-",rownames(sce))))

libsize_drop <- isOutlier(
  metric = sce$total,
  nmads = 2,
  type = "lower",
  log = TRUE)
colData(sce)$libsize_drop<-libsize_drop

mito_drop <- isOutlier(
  metric = colData(sce)$subsets_Mito_percent,
  nmads = 3,
  type = "higher")
colData(sce)$mito_drop<-mito_drop

plot_df <- data.frame(logtotal=log(sce$total),
                      libsize_drop=factor(libsize_drop),
                      mito_drop=factor(mito_drop),
                      logdetected=log(sce$detected))
plot_df$mito_drop<-relevel(plot_df$mito_drop, "TRUE")
plot_df$libsize_drop<-relevel(plot_df$libsize_drop, "TRUE")

p <- plot_df %>%
  ggplot( aes(x=logtotal, fill=libsize_drop)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity',bins=30) +
  scale_fill_manual(values=c( "#69b3a2","#404080")) +
  theme_ipsum() +
  labs(fill="") +
  xlab("Cell level log total count")+
  ylab("Frequency")
p+ggtitle('The distribution of cell-level log total counts\n  flagging cells with low library size')

p2 <- plot_df %>%
  arrange(desc(mito_drop)) %>%
  ggplot( aes(x=logtotal,y=logdetected, fill=mito_drop)) +
  geom_point(
    mapping = aes(colour = mito_drop, shape = libsize_drop),size = 2,alpha = 5 / 6) +     scale_color_manual(values=c("#69b3a2","#404080"))+
  xlab("Log total cell-level count")+ylab("Log number of cell-level detected genes")
p2+ggtitle('The cell-level number of detected genes vs total cell-level count (on log scale)')


sce <- addPerFeatureQCMetrics(x = sce)

#Remove genes with zero counts for each gene
sce <- subset(sce, rowData(sce)$mean>0 )

# detect low abundant genes
lowcount_drop <-log(rowData(sce)$mean)< -5

#mean count of genes across all cells
plot_df2 <- data.frame(mean_genecount=log(rowData(sce)$mean), lowcount_drop=factor(lowcount_drop))
p <- plot_df2 %>%
  arrange(desc(lowcount_drop)) %>%
  ggplot( aes(x=mean_genecount, fill=lowcount_drop)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity',bins=50) +
  scale_fill_manual(values=c("#404080", "#69b3a2")) +
  theme_ipsum() +
  labs(fill="") +xlab("Log mean count across all cells")+ylab("Frequency")
p+ggtitle('The distribution of log mean count across\n all cells, flagging those with a low mean count')

sce <- sce[!(lowcount_drop), !(libsize_drop | mito_drop)]

# Reading in the control genes
data(Hs.schk)
ctl <- as.character(Hs.schk)
# Creating a logical vector to idetify control genes
rowData(sce)$ctlLogical<-rownames(assays(sce)$counts) %in% ctl

# Construct the replicate matrix M using pseudo-replicates identified using initial clustering
M <- matrix(0,ncol(assays(sce)$counts),length(unique(sce$cluster_name_slim)))
colnames(M) <- unique(sce$cluster_name_slim)
cl<- unique(sce$cluster_name_slim)
for(CL in cl)
  M[which(sce$cluster_name_slim==CL),CL] <- 1


#RUV-III-NB code
ruv3nb_out <- fastruvIII.nb(Y=DelayedArray(assays(sce)$counts), # count matrix with genes as rows and cells as columns
                          M=M, #Replicate matrix constructed as above
                          ctl=rowData(sce)$ctlLogical, #A vector denoting control genes
                          k=2, # dimension of unwanted variation factors
                          ncores = 12,
                          batch = as.numeric(colData(sce)$cohort))

qsave(ruv3nb_out,file="ruv3nb_out.qs")

data <- assays(ruv3nb_out)$logPAC
total_ADSC <- CreateSeuratObject(counts = as.sparse(data),
                                 meta.data = as.data.frame(colData(sce)))

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
  
  for(y in c("PVAT","BAT","SAT","VAT")){
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
    
    colnames(exp) <- c("BAT","Inhouse","PVAT","SAT","VAT")
    
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

res <- data.frame(V1=names(diff_gene),
                  Num = lapply(diff_gene,length) %>% unlist())
write.csv(res,file = "DEG_number_in each group.csv",row.names = F,col.names = T)

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
                      fdrThr = 1,
                      isOutput = F)
  GOBP <- GOBP[GOBP$FDR<0.05&GOBP$overlap>2,]
  
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


meta <- total_ADSC@meta.data[,c("sample","cohort","group")] %>% distinct()
meta$sample <- paste(meta$sample,meta$group,sep = "_")

for(x in c("ADSC_CXCL14","ADSC_CD55")){
  
  pb2 <- qread(paste0("PseudoBulk_",x,".qs"))
  tpm <- count2tpm(pb2,idType = "Symbol")
  
  for(y in c("PVAT","SAT","VAT","BAT")){

    gene_diff <- c(diff_gene[[paste(x,y,"up",sep = "_")]],diff_gene[[paste(x,y,"down",sep = "_")]])
    
    meta_sub <- filter(meta,group %in% c("InHouse",y))
    
    tpm_plot <- tpm[gene_diff,intersect(meta_sub$sample,colnames(tpm))]
    
    tpm_plot <- log2(tpm_plot+1)
    
    meta_sub <- meta_sub[meta_sub$sample %in% intersect(meta_sub$sample,colnames(tpm)),]
    
    anno_col <- data.frame(group=meta_sub$group,
                           row.names = colnames(tpm_plot))
    
    anno_col$group <- factor(anno_col$group,levels = c("InHouse",y))
    
    pheatmap::pheatmap(na.omit(tpm_plot),cluster_rows = T,cluster_cols = F,
                       scale = "row",
                       gaps_col = c(8),
                       annotation_col = anno_col,
                       show_rownames = F,
                       show_colnames = F,
                       width = 6,height = 6,
                       filename = paste0("Heatmap_DEG_",x,"_",y,".pdf"))
    
  }
}