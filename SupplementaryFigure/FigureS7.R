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
  cols <- tryCatch(ncol(x), error = function(e) NULL)
  ifelse(!is.null(cols) && cols >= 40, TRUE, FALSE)
})

data.list <- data.list[keep_elements]

data.list <- lapply(data.list,FUN = function(x){DefaultAssay(x) <- "RNA";x <- SCTransform(x)})

options(future.globals.maxSize= 100000000000000)

data.features <- SelectIntegrationFeatures(object.list = data.list)

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

