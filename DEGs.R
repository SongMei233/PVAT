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
library(SingleCellExperiment)
library(scuttle)
library(Matrix.utils)
library(IOBR)
library(DESeq2)
library(sva)
library(UpSetR)
library(limma)

total_ADSC <- qread("total_ADSC_mergePublic_V2.qs")

total_ADSC$sample <- paste(total_ADSC$sample,total_ADSC$group,sep = "_")

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

meta <- filter(meta,cohort != "Nature2025")

pb2 <- pseudoBulk(total_ADSC_CD55)

qsave(pb2,file = "PseudoBulk_ADSC_CXCL14.qs")

pb2 <- qread("PseudoBulk_ADSC_CD55.qs")

for(x in c("PVAT","SAT","VAT")){
  
  meta_sub <- filter(meta,group %in% c("InHouse",x))
  pb2_sub <- pb2[,row.names(meta_sub)]
  
  meta_sub$group[meta_sub$group!="InHouse"] <- "Control"
  meta_sub$cohort <- as.character(meta_sub$cohort)
  
  pb2_sub <- pb2_sub[row.names(pb2_sub) %in% proCod,]

  data.features <- row.names(pb2_sub)
  RB_genes = c(grep("^RPS", data.features, v=T),grep("^RPL", data.features, v=T))
  MT_genes = grep("^MT-", data.features, v=T)
  remove <- c(RB_genes,MT_genes)
  pb2_sub <- pb2_sub[!(row.names(pb2_sub) %in% remove),]

  pb2_sub <- pb2_sub[rowSums(pb2_sub) > 0, ]
  
  # index1 <- meta_sub[meta_sub$group=="InHouse","sample"]
  # index2 <- meta_sub[meta_sub$group=="Control","sample"]
  # 
  # pb2_sub1 <- pb2[,index1]
  # pb2_sub2 <- pb2[,index2]
  # 
  # gene1 <- pb2_sub1[rowSums(pb2_sub1 > 0) > ncol(pb2_sub1)/2, ] %>% row.names()
  # gene2 <- pb2_sub2[rowSums(pb2_sub2 > 0) > ncol(pb2_sub2)/2, ] %>% row.names()
  # 
  # gene_sel <- intersect(gene1,gene2)
  
  # pb2_sub <- pb2_sub[row.names(pb2_sub) %in% gene_sel,]
  
  tpm <- count2tpm(pb2_sub,idType = "Symbol")
  
  combat_corrected_data <- ComBat(
    dat = tpm,
    batch = meta_sub$cohort
  )
  
  design <- model.matrix(~0 + group, data = meta_sub)
  colnames(design) <- gsub("group", "", colnames(design))
  
  fit <- lmFit(combat_corrected_data, design)
  
  cont.matrix <- makeContrasts(
    Case_vs_Control = InHouse - Control,
    levels = design
  )
  
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2, trend = TRUE)
  results <- topTable(fit2, coef = 1, adjust.method = "BH", number = Inf)

  write.csv(results,file = paste0("DEG_ADSC_CD55_",x,".csv"),row.names = T,col.names = T)
  
}


i=1
diff_gene <- list()
for(x in c("ADSC_CXCL14","ADSC_CD55")){
  
  for(y in c("PVAT","SAT","VAT")){
    marker <- read.csv(paste0("DEG_",x,"_",y,".csv"),row.names = 1,header = T)
    
    gene_up <- filter(marker,P.Value < 0.05,logFC > 1) %>% row.names()
    gene_down <- filter(marker,P.Value < 0.05,logFC < (-1)) %>% row.names()
    
    diff_gene[[i]] <- gene_up
    diff_gene[[i+1]] <- gene_down
    names(diff_gene)[i] <- paste(x,y,"up",sep = "_")
    names(diff_gene)[i+1] <- paste(x,y,"down",sep = "_")
    i=i+2
    
    data_vol <- subset(marker,select=c(4,1))
    colnames(data_vol) <- c("pvalue","log2FoldChange")
    data_vol <- na.omit(data_vol)
    data_vol$sig[((data_vol$log2FoldChange<1)&(data_vol$log2FoldChange>-1))|(data_vol$pvalue>0.05)] <- "no"
    data_vol$sig[(data_vol$log2FoldChange>=1)&(data_vol$pvalue<=0.05)] <- "up"
    data_vol$sig[(data_vol$log2FoldChange<=-1)&(data_vol$pvalue<=0.05)] <- "down"
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
      #geom_text_repel(aes(log2FoldChange,-1*log10(pvalue),label=label),size=3)+
      labs(x="log2(FoldChange)",y="-log10(FDR)") + 
      scale_color_manual(values = cols)+
      geom_hline(yintercept=-log10(0.05),linetype=4)+
      geom_vline(xintercept=c(-1,1),linetype=4)+
      theme_bw()+
      xlim(-10,10)+
      theme(panel.grid.major =element_blank(),
            panel.grid.minor =element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(color="black"),
            strip.background.x = element_rect(colour = "white"),
            axis.text=element_text(size=20),
            axis.title=element_text(size=20),
            legend.position = "none")
    
    ggsave(filename = paste0("./Volcano_",x,"_",y,".pdf"),p,device = "pdf",width = 6,height = 6)

  }
}

res <- data.frame(V1=names(diff_gene),
                  Num = lapply(diff_gene,length) %>% unlist())
write.csv(res,file = "DEG_number_in each group.csv",row.names = F,col.names = T)


library(WebGestaltR)
library(Hmisc)

for(i in 9:length(diff_gene)){
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