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
library(CellMixS)
library(SingleCellExperiment)

total_ADSC <- qread("total_ADSC_mergePublic.qs")

sce_ADSC <- as.SingleCellExperiment(total_ADSC,assay = "RNA")

sce_cms <- cms(sce_ADSC, k = 70, group = "cohort")
sce_cms <- cms(sce_ADSC, k = 70, group = "cohort",dim_red = "UMAP",res_name = "CCA",n_dim = 2)
sce_ldf <- ldfDiff(sce_pre_list, sce_combined, group = "batch", k = 70)

visHist(sce_cms)

visHist(sce_cms, metric = "cms.",  n_col = 2)

visMetric(sce_cms, metric_var = "cms_smooth.CCA",dim_red = "UMAP")
visGroup(sce_cms, group = "cohort",dim_red = "UMAP")

##############kBET##############

library(kBET)

#10X-BGI

list1 <- paste0("PVAT",1:6)
list2 <- c(paste0("CAS",5:7),"CBT4")
total$batch <- "10X"
total$batch[total$orig.ident %in% list2] <- "BGI"
total$batch <- factor(total$batch,levels = c("10X","BGI"))

pca_mat <- Embeddings(total, reduction = "pca")[, 1:30]  
batch <- total$batch

batch.estimate <- kBET(pca_mat, batch,do.pca = FALSE,k0 = 10)

plot.data <- data.frame(class=rep(c('observed', 'expected'), 
                                  each=length(batch.estimate$stats$kBET.observed)), 
                        data =  c(batch.estimate$stats$kBET.observed,
                                  batch.estimate$stats$kBET.expected))
g <- ggplot(plot.data, aes(class, data)) + geom_boxplot() + 
  labs(x='Test', y='Rejection rate',title='kBET test results') +
  theme_bw() +  
  scale_y_continuous(limits=c(0,1));g


pvals <- batch.estimate$results$kBET.pvalue.test
pvals <- p.adjust(pvals,method = "BH")
reject <- ifelse(pvals < 0.05, "Unmixed", "Mixed-well")
total$kBET <- reject

p <- DimPlot(total, group.by = "kBET",cols = c("steelblue", "firebrick"),raster = F,pt.size = 0.1);p
ggsave("FigureS1_Batch_evaluate.pdf",p,width = 8,height = 6)


#ADSC-Public

total <- qread("total_ADSC_mergePublic.qs")

pca_mat <- Embeddings(total, reduction = "pca")[, 1:30]  
batch <- total$cohort

batch.estimate <- kBET(pca_mat, batch,do.pca = FALSE,k0 = 10)

pvals <- batch.estimate$results$kBET.pvalue.test
pvals <- p.adjust(pvals,method = "BH")
reject <- ifelse(pvals < 0.05, "Unmixed", "Mixed-well")
total$kBET <- reject

p <- DimPlot(total, group.by = "kBET",cols = c("steelblue", "firebrick"),raster = F,pt.size = 0.1);p
ggsave("FigureS7_Batch_evaluate_ADSC.pdf",p,width = 8,height = 6)
