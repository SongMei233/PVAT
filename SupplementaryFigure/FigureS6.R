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
ggsave("FigureS6A_VlnPlot_remodeling_secreted.pdf",width = 12,height = 9)