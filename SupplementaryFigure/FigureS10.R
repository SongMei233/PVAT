####################Integrate public and in-house data (T cell)###########################

total_T <- readRDS("total_T.rds")

DimPlot(total_T)

DefaultAssay(total_T) <- "RNA"
total_T[["integrated"]] <- NULL

total_T$cohort <- "InHouse"

data.list1 <- SplitObject(total_T,split.by = "sample")

#public data

#ATVB_2024

data_AVTB2024 <- qread("../PublicSCData/ATVB_2024/total_T.qs")
data_AVTB2024$cohort <- "AVTB2024"
data.list2 <- SplitObject(data_AVTB2024,split.by = "sample")

#Nature_2020

data_Nature2020 <- qread("../PublicSCData/Nature_2020/total_T.qs")
data_Nature2020$cohort <- paste0("Nature2020_",data_Nature2020$group)
data.list3 <- SplitObject(data_Nature2020,split.by = "sample")

#NM_2020

data_NM2020 <- qread("../PublicSCData/NatureMeta_2020/total_T.qs")
data_NM2020$SampleName <- droplevels(data_NM2020$SampleName)
data_NM2020$sample <- data_NM2020$SampleName
data_NM2020$cohort <- paste0("NM2020_",data_NM2020$Tissue)
data.list4 <- SplitObject(data_NM2020,split.by = "sample")

#NM_2021

data_NM2021 <- qread("../PublicSCData/NatureMeta_2021/total_T.qs")
data_NM2021$cohort <- "NM2021"
data.list5 <- SplitObject(data_NM2021,split.by = "sample")


data.list <- c(data.list1,data.list2,data.list3,data.list4,data.list5)

keep_elements <- sapply(data.list, function(x) {
  # Check if the element has a number of columns attribute (such as matrices, data frames)
  cols <- tryCatch(ncol(x), error = function(e) NULL)
  ifelse(!is.null(cols) && cols >= 40, TRUE, FALSE)
})

data.list <- data.list[keep_elements]

rm(list = setdiff(ls(),"data.list"))

data.list <- lapply(data.list,FUN = function(x){DefaultAssay(x) <- "RNA";x <- SCTransform(x)})

options(future.globals.maxSize= 100000000000000)

data.features <- SelectIntegrationFeatures(object.list = data.list)
data.list <- PrepSCTIntegration(object.list = data.list,assay = "SCT", anchor.features = data.features,verbose = FALSE)
data.anchors <- FindIntegrationAnchors(object.list = data.list, normalization.method = "SCT", anchor.features = data.features,k.filter = 40, verbose = FALSE)
total_T <- IntegrateData(anchorset = data.anchors, normalization.method = "SCT",k.weight = 40, verbose = FALSE)

total_T <- Reduct(total_T)

Idents(total_T) <- "integrated_snn_res.0.5"
DimPlot(total_T,label = T)
DimPlot(total_T,label = T,group.by = "cluster_name_slim")

qsave(total_T,file = "total_T_mergePublic.qs")

#harmony

DefaultAssay(total_T) <- "RNA"

total_T <- NormalizeData(total_T)
total_T <- FindVariableFeatures(total_T,selection.method = "vst", nfeatures = 1000)
total_T <- ScaleData(total_T)
total_T <- RunPCA(total_T, verbose = FALSE)

total_T@meta.data$sample <- as.factor(total_T@meta.data$sample)
total_T <- RunHarmony(total_T,group.by.vars="sample",plot_convergence = TRUE)
ElbowPlot(total_T,ndims = 50)
total_T <- RunUMAP(total_T,reduction = "harmony", dims = 1:30)
total_T <- FindNeighbors(total_T,reduction = "harmony", dims = 1:30)
total_T <- FindClusters(total_T, resolution = 0.5)

DimPlot(total_T, reduction = "umap", label=T,group.by = "RNA_snn_res.0.5")
DimPlot(total_T,label = T,group.by = "cluster_name_slim")

#FastMNN
library(SeuratWrappers)

for(i in 1:length(data.list)){
  DefaultAssay(data.list[[i]]) <- "RNA"
  data.list[[i]][["SCT"]] <- NULL
  data.list[[i]] <- NormalizeData(data.list[[i]])
  data.list[[i]] <- FindVariableFeatures(data.list[[i]])
}

total_T <- RunFastMNN(object.list = data.list)
total_T <- RunUMAP(total_T, reduction = "mnn", dims = 1:30)
total_T <- FindNeighbors(total_T, reduction = "mnn", dims = 1:30) 
total_T <- FindClusters(total_T,resolution = 0.5)

DimPlot(total_T, reduction = "umap", label=T,group.by = "RNA_snn_res.0.5")
DimPlot(total_T,label = T,group.by = "cluster_name_slim")

clu13.marker <- FindMarkers(total_T, ident.1 = "13", logfc.threshold = 0.2, min.pct = 0.5, min.diff.pct = 0.2,only.pos = T)

Plot_Density_Custom(total_T,features = c("PDCD1","HSPA1A"))
Plot_Density_Custom(total_T,features = c("NCAM1","FCGR3A"))

total_T$cluster_name_slim <- "CD4/CD8Tmemory"
total_T$cluster_name_slim[total_T$RNA_snn_res.0.5 %in% c(5,9)] <- "CD4Treg"
total_T$cluster_name_slim[total_T$RNA_snn_res.0.5 == 4] <- "CD8Teffector"
total_T$cluster_name_slim[total_T$RNA_snn_res.0.5 %in% c(3,6)] <- "Remove"
total_T$cluster_name_slim[total_T$RNA_snn_res.0.5 %in% c(13)] <- "NK"

total_NK <- subset(total_T,subset = cluster_name_slim == "NK")

total_T <- subset(total_T,subset = cluster_name_slim %in% c("CD4/CD8Tmemory","CD4Treg","CD8Teffector"))

total_T$cluster_name_slim <- factor(total_T$cluster_name_slim, levels = c("CD4/CD8Tmemory","CD4Treg","CD8Teffector"))

total_T <- subset(total_T,subset = cohort %in% c("InHouse","Nature2020_BAT","AVTB2024","NM2020_SAT"))

total_T$cohort <- factor(total_T$cohort,levels = c("InHouse","Nature2020_BAT","AVTB2024",
                                                   "NM2020_SAT"))

Idents(total_T) <- "cluster_name_slim"
pdf(file = "Dimplot_T_mergePublic.pdf",width = 15,height = 5)
DimPlot(total_T,reduction = "umap",label = F,pt.size = 0.15,split.by = "cohort")+theme(legend.position = "bottom")
dev.off()

#cell proportion

cell.prop <- table(total_T$cluster_name_slim,total_T$cohort) %>% as.data.frame.matrix()
#cell.prop <- t(t(cell.prop)/colSums(cell.prop)) %>% as.data.frame()
cell.prop$cell <- row.names(cell.prop)
cell.prop <- pivot_longer(cell.prop,cols = setdiff(colnames(cell.prop),"cell"),names_to = "cohort",values_to = "num")
cell.prop = ddply(cell.prop, 'cohort', transform, percent=num/sum(num)*100)

cell.prop$cell <- factor(cell.prop$cell,levels = rev(levels(total_T$cluster_name_slim)))
cell.prop$cohort <- factor(cell.prop$cohort,levels = levels(total_T$cohort))

p <- ggplot(cell.prop,aes(x=cohort,y=percent,fill=cell))+
  geom_bar(stat = 'identity',width = 0.8,colour='black')+
  scale_fill_manual(values = rev(hue_pal()(3)))+
  scale_y_continuous(expand = c(0,0))+
  theme_classic()+
  theme(axis.text = element_text(size = 14,colour = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 16,colour = "black"))+
  ylab("Cell ratio (%)")+xlab("");p
ggsave("Cell_proportion_bar_T_mergePublic.pdf",p,width = 5,height = 3.5)

qsave(total_T,file = "total_T_mergePublic.qs")

qsave(total_NK,file = "total_NK_Public.qs")

####################Integrate public and in-house data (Myeloid cell)###########################

total_Myeloid <- readRDS("total_myeloid.rds")

Idents(total_Myeloid) <- "cluster_name_slim"

DimPlot(total_Myeloid)

DefaultAssay(total_Myeloid) <- "RNA"
total_Myeloid[["integrated"]] <- NULL

total_Myeloid$cohort <- "InHouse"

data.list1 <- SplitObject(total_Myeloid,split.by = "sample")

#public data

#ATVB_2024

data_AVTB2024 <- qread("./PublicSCData/ATVB_2024/total_Mye.qs")
data_AVTB2024$cohort <- "AVTB2024"
data.list2 <- SplitObject(data_AVTB2024,split.by = "sample")

#Nature_2020

data_Nature2020 <- qread("./PublicSCData/Nature_2020/total_Mye.qs")
data_Nature2020$cohort <- paste0("Nature2020_",data_Nature2020$group)
data.list3 <- SplitObject(data_Nature2020,split.by = "sample")

#NM_2020

data_NM2020 <- qread("./PublicSCData/NatureMeta_2020/total_Mye.qs")
data_NM2020$SampleName <- droplevels(data_NM2020$SampleName)
data_NM2020$sample <- data_NM2020$SampleName
data_NM2020$cohort <- paste0("NM2020_",data_NM2020$Tissue)
data.list4 <- SplitObject(data_NM2020,split.by = "sample")

#NM_2021

data_NM2021 <- qread("./PublicSCData/NatureMeta_2021/total_Mye.qs")
data_NM2021$cohort <- "NM2021"
data.list5 <- SplitObject(data_NM2021,split.by = "sample")


data.list <- c(data.list1,data.list2,data.list3,data.list4,data.list5)

keep_elements <- sapply(data.list, function(x) {
  cols <- tryCatch(ncol(x), error = function(e) NULL)
  ifelse(!is.null(cols) && cols >= 40, TRUE, FALSE)
})

data.list <- data.list[keep_elements]

rm(list = setdiff(ls(),"data.list"))

data.list <- lapply(data.list,FUN = function(x){DefaultAssay(x) <- "RNA";x <- SCTransform(x)})

options(future.globals.maxSize= 100000000000000)

data.features <- SelectIntegrationFeatures(object.list = data.list)
data.list <- PrepSCTIntegration(object.list = data.list,assay = "SCT", anchor.features = data.features,verbose = FALSE)
data.anchors <- FindIntegrationAnchors(object.list = data.list, normalization.method = "SCT", anchor.features = data.features,k.filter = 40, verbose = FALSE)
total_Myeloid <- IntegrateData(anchorset = data.anchors, normalization.method = "SCT",k.weight = 40, verbose = FALSE)

total_Myeloid <- Reduct(total_Myeloid)

Idents(total_Myeloid) <- "integrated_snn_res.0.5"
DimPlot(total_Myeloid,label = T)
DimPlot(total_Myeloid,label = T,group.by = "cluster_name_slim")

markers <- c("NCAM1",'FCGR3A',"CD69","NKG7","XCR1","CLEC9A","CD1C","CLEC10A","LAMP3","CCL19","CDKN1C",'LILRB2',"FCN1","S100A8","S100A9","CSF3R","CD163","C1QA","MRC1","IL3RA","LILRA4","TCF4")
Idents(total_Myeloid) <- "integrated_snn_res.0.5"
DefaultAssay(total_Myeloid) <- "SCT"
DotPlot(total_Myeloid, features = markers)  + theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) + RotatedAxis()

#FastMNN
library(SeuratWrappers)

count <- total_Myeloid@assays$RNA@counts
meta <- total_Myeloid@meta.data
total_Myeloid <- CreateSeuratObject(counts = count,meta.data = meta)

data.list <- SplitObject(total_Myeloid,split.by = "sample")

for(i in 1:length(data.list)){
  DefaultAssay(data.list[[i]]) <- "RNA"
  data.list[[i]] <- NormalizeData(data.list[[i]])
  data.list[[i]] <- FindVariableFeatures(data.list[[i]])
}

total_Myeloid <- RunFastMNN(object.list = data.list)
total_Myeloid <- RunUMAP(total_Myeloid, reduction = "mnn", dims = 1:30)
total_Myeloid <- FindNeighbors(total_Myeloid, reduction = "mnn", dims = 1:30) 
total_Myeloid <- FindClusters(total_Myeloid,resolution = 0.5)

DimPlot(total_Myeloid, reduction = "umap", label=T,group.by = "RNA_snn_res.0.5")
DimPlot(total_Myeloid,label = T,group.by = "cluster_name_slim")

markers <- c("PTPRC","NCAM1","FCGR3A","ITGA1","B3GAT1","KLRD1","GNLY","NKG7","IL1B","S100A8","S100A9","CSF3R","CD68","C1QA","LYVE1","XCR1","CLEC9A","CD1C","CD1E","LAMP3","FSCN1","PTGDS","IL3RA","LILRA4")
Idents(total_Myeloid) <- "RNA_snn_res.0.5"
DotPlot(total_Myeloid, features = markers, dot.scale = 6)  + theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) + RotatedAxis()

clu12.marker <- FindMarkers(total_Myeloid, ident.1 = "12", logfc.threshold = 0.2, min.pct = 0.5, min.diff.pct = 0.2,only.pos = T)

total_Myeloid$cluster_name_slim <- "Macrophage"
total_Myeloid$cluster_name_slim[total_Myeloid$RNA_snn_res.0.5 %in% c(4)] <- "Neutrophil/Monocyte"
total_Myeloid$cluster_name_slim[total_Myeloid$RNA_snn_res.0.5 == 8] <- "cDC2"
total_Myeloid$cluster_name_slim[total_Myeloid$RNA_snn_res.0.5 %in% c(10)] <- "pDC"
total_Myeloid$cluster_name_slim[total_Myeloid$RNA_snn_res.0.5 %in% c(11)] <- "cDC1"
total_Myeloid$cluster_name_slim[total_Myeloid$RNA_snn_res.0.5 %in% c(13)] <- "DC_LAMP3+"
total_Myeloid$cluster_name_slim[total_Myeloid$RNA_snn_res.0.5 %in% c(7,12)] <- "Remove"

total_Myeloid <- subset(total_Myeloid,subset = cluster_name_slim != "Remove")

total_Myeloid$cluster_name_slim <- factor(total_Myeloid$cluster_name_slim,
                                          levels = c("Macrophage","Neutrophil/Monocyte","cDC1","cDC2","DC_LAMP3+","pDC"))


total_Myeloid$cohort <- factor(total_Myeloid$cohort,levels = c("InHouse","Nature2020_BAT","NM2021","AVTB2024",
                                                               "NM2020_SAT",
                                                               "NM2020_VAT"))

Idents(total_Myeloid) <- "cluster_name_slim"
pdf(file = "Dimplot_Myeloid_mergePublic.pdf",width = 16,height = 5)
DimPlot(total_Myeloid,reduction = "umap",label = F,pt.size = 0.15,split.by = "cohort")+theme(legend.position = "bottom")
dev.off()

#cell proportion

cell.prop <- table(total_Myeloid$cluster_name_slim,total_Myeloid$cohort) %>% as.data.frame.matrix()
#cell.prop <- t(t(cell.prop)/colSums(cell.prop)) %>% as.data.frame()
cell.prop$cell <- row.names(cell.prop)
cell.prop <- pivot_longer(cell.prop,cols = setdiff(colnames(cell.prop),"cell"),names_to = "cohort",values_to = "num")
cell.prop = ddply(cell.prop, 'cohort', transform, percent=num/sum(num)*100)

cell.prop$cell <- factor(cell.prop$cell,levels = rev(levels(total_Myeloid$cluster_name_slim)))
cell.prop$cohort <- factor(cell.prop$cohort,levels = levels(total_Myeloid$cohort))

p <- ggplot(cell.prop,aes(x=cohort,y=percent,fill=cell))+
  geom_bar(stat = 'identity',width = 0.8,colour='black')+
  scale_fill_manual(values = rev(hue_pal()(6)))+
  scale_y_continuous(expand = c(0,0))+
  theme_classic()+
  theme(axis.text = element_text(size = 14,colour = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 16,colour = "black"))+
  ylab("Cell ratio (%)")+xlab("");p
ggsave("Cell_proportion_bar_Myeloid_mergePublic.pdf",p,width = 6,height = 3.5)

qsave(total_Myeloid,file = "total_Myeloid_mergePublic.qs")

####################Integrate public and in-house data (B cell)###########################

total_B <- readRDS("total_B.rds")

Idents(total_B) <- "cluster_name_slim"

DimPlot(total_B)

DefaultAssay(total_B) <- "RNA"
total_B[["integrated"]] <- NULL

total_B$cohort <- "InHouse"

data.list1 <- SplitObject(total_B,split.by = "sample")

#public data

#ATVB_2024

data_AVTB2024 <- qread("./PublicSCData/ATVB_2024/total_B.qs")
data_AVTB2024$cohort <- "AVTB2024"
data.list2 <- SplitObject(data_AVTB2024,split.by = "sample")

#Nature_2020

data_Nature2020 <- qread("./PublicSCData/Nature_2020/total_B.qs")
data_Nature2020$cohort <- paste0("Nature2020_",data_Nature2020$group)
data.list3 <- SplitObject(data_Nature2020,split.by = "sample")

#NM_2020

data_NM2020 <- qread("./PublicSCData/NatureMeta_2020/total_B.qs")
data_NM2020$SampleName <- droplevels(data_NM2020$SampleName)
data_NM2020$sample <- data_NM2020$SampleName
data_NM2020$cohort <- paste0("NM2020_",data_NM2020$Tissue)
data.list4 <- SplitObject(data_NM2020,split.by = "sample")

data.list <- c(data.list1,data.list2,data.list3,data.list4)

keep_elements <- sapply(data.list, function(x) {
  cols <- tryCatch(ncol(x), error = function(e) NULL)
  ifelse(!is.null(cols) && cols >= 40, TRUE, FALSE)
})

data.list <- data.list[keep_elements]

rm(list = setdiff(ls(),"data.list"))

data.list <- lapply(data.list,FUN = function(x){DefaultAssay(x) <- "RNA";x <- SCTransform(x)})

options(future.globals.maxSize= 100000000000000)

data.features <- SelectIntegrationFeatures(object.list = data.list)
data.list <- PrepSCTIntegration(object.list = data.list,assay = "SCT", anchor.features = data.features,verbose = FALSE)
data.anchors <- FindIntegrationAnchors(object.list = data.list, normalization.method = "SCT", anchor.features = data.features,k.filter = 40, verbose = FALSE)
total_B <- IntegrateData(anchorset = data.anchors, normalization.method = "SCT",k.weight = 40, verbose = FALSE)

total_B <- Reduct(total_B)

Idents(total_B) <- "integrated_snn_res.0.5"
DimPlot(total_B,label = T)
DimPlot(total_B,label = T,group.by = "cluster_name_slim")

qsave(total_B,file = "total_B_mergePublic.qs")

#FastMNN
library(SeuratWrappers)

count <- total_B@assays$RNA@counts
meta <- total_B@meta.data
total_B <- CreateSeuratObject(counts = count,meta.data = meta)

data.list <- SplitObject(total_B,split.by = "sample")

for(i in 1:length(data.list)){
  DefaultAssay(data.list[[i]]) <- "RNA"
  data.list[[i]] <- NormalizeData(data.list[[i]])
  data.list[[i]] <- FindVariableFeatures(data.list[[i]])
}

total_B <- RunFastMNN(object.list = data.list)
total_B <- RunUMAP(total_B, reduction = "mnn", dims = 1:30)
total_B <- FindNeighbors(total_B, reduction = "mnn", dims = 1:30) 
total_B <- FindClusters(total_B,resolution = 0.5)

DimPlot(total_B, reduction = "umap", label=T,group.by = "RNA_snn_res.0.5")
DimPlot(total_B,label = T,group.by = "cluster_name_slim")

total_B$cluster_name_slim <- "B_naive"
total_B$cluster_name_slim[total_B$RNA_snn_res.0.5 %in% c(3,4,5,8)] <- "B_memory(unswitched)"
total_B$cluster_name_slim[total_B$RNA_snn_res.0.5 %in% c(1,2,7)] <- "B_memory(switched)"
total_B$cluster_name_slim[total_B$RNA_snn_res.0.5 %in% c(9,10,11,12)] <- "Plasma cell"

total_B$cluster_name_slim <- factor(total_B$cluster_name_slim,
                                    levels = c("B_naive","B_memory(unswitched)","B_memory(switched)","Plasma cell"))
Idents(total_B) <- "cluster_name_slim"

total_B$cohort <- factor(total_B$cohort,levels = c("InHouse","Nature2020_BAT","AVTB2024",
                                                   "NM2020_SAT"))

Idents(total_B) <- "cluster_name_slim"
pdf(file = "Dimplot_B_mergePublic.pdf",width = 15,height = 5)
DimPlot(total_B,reduction = "umap",label = F,pt.size = 0.15,split.by = "cohort",cols = pal[1:10])+theme(legend.position = "bottom")
dev.off()

#cell proportion

cell.prop <- table(total_B$cluster_name_slim,total_B$cohort) %>% as.data.frame.matrix()
#cell.prop <- t(t(cell.prop)/colSums(cell.prop)) %>% as.data.frame()
cell.prop$cell <- row.names(cell.prop)
cell.prop <- pivot_longer(cell.prop,cols = setdiff(colnames(cell.prop),"cell"),names_to = "cohort",values_to = "num")
cell.prop = ddply(cell.prop, 'cohort', transform, percent=num/sum(num)*100)

cell.prop$cell <- factor(cell.prop$cell,levels = rev(levels(total_B$cluster_name_slim)))
cell.prop$cohort <- factor(cell.prop$cohort,levels = levels(total_B$cohort))

p <- ggplot(cell.prop,aes(x=cohort,y=percent,fill=cell))+
  geom_bar(stat = 'identity',width = 0.8,colour='black')+
  scale_fill_manual(values = rev(pal[1:4]))+
  scale_y_continuous(expand = c(0,0))+
  theme_classic()+
  theme(axis.text = element_text(size = 14,colour = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 16,colour = "black"))+
  ylab("Cell ratio (%)")+xlab("");p
ggsave("Cell_proportion_bar_B_mergePublic.pdf",p,width = 6,height = 3.5)

qsave(total_B,file = "total_B_mergePublic.qs")

####################Integrate public and in-house data (NK cell)###########################

total_NK <- readRDS("total_NK.rds")

Idents(total_NK) <- "cluster_name_slim"

DimPlot(total_NK)

DefaultAssay(total_NK) <- "RNA"
total_NK[["integrated"]] <- NULL

total_NK$cohort <- "InHouse"

data.list1 <- SplitObject(total_NK,split.by = "sample")

data_public <- qread("total_NK_Public.qs")

data.list2 <- SplitObject(data_public,split.by = "sample")

data.list <- c(data.list1,data.list2)

keep_elements <- sapply(data.list, function(x) {
  cols <- tryCatch(ncol(x), error = function(e) NULL)
  ifelse(!is.null(cols) && cols >= 40, TRUE, FALSE)
})

data.list <- data.list[keep_elements]

for(i in 1:length(data.list)){
  DefaultAssay(data.list[[i]]) <- "RNA"
  data.list[[i]] <- NormalizeData(data.list[[i]])
  data.list[[i]] <- FindVariableFeatures(data.list[[i]])
}

total_NK <- RunFastMNN(object.list = data.list)
total_NK <- RunUMAP(total_NK, reduction = "mnn", dims = 1:30)
total_NK <- FindNeighbors(total_NK, reduction = "mnn", dims = 1:30) 
total_NK <- FindClusters(total_NK,resolution = 0.5)

DimPlot(total_NK, reduction = "umap", label=T,group.by = "RNA_snn_res.0.5")
DimPlot(total_NK,label = T,group.by = "cluster_name_slim")

total_NK$cluster_name_slim <- "NK_CD56bright"
total_NK$cluster_name_slim[total_NK$RNA_snn_res.0.5 %in% c(0)] <- "NK_CD56dim"
total_NK$cluster_name_slim[total_NK$RNA_snn_res.0.5 %in% c(4)] <- "NK_tr"

total_NK$cluster_name_slim <- factor(total_NK$cluster_name_slim,
                                     levels = c("NK_CD56bright","NK_CD56dim","NK_tr"))

total_NK$cohort <- factor(total_NK$cohort,levels = c("InHouse","AVTB2024"))

Idents(total_NK) <- "cluster_name_slim"
pdf(file = "Dimplot_NK_mergePublic.pdf",width = 10,height = 5)
DimPlot(total_NK,reduction = "umap",label = F,pt.size = 0.15,split.by = "cohort",cols = pal[1:3])+theme(legend.position = "bottom")
dev.off()

#cell proportion

cell.prop <- table(total_NK$cluster_name_slim,total_NK$cohort) %>% as.data.frame.matrix()
#cell.prop <- t(t(cell.prop)/colSums(cell.prop)) %>% as.data.frame()
cell.prop$cell <- row.names(cell.prop)
cell.prop <- pivot_longer(cell.prop,cols = setdiff(colnames(cell.prop),"cell"),names_to = "cohort",values_to = "num")
cell.prop = ddply(cell.prop, 'cohort', transform, percent=num/sum(num)*100)

cell.prop$cell <- factor(cell.prop$cell,levels = rev(levels(total_NK$cluster_name_slim)))
cell.prop$cohort <- factor(cell.prop$cohort,levels = levels(total_NK$cohort))

p <- ggplot(cell.prop,aes(x=cohort,y=percent,fill=cell))+
  geom_bar(stat = 'identity',width = 0.8,colour='black')+
  scale_fill_manual(values = rev(pal[1:3]))+
  scale_y_continuous(expand = c(0,0))+
  theme_classic()+
  theme(axis.text = element_text(size = 14,colour = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 16,colour = "black"))+
  ylab("Cell ratio (%)")+xlab("");p
ggsave("Cell_proportion_bar_NK_mergePublic.pdf",p,width = 4,height = 3.5)

qsave(total_NK,file = "total_NK_mergePublic.qs")

