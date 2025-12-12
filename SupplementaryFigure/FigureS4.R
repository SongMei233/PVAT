######################Merge pla data######################

x <- "Patient 1"
count1 <- Read10X(data.dir=paste0(x,"/C230517003-pla"),gene.column=1) 
count1 <- CreateSeuratObject(counts = as.sparse(count1),project=x)
count2 <- Read10X(data.dir=paste0(x,"/C230517004-pla"),gene.column=1)
count2 <- CreateSeuratObject(counts = as.sparse(count2),project=x)
count3 <- Read10X(data.dir=paste0(x,"/C230517011-pla"),gene.column=1)
count3 <- CreateSeuratObject(counts = as.sparse(count3),project=x)
count4 <- Read10X(data.dir=paste0(x,"/C230517012-pla"),gene.column=1)
count4 <- CreateSeuratObject(counts = as.sparse(count4),project=x)

count <- merge(count1,list(count2,count3,count4))
count <- count@assays$RNA$counts

index <- paste(strsplit(colnames(count),"_") %>% lapply(.,function(x){x[1]}) %>% unlist(),
               strsplit(colnames(count),"_") %>% lapply(.,function(x){x[2]}) %>% unlist(),
               sep = "_")
count <- as.data.frame(count) %>% t() %>% as.data.frame()
count <- as.sparse(count)
count <- aggregate(count,by = list(index),FUN = sum)
count <- column_to_rownames(count,var = "Group.1")
count <- t(count)
data <- CreateSeuratObject(counts = as.sparse(count),project=x,min.cells=10)
data <- RenameCells(object=data,add.cell.id="Patient1")
data$sample <- "Patient1"
data <- PercentageFeatureSet(data, pattern = "^MT-", col.name = "percent.mt")
data <- subset(data, subset = nFeature_RNA > 400 & nFeature_RNA < 8000 & percent.mt < 15)
data <- SCTransform(data)

saveRDS(data,file = "Patient1.rds")

x <- "Patient 2"
count1 <- Read10X(data.dir=paste0(x,"/C231227001-pla"),gene.column=1) 
count1 <- CreateSeuratObject(counts = as.sparse(count1),project=x)
count2 <- Read10X(data.dir=paste0(x,"/C231227002-pla"),gene.column=1)
count2 <- CreateSeuratObject(counts = as.sparse(count2),project=x)
count3 <- Read10X(data.dir=paste0(x,"/C231227007-pla"),gene.column=1)
count3 <- CreateSeuratObject(counts = as.sparse(count3),project=x)
count4 <- Read10X(data.dir=paste0(x,"/C231227008-pla"),gene.column=1)
count4 <- CreateSeuratObject(counts = as.sparse(count4),project=x)

count <- merge(count1,list(count2,count3,count4))
count <- count@assays$RNA$counts

index <- paste(strsplit(colnames(count),"_") %>% lapply(.,function(x){x[1]}) %>% unlist(),
               strsplit(colnames(count),"_") %>% lapply(.,function(x){x[2]}) %>% unlist(),
               sep = "_")
count <- as.data.frame(count) %>% t() %>% as.data.frame()
count <- as.sparse(count)
count <- aggregate(count,by = list(index),FUN = sum)
count <- column_to_rownames(count,var = "Group.1")
count <- t(count)
data <- CreateSeuratObject(counts = as.sparse(count),project=x,min.cells=10)
data <- RenameCells(object=data,add.cell.id="Patient2")
data$sample <- "Patient2"
data <- PercentageFeatureSet(data, pattern = "^MT-", col.name = "percent.mt")
data <- subset(data, subset = nFeature_RNA > 400 & nFeature_RNA < 8000 & percent.mt < 15)
data <- SCTransform(data)

saveRDS(data,file = "Patient2.rds")

x <- "Patient 3"
count1 <- Read10X(data.dir=paste0(x,"/C231122001-pla"),gene.column=1) 
count1 <- CreateSeuratObject(counts = as.sparse(count1),project=x)
count2 <- Read10X(data.dir=paste0(x,"/C231122002-pla"),gene.column=1)
count2 <- CreateSeuratObject(counts = as.sparse(count2),project=x)


count <- merge(count1,count2)
count <- count@assays$RNA$counts

index <- paste(strsplit(colnames(count),"_") %>% lapply(.,function(x){x[1]}) %>% unlist(),
               strsplit(colnames(count),"_") %>% lapply(.,function(x){x[2]}) %>% unlist(),
               sep = "_")
count <- as.data.frame(count) %>% t() %>% as.data.frame()
count <- as.sparse(count)
count <- aggregate(count,by = list(index),FUN = sum)
count <- column_to_rownames(count,var = "Group.1")
count <- t(count)
data <- CreateSeuratObject(counts = as.sparse(count),project=x,min.cells=10)
data <- RenameCells(object=data,add.cell.id="Patient3")
data$sample <- "Patient3"
data <- PercentageFeatureSet(data, pattern = "^MT-", col.name = "percent.mt")
data <- subset(data, subset = nFeature_RNA > 400 & nFeature_RNA < 8000 & percent.mt < 15)
data <- SCTransform(data)

saveRDS(data,file = "Patient3.rds")

data.list <- list()
data.list[[1]] <- readRDS("Patient1.rds")
data.list[[2]] <- readRDS("Patient2.rds")
data.list[[3]] <- readRDS("Patient3.rds")

total2 <- readRDS("total_final.2.rds")

options(future.globals.maxSize= 100000000000000)

data.features <- SelectIntegrationFeatures(object.list = data.list)
data.list <- PrepSCTIntegration(object.list = data.list,assay = "SCT", anchor.features = data.features,verbose = FALSE)
data.anchors <- FindIntegrationAnchors(object.list = data.list, normalization.method = "SCT", anchor.features = data.features, verbose = FALSE)
total <- IntegrateData(anchorset = data.anchors, normalization.method = "SCT", verbose = FALSE)

total <- Reduct(total)

DimPlot(total, reduction = "umap", label=T,group.by = "integrated_snn_res.1")

DefaultAssay(total) <- "SCT"
Idents(total) <- "integrated_snn_res.1"
markers <- c("PTPRC","MS4A1","CD19","BANK1","JCHAIN","SDC1","IGKC","CD3D","CD4","CCR7","SELL","CD8A","GZMK","CCL5","KLRD1","GNLY","NKG7","MKI67","STMN1","TUBA1B","IL1B","S100A8","S100A9","CD68","C1QA","LYVE1","CD1E","LAMP3","FSCN1","PTGDS","IL3RA","LILRA4","F3","ITGB1","CD34","COL1A1","CXCL14","CD55","HBB","HBA1","VWF","PECAM1")
DotPlot(total, features = markers, dot.scale = 6)  + theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) + RotatedAxis()

total <- PrepSCTFindMarkers(total)
clu6.markers <- FindMarkers(total,ident.1 = "Pericyte",min.pct = 0.5, min.diff.pct = 0.2, logfc.threshold =0.25, only.pos = TRUE, test.use = "wilcox")

total$cluster_name_slim <- "CD4T cell"
total$cluster_name_slim[total$integrated_snn_res.0.5 %in% c(12)] <- "B cell"
total$cluster_name_slim[total$integrated_snn_res.0.5 %in% c(18)] <- "Plasma cell"
total$cluster_name_slim[total$integrated_snn_res.0.5 %in% c(2,3)] <- "CD8T cell"
total$cluster_name_slim[total$integrated_snn_res.1 %in% c(12)] <- "NK"
total$cluster_name_slim[total$integrated_snn_res.0.5 %in% c(15)] <- "Proliferating Lymphocyte"
total$cluster_name_slim[total$integrated_snn_res.0.5 %in% c(5)] <- "Monocyte/Neutrophil"
total$cluster_name_slim[total$integrated_snn_res.0.5 %in% c(4,7)] <- "Macrophage"
total$cluster_name_slim[total$integrated_snn_res.0.5 %in% c(9)] <- "DC"
total$cluster_name_slim[total$integrated_snn_res.0.5 %in% c(17)] <- "pDC"
total$cluster_name_slim[total$integrated_snn_res.0.5 %in% c(0,11)] <- "Fibroblast"
total$cluster_name_slim[total$integrated_snn_res.0.5 %in% c(6)] <- "Endothelial cell"
total$cluster_name_slim[total$integrated_snn_res.0.5 %in% c(10)] <- "Mast cell"
total$cluster_name_slim[total$integrated_snn_res.0.5 %in% c(8,13,14,16)] <- "Doublets"

total <- subset(total,subset = cluster_name_slim != "Doublets")

total$cluster_name_slim <- factor(total$cluster_name_slim,
                                  levels = c("B cell","Plasma cell","CD4T cell","CD8T cell","NK","Proliferating Lymphocyte","Monocyte/Neutrophil","Macrophage","DC","pDC","Mast cell","Fibroblast","Endothelial cell"))

DefaultAssay(total) <- "integrated"
total <- RunUMAP(total, dims = 1:30)

DimPlot(total,group.by = "cluster_name_slim",label=T)

qsave(total,file = "total_pla.qs")

Idents(total) <- "cluster_name_slim"

pal = DiscretePalette_scCustomize(num_colors = 14, palette = "ditto_seq")[-c(8)]

pdf(file = "Dimplot_mergePLA.pdf",width = 8,height = 6)
DimPlot_scCustom(total,colors_use = pal[1:13],figure_plot = T,label = T,label.size = 5,pt.size = 0.5)
dev.off()

markers <- c("PTPRC","MS4A1","CD19","BANK1","JCHAIN","SDC1","IGKC","CD3D","CD4","CCR7","SELL","CD8A","GZMK","CCL5","KLRD1","GNLY","NKG7","MKI67","STMN1","TUBA1B","IL1B","S100A8","S100A9","CD68","C1QA","LYVE1","CD1E","LAMP3","FSCN1","PTGDS","IL3RA","LILRA4","TPSB2","TPSAB1","CPA3","ITGB1","COL1A1","FN1","CD34","PECAM1","VWF")
Idents(total) <- "cluster_name_slim"
DefaultAssay(total) <- "SCT"
pdf(file = "DotPlot_mergePLA.pdf",width = 14,height = 5)
DotPlot(total, features = markers, dot.scale = 6)  + theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) + RotatedAxis()
dev.off()