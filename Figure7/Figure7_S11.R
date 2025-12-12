######################Figure 7 CellChat##############################
library(CellChat)

CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")

#total cells
total <- readRDS("total_final.2.rds")

data.input <- as.matrix(total@assays$SCT@data)
meta <- total@meta.data

for(x in c("Ctl","CAS")){
  cell.use <- rownames(meta)[meta$disease == x]
  data.input_2 <- data.input[, cell.use]
  meta_2 = meta[cell.use, ]
  
  cellchat <- createCellChat(object = data.input_2, meta = meta_2, group.by = "cluster_name_slim")
  groupSize <- as.numeric(table(cellchat@idents))
  cellchat@DB <- CellChatDB.use
  
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)
  cellchat <- computeCommunProb(cellchat, raw.use = T,type = "truncatedMean",trim = 0.2,population.size = F)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  df.net <- subsetCommunication(cellchat)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  assign(paste0("cellchat_",x),cellchat)
}


object.list <- list(Ctl = cellchat_Ctl, CAS = cellchat_CAS)
saveRDS(object.list, file="cellchat.object.list.rds")

cellchat <- mergeCellChat(object.list, add.names = names(object.list))

weight.max <- getMaxWeight(object.list, attribute = c("idents","weight"))
pdf("FigureS11A_InteractionWeights_ADSC.pdf",width = 10,height = 5)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$weight,idents.use = c("ADSC"), weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Interaction weights - ", names(object.list)[i]))
}
dev.off()

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(unlist(num.link)), max(unlist(num.link))) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]], slot.name = "netP")
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)+xlim(0,3.5)+ylim(0,2)
}
pdf("FigureS11B_OutIncomeStrength.pdf",width = 10,height = 5)
patchwork::wrap_plots(plots = gg)
dev.off()

pdf("FigureS11C_InteractionStrengthDiff.pdf",width = 5,height = 4)
netVisual_heatmap(cellchat, measure = "weight")
dev.off()

#ADSC and innate immune cell
total_immune <- readRDS("total_innate_immune.rds")
total_ADSC <- readRDS("total_ADSC.rds")
total <- merge(total_ADSC,total_immune)

data.input <- as.matrix(total@assays$SCT@data)
meta <- total@meta.data

for(x in c("Ctl","CAS")){
  cell.use <- rownames(meta)[meta$disease == x]
  data.input_2 <- data.input[, cell.use]
  meta_2 = meta[cell.use, ]
  
  cellchat <- createCellChat(object = data.input_2, meta = meta_2, group.by = "cluster_name_slim")
  groupSize <- as.numeric(table(cellchat@idents))
  cellchat@DB <- CellChatDB.use
  
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)
  cellchat <- computeCommunProb(cellchat, raw.use = T,type = "truncatedMean",trim = 0.2,population.size = F)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  df.net <- subsetCommunication(cellchat)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  assign(paste0("cellchat_",x),cellchat)
}


object.list <- list(Ctl = cellchat_Ctl, CAS = cellchat_CAS)
saveRDS(object.list, file="cellchat.object.list.innate.rds")

object.list <- readRDS("cellchat.object.list.innate.rds")
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

weight.max <- getMaxWeight(object.list, attribute = c("idents","weight"))
pdf("Figure7D_InteractionWeights_ADSC_innate.pdf",width = 10,height = 5)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$weight,idents.use = c("ADSC_CXCL14","ADSC_CD55"), weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Interaction weights - ", names(object.list)[i]))
}
dev.off()

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(unlist(num.link)), max(unlist(num.link))) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]], slot.name = "netP")
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)+xlim(0,6)+ylim(0,3)
}
pdf("Figure7E_OutIncomeStrength_innate.pdf",width = 10,height = 5)
patchwork::wrap_plots(plots = gg)
dev.off()

pdf("Figure7F_InteractionStrengthDiff_innate.pdf",width = 5,height = 4)
netVisual_heatmap(cellchat, measure = "weight")
dev.off()

p1 <- netVisual_bubble(cellchat,sources.use = c(1,2),targets.use = c(3:12),comparison = c(1, 2),max.dataset = 2, title.name = "Increased signaling in CAS",remove.isolate = T,color.text = c("#00BFC4","#F8766D"))
p2 <- netVisual_bubble(cellchat,sources.use = c(1,2),targets.use = c(3:12),comparison = c(1, 2),max.dataset = 1, title.name = "Decreased signaling in CAS",remove.isolate = T,color.text = c("#00BFC4","#F8766D"))
pdf("Figure7H_Bubble_ADSC_to_Myeloid.pdf",width = 15,height = 7)
p1+p2
dev.off()

p1 <- netVisual_bubble(cellchat,sources.use = c(3:12),targets.use = c(1,2),comparison = c(1, 2),max.dataset = 2, title.name = "Increased signaling in CAS",remove.isolate = T,color.text = c("#00BFC4","#F8766D"))
p2 <- netVisual_bubble(cellchat,sources.use = c(3:12),targets.use = c(1,2),comparison = c(1, 2),max.dataset = 1, title.name = "Decreased signaling in CAS",remove.isolate = T,color.text = c("#00BFC4","#F8766D"))
pdf("FigureS11E_Bubble_Myeloid_to_ADSC.pdf",width = 15,height = 7)
p1+p2
dev.off()


df.net <- subsetCommunication(object.list[[2]])

pathways.show <- c("PTN") 
LR.show <- c("PTN_NCL")
pdf("Pathway_ChordPlot_innate_PTN.pdf",width = 5,height = 5)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_individual(object.list[[i]], signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
}
dev.off()

pathways.show <- c("ANNEXIN") 
LR.show <- c("ANXA1_FPR1")
pdf("Pathway_ChordPlot_innate_ANXA1.pdf",width = 5,height = 5)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_individual(object.list[[i]], signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
}
dev.off()

pathways.show <- c("COMPLEMENT") 
LR.show <- c("C3_ITGAX_ITGB2")
pdf("Pathway_ChordPlot_innate_C3.pdf",width = 5,height = 5)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_individual(object.list[[i]], signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
}
dev.off()

pathways.show <- c("MIF") 
LR.show <- c("MIF_CD74_CXCR4")
pdf("Pathway_ChordPlot_innate_MIF.pdf",width = 5,height = 5)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_individual(object.list[[i]], signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
}
dev.off()

pathways.show <- c("VEGF") 
LR.show <- c("VEGFA_VEGFR1")
pdf("Pathway_ChordPlot_innate_VEGFA.pdf",width = 5,height = 5)
netVisual_individual(object.list[[2]], signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
dev.off()

pathways.show <- c("VEGF") 
LR.show <- c("VEGFB_VEGFR1")
pdf("Pathway_ChordPlot_innate_VEGFB.pdf",width = 5,height = 5)
netVisual_individual(object.list[[2]], signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
dev.off()

list <- c("LGALS9_CD44","EREG_EGFR","AREG_EGFR","MIF_ACKR3","NAMPT_ITGA5_ITGB1","TNFSF12_TNFRSF12A",
          "GAS6_AXL","CXCL12_ACKR3","CXCL12_CXCR4","MIF_CD74_CXCR4","C3_ITGAX_ITGB2","ANXA1_FPR1",
          "PTN_NCL","ANGPTL1_PIRB","CSF1_CSF1R","FGF2_FGFR1","RARRES2_CMKLR1","SEMA3C_PLXND1")

df.net.sub <- df.net[df.net$interaction_name %in% list,c("interaction_name","pathway_name")] %>% distinct()

for(x in 1:nrow(df.net.sub)){
  pathways.show <- df.net.sub$pathway_name[x] 
  LR.show <- df.net.sub$interaction_name[x]
  pdf(paste0("Figure7_Pathway_ChordPlot_innate_",LR.show,".pdf"),width = 5,height = 5)
  par(mfrow = c(1,2), xpd=TRUE)
  for (i in 1:length(object.list)) {
    netVisual_individual(object.list[[i]], signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
  }
  dev.off()
}



#ADSC and T cell
total_T <- readRDS("total_T.rds")
total_ADSC <- readRDS("total_ADSC.rds")

total <- merge(total_ADSC,total_T)

data.input <- as.matrix(total@assays$SCT@data)
meta <- total@meta.data

for(x in c("Ctl","CAS")){
  cell.use <- rownames(meta)[meta$disease == x]
  data.input_2 <- data.input[, cell.use]
  meta_2 = meta[cell.use, ]
  
  cellchat <- createCellChat(object = data.input_2, meta = meta_2, group.by = "cluster_name_slim")
  groupSize <- as.numeric(table(cellchat@idents))
  cellchat@DB <- CellChatDB.use
  
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)
  cellchat <- computeCommunProb(cellchat, raw.use = T,type = "truncatedMean",trim = 0.2,population.size = F)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  df.net <- subsetCommunication(cellchat)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  assign(paste0("cellchat_",x),cellchat)
}


object.list <- list(Ctl = cellchat_Ctl, CAS = cellchat_CAS)
saveRDS(object.list, file="cellchat.object.list.T.rds")

object.list <- readRDS("cellchat.object.list.T.rds")
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

weight.max <- getMaxWeight(object.list, attribute = c("idents","weight"))
pdf("Figure7A_InteractionWeights_ADSC_T.pdf",width = 10,height = 5)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$weight,idents.use = c("ADSC_CXCL14","ADSC_CD55"), weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Interaction weights - ", names(object.list)[i]))
}
dev.off()

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(unlist(num.link)), max(unlist(num.link))) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]], slot.name = "netP")
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)+xlim(0,3)+ylim(0,2)
}
pdf("Figure7B_OutIncomeStrength_T.pdf",width = 10,height = 5)
patchwork::wrap_plots(plots = gg)
dev.off()

pdf("Figure7C_InteractionStrengthDiff_T.pdf",width = 5,height = 4)
netVisual_heatmap(cellchat, measure = "weight")
dev.off()

p1 <- netVisual_bubble(cellchat,sources.use = c(1,2),targets.use = c(3:12),comparison = c(1, 2),max.dataset = 2, title.name = "Increased signaling in CAS",remove.isolate = T,color.text = c("#00BFC4","#F8766D"))
p2 <- netVisual_bubble(cellchat,sources.use = c(1,2),targets.use = c(3:12),comparison = c(1, 2),max.dataset = 1, title.name = "Decreased signaling in CAS",remove.isolate = T,color.text = c("#00BFC4","#F8766D"))
pdf("Figure7G_Bubble_ADSC_to_T.pdf",width = 15,height = 3.5)
p1+p2
dev.off()

p1 <- netVisual_bubble(cellchat,sources.use = c(3:12),targets.use = c(1,2),comparison = c(1, 2),max.dataset = 2, title.name = "Increased signaling in CAS",remove.isolate = T,color.text = c("#00BFC4","#F8766D"))
p2 <- netVisual_bubble(cellchat,sources.use = c(3:12),targets.use = c(1,2),comparison = c(1, 2),max.dataset = 1, title.name = "Decreased signaling in CAS",remove.isolate = T,color.text = c("#00BFC4","#F8766D"))
pdf("FigureS11D_Bubble_T_to_ADSC.pdf",width = 15,height = 6)
p1+p2
dev.off()

df.net <- subsetCommunication(object.list[[2]])

pathways.show <- c("CXCL") 
LR.show <- c("CXCL12_CXCR4")
pdf("Figure7I_Pathway_CirclePlot_T_CXCL.pdf",width = 5,height = 5)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_individual(object.list[[i]], signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
}
dev.off()

pathways.show <- c("PTN") 
LR.show <- c("PTN_NCL")
pdf("Pathway_CirclePlot_T_PTN.pdf",width = 5,height = 5)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_individual(object.list[[i]], signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
}
dev.off()

pathways.show <- c("MIF") 
LR.show <- c("MIF_CD74_CXCR4")
pdf("Pathway_CirclePlot_T_MIF.pdf",width = 5,height = 5)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_individual(object.list[[i]], signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
}
dev.off()

list <- c("MIF_ACKR3","MIF_CD74_CD44","IFNG_IFNGR1_IFNGR2","CD40LG_ITGA5_ITGB1",
          "AREG_EGFR","CD70_CD27","PTN_NCL","BTLA_TNFRSF14")

df.net.sub <- df.net[df.net$interaction_name %in% list,c("interaction_name","pathway_name")] %>% distinct()

for(x in c(2,3,4,6,7,8)){
  pathways.show <- df.net.sub$pathway_name[x] 
  LR.show <- df.net.sub$interaction_name[x]
  pdf(paste0("Figure7_Pathway_ChordPlot_T_",LR.show,".pdf"),width = 5,height = 5)
  par(mfrow = c(1,2), xpd=TRUE)
  for (i in 1:length(object.list)) {
    netVisual_individual(object.list[[i]], signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
  }
  dev.off()
}

x=1
pathways.show <- df.net.sub$pathway_name[x] 
LR.show <- df.net.sub$interaction_name[x]
pdf(paste0("Figure6_Pathway_ChordPlot_T_",LR.show,"_onlyCAS.pdf"),width = 5,height = 5)
i=2
netVisual_individual(object.list[[i]], signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
dev.off()

x=5
pathways.show <- df.net.sub$pathway_name[x] 
LR.show <- df.net.sub$interaction_name[x]
pdf(paste0("Figure6_Pathway_ChordPlot_T_",LR.show,"_onlyCAS.pdf"),width = 5,height = 5)
i=2
netVisual_individual(object.list[[i]], signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
dev.off()