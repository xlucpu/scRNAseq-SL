##### Figure plotting #####
## Figure 7A ## 
library(Seurat)
library(ggplot2)
library(dplyr)
setwd('D:/2021scRNA-seq')
load("seurat_Stromal.Rdat")
table(seurat$group)
###统计各cluster比例
prop.group<-as.data.frame(table(seurat$group))
table(seurat$patient)
prop.cluster<-as.data.frame(table(seurat@meta.data$seurat_clusters,seurat$group))
colnames(prop.cluster)<-c('cluster','group','freq.cluster')
prop_data <- merge(prop.group,prop.cluster,by.x=c("Var1"),by.y=c("group"))
prop_data$clusterprop<-100*prop_data$freq.cluster/prop_data$Freq
arrange(prop_data,cluster)
Stromal_metadata<-seurat@meta.data


seurat  <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 4000)
top10<-head(VariableFeatures(seurat),10)
top10
scale.genes <-  rownames(seurat)
seurat <- ScaleData(seurat, features = scale.genes)
seurat <- RunPCA(seurat , features = VariableFeatures(object = seurat ))
plot2 <- ElbowPlot(seurat, ndims=20, reduction="pca") 
plot2
pc.num=1:9
seurat <- FindNeighbors(seurat, dims = pc.num)
seurat <- FindClusters(seurat, resolution = 0.5)
table(seurat@meta.data$group)
table(seurat@meta.data$seurat_clusters,seurat@meta.data$group)

metadata <- seurat@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
write.csv(cell_cluster,'Stromal cluster/cell_cluster nfeatures4000 20210407.csv',row.names = F)

seurat <- RunUMAP(seurat, n.neighbors = 50,n.epochs=50,dims = pc.num)
embed_umap <- Embeddings(seurat, 'umap')
write.csv(embed_umap,'B cluster/embed_umap n.neighbors = 50,n.epochs=50 20210407.csv') 

#对细胞群添加注释信息
seurat@meta.data$celltypes<-seurat@meta.data$seurat_clusters
Idents(seurat) <- "celltypes" 
#DimPlot改颜色
load("Stromal cluster\\stromal.col.Rdata")
cell.col
seurat<- RenameIdents(seurat, '0'='PDGFRA+ fibroblasts','2'='Stromal-2','12'='EC-2','6'='Stromal-3',
                      '7'='Enteric glial cells','10'='Stalk-like ECs','13'='Myofibroblasts','3'='Stromal-1',
                      '4'='Tip-like ECs','5'='Stromal-1','11'='Vascular SMCs','8'='Stromal-1',
                      '9'='Pericytes','1'='Stromal-2','14'='EC-1',
                      '15'='EPCs')
seurat$celltypes<-Idents(seurat)
plot2<-DimPlot(seurat, reduction = "umap",  cols=c("#A60000","#ABDDA4","#FFFFBF","#531f7a","#5E4FA2",
                                                   "#FEE08B","#9E0142","#FDAE61","#66C2A5","#F46D43","#E6F598","#D53E4F","#3288BD"),
               group.by ='celltypes',label = T)
ggsave("Stromal cluster/UMAP  20211115.png", plot = plot2, width = 8/1.2, height = 5.5/1.2)



## Figure 7B ## 
load("seurat_Stromal.Rdat")
marker <- list("PDGRFA+ fibroblasts" = c("COL3A1","PDGFRA","POSTN","BMP2","BMP5","WNT5A","WNT5B"),
               "Stromal-1" = c("SFRP2"),
               "Stromal-2" = c("ADAMDEC1","APOE","MDK"),
               "Stromal-3" = c("CCL5","CCL4","IL7R"),
               "Myofibroblasts" = c("TAGLN","ACTA2","MYL9"),
               "Stalk-like ECs" = c("ACKR1","SELP"),
               "Tip-like ECs" = c("RAMP3","RGCC"),
               "EC-1" = c("PECAM1","ENG","ICAM1"),
               "EC-2" = c("CD36"),
               "EPCs" = c("HMGB2","CKS2","TYROBP","C1QB","RAC2"),
               "Vascular SMCs" = c("CNN1","SYNPO2","DES"),
               "Pericytes" = c("RGS5","PDGFRB","CSPG4","ABCC9"),
               "Enteric glial cells" = c("SOX10","S100B"))

module.gene <- sapply(marker,length)
cellLevel <- names(marker)
marker <- unlist(marker[cellLevel])
annotation <- data.frame(cell = seurat$orig.ident,
                         cluster = as.numeric(as.character(seurat$seurat_clusters)),
                         celltype = seurat$celltypes)

# extract averaged expression
marker.expr.celltype <- AverageExpression(seurat,assays = "RNA",features = marker,group.by = "celltypes",verbose = TRUE)
marker.expr.celltype <- as.data.frame(marker.expr.celltype$RNA)
marker.expr.celltype <- marker.expr.celltype[,cellLevel]
marker.expr.celltype <- marker.expr.celltype[marker,]

#load("heatmap_AveE_Stromal.Rdat")
#marker.expr.celltype <- heatmap_AveE[marker,]
annCol <- data.frame(celltype = factor(names(module.gene),levels = names(module.gene)),
                     row.names = names(module.gene),
                     stringsAsFactors = F)
annColors <- list()
annColors[["celltype"]] <- cell.col

plotdata <- standarize.fun(marker.expr.celltype,halfwidth = 1)
hm <- pheatmap(as.matrix(plotdata),
               border_color = NA,
               color = NMF:::ccRamp(heatmap.BlBkRd,64),
               #color = viridis(64),
               cluster_rows = F,
               cluster_cols = F,
               show_rownames = T,
               show_colnames = T,
               gaps_col = c(1,4,5,9,10,11,12),
               gaps_row = c(7,14,17,25,30,33,37,39),
               cellwidth = 13,
               cellheight = 10,
               name = "AvgExpr.\n[z-scored]",
               annotation_col = annCol,
               annotation_colors = annColors["celltype"])
pdf("heatmap of marker in stromal cells2.pdf",width = 6,height = 8)
draw(hm)
invisible(dev.off())



## Figure 7D ## 
path <- "~/plot/VolcanicPlot";setwd(path)
data <- openxlsx::read.xlsx("PDGFRA+ fibroblasts20210828.xlsx")

library(ggplot2)
# jco <- c("#2874C5","#EABF00","#868686","#C6524A","#80A7DE")

plotdata <- data.frame(
  "gene" = data$X1,
  "log2FC" = data$avg_log2FC,
  "logPadj" = (-log10(data$p_val_adj + 2.225074e-308))
  # "loglog" = log10(-log10(data$p_val_adj + 2.225074e-308)+1)
)
plotgene <- c("PDGFRA", "SOX6", "POSTN", "PDGFD", "IGFBP3", "BMP5", "WNT5B",
              "MMP11", "BMP4", "TSLP", "BMP2", "CCL2", "SEMA4D", "APOE", "ITGA1", "GPNMB", "WNT5A", "MMP1")

plotdata$label[plotdata$gene %in% plotgene] <- plotdata$gene[plotdata$gene %in% plotgene]

logFCcut <- c(-2, -0.5, 0.5, 2)
pvalcut <- c(50, 100)
plotdata$color <- "#868686"
plotdata$color[plotdata$log2FC < logFCcut[2] & plotdata$logPadj > pvalcut[1]] <- "#4966AE"
plotdata$color[plotdata$log2FC > logFCcut[3] & plotdata$logPadj > pvalcut[1]] <- "#E9281E"

plotdata$size <- abs(plotdata$log2FC * plotdata$logPadj)
plotdata$size <- plotdata$size*7/max(plotdata$size)+1


#################################
# c("#3B00FF", "#9376F2", "#A9A9A9", "#F5897E", "#F61600")
col <- names(table(plotdata$color))
names(col) <- col

ggplot(data = plotdata, aes(x = log2FC, y = logPadj, color = color, size = size)) + 
  geom_point(alpha = 0.6) +
  scale_color_manual(values = col) +
  scale_size_continuous(range=c(1, 8)) +
  # ggrepel::geom_text_repel(aes(x = log2FC, y = logPadj, label = label,
  #                              box.padding = unit(0.35, "lines"),  point.padding = unit(0.3, "lines")))+
  ggrepel::geom_label_repel(aes(label = label), box.padding = unit(0.35, "lines"), 
                            point.padding = unit(0.3, "lines"), size = 4)+
  geom_vline(xintercept = logFCcut[c(2,3)], color="grey40", linetype="longdash", lwd = 0.25) +
  geom_hline(yintercept = pvalcut[1], color="grey40", linetype="longdash", lwd = 0.25) +
  ggtitle("PDGFRA+ fibroblasts") +
  theme_bw() +
  theme(panel.grid=element_blank(), plot.title = element_text(hjust = 0.5), legend.position = "none") +
  geom_point(data = subset(plotdata, !is.na(label)),  alpha = 1, shape = 1, 
             stroke = 0.5, #圈粗细
             color = "black")
ggsave("PDGFRA+ fibroblasts0922.pdf", width = 793/100, height = 556/100)



## Figure 7E ## 
library(CellChat)
library(patchwork)
library(openxlsx)
expr <- readRDS("i:/genomicdata/jiaoda/raw/expr0617.rds")
cell_info <- readRDS("i:/genomicdata/jiaoda/raw/cell_info.rds")
load("i:/genomicdata/jiaoda/raw/CD8T_metadata.Rdata")
expr <- expr[,rownames(cell_info)]
all(rownames(cell_info)%in%row)
# seu <- CreateSeuratObject(expr)
# seu <- NormalizeData(seu)
# Norm.expr <- GetAssayData(seu, slot = "data")
B <- table(subset(cell_info,type == "B")["celltypes"])
B <- names(B[B>0])
epi <- table(subset(cell_info,type == "epithelic")["celltypes"])
epi <- names(epi[epi>0])
mye <- table(subset(cell_info,type == "myeloid")["celltypes"])
mye <- names(mye[mye>0])
stro <- table(subset(cell_info,type == "stromal")["celltypes"])
stro <- names(stro[stro>0])
CD8T <- table(subset(cell_info,division == "CD8+ T cells")["celltypes"])
CD8T <- names(CD8T[CD8T>0])
panT <- table(subset(cell_info,type == "panT")["celltypes"])
panT <- names(panT[panT>0])

# dotplot for one pair in multiple celltypes
ArrayToDataFrame <- function(array,pair = NULL){
  source = dimnames(array)[[1]]
  target = dimnames(array)[[2]]
  if (is.null(pair)) pair = dimnames(array)[[3]]
  data = lapply(source,function(i){
    lapply(target,function(j){
      t(as.matrix(sapply(pair,function(k){
        c(i, j, k, as.numeric(array[i,j,k]))#source, target, pair, value
      })))
    })
  })
  data = lapply(data,function(x) do.call(rbind,x))
  data = do.call(rbind,data)
  data = as.data.frame(data)
  rownames(data) = NULL
  colnames(data) = c("source", "target", "pair", "value")
  return(data)
}
ExtractCellChat <- function(cellchat,pair){
  prob <- cellchat@net$prob
  prob <- ArrayToDataFrame(prob,pair)
  pval <- cellchat@net$pval
  pval <- ArrayToDataFrame(pval,pair)
  if (all(prob[,1:3]==pval[,1:3])){
    res = data.frame(
      "source" = prob$source,
      "target" = prob$target,
      "pair" = prob$pair,
      "prob" = prob$value,
      "pval" = pval$value
    )
  }else{print("index not match")}
}
FetchInteractionScoreForDataFrame <- function(data){
  data = split(data,f = paste0(data$source, data$target, data$pair))
  data = lapply(data, function(x){
    source = as.character(x[1]); target = as.character(x[2])
    ligand = strsplit(as.character(x[3]), "_")[[1]][1]
    receptor = strsplit(as.character(x[3]), "_")[[1]][2]
    source.expr = expr[ligand,cell_info$celltypes%in%source]
    target.expr = expr[receptor,cell_info$celltypes%in%target]
    score = (mean(source.expr)+mean(target.expr))/2
    x["score"] = score
    return(x)
  })
  data = as.data.frame(do.call(rbind,data))
  rownames(data) = NULL
  colnames(data) = c("source", "target", "pair", "prob", "pval", "score")
  return(data)
}

path <- "i:/genomicdata/jiaoda/cellchat/old";setwd(path)
project <- c("stro_panT_all","stro_B_all","stro_mye_all")
cellchatdata <- lapply(project, function(x){
  cellchat = readRDS(file.path(path, paste0(x,".rds")))
  data = ExtractCellChat(cellchat,"CXCL12_CXCR4")
  # data = ExtractCellChat(cellchat,"CCL5_CCR1")
})
cellchatdata <- do.call(rbind,cellchatdata)
cellchatdata = subset(cellchatdata, source%in%stro & target%in%c(panT,B,mye))
cellchatdata <- FetchInteractionScoreForDataFrame(cellchatdata)

plotdata <- cellchatdata
plotdata$source <- factor(plotdata$source,
                          levels = c("PDGRFA+ fibroblasts", "Stromal-1", "Stromal-2",
                                     "Stromal-3", "Myofibroblasts", "Stalk-like ECs", 
                                     "Tip-like ECs", "EC-1", "EC-2", "EPCs", "Vascular SMCs",
                                     "Pericytes", "Enteric glial cells"))
plotdata$target <- factor(plotdata$target, levels = c(panT, B, "FCN1+ inflammatory macrophages","Anti-inflammatory macrophages","Conventional DCs","Plasmacytoid DCs"))
plotdata$P <- cut(as.numeric(plotdata$pval),c(1,0.05,0.01,0),include.lowest = T)
plotdata$P <- factor(plotdata$P,levels = c("(0.05,1]","(0.01,0.05]","[0,0.01]"))
plotdata$score[plotdata$score>5]=5

CXCL12=plotdata
CXCL12 = subset(CXCL12, !target=="CD8+ Naïve")
CXCL12$target <- factor(CXCL12$target, 
                        levels = c(panT, B, c("FCN1+ inflammatory macrophages","Anti-inflammatory macrophages","Conventional DCs","Plasmacytoid DCs"))) 
grDevices::cairo_pdf(file.path(path,"CXCL12-CXCR4_1211.pdf"), width = 769/100, height = 472/100)
ggplot(CXCL12, aes(x = target, y = source, col = score, size = P))+
  geom_point()+
  scale_colour_gradientn(colors = colorRampPalette(c("darkblue","yellow","red"))(99), na.value = "white") +
  scale_y_discrete(limits = rev)+
  theme_linedraw() + theme(panel.grid.major = element_blank())+
  theme(axis.text.y = element_text(color = "black", size = 10), axis.title = element_blank(),
        axis.text.x = element_text(color = "black", size = 10, angle = 45, hjust = 1, vjust = 1))+
  geom_vline(xintercept = c(10.5,14.5))+
  ggtitle("CXCL12-CXCR4")
invisible(dev.off())


