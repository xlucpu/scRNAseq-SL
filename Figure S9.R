####Figure S9C
path <- "i:/genomicdata/jiaoda/SCENIC/output/";setwd(path)
raw.path <- "i:/genomicdata/jiaoda/raw"
load(file.path(raw.path, "seurat_Myeloid.Rdata"))
# load(file.path(raw.path, "seurat_B.Rdata"))
# count <- readRDS("~/genomicdata/jiaoda/expr0617.rds")
# auc.seu <- readRDS("~/genomicdata/jiaoda/SCENIC/auc.seu.rds")
aucmtx <- read.csv("auc_mtx.csv")
bimtx <- read.table("~/genomicdata/jiaoda/SCENIC/output/binary_info.txt",sep = "\t", header = T, row.names = 1)
colnames(aucmtx) <- substr(colnames(aucmtx),1,nchar(colnames(aucmtx))-3)
rownames(aucmtx) <- aucmtx[,1]
rownames(aucmtx) <- gsub("\\.","-",rownames(aucmtx))
info <- as.data.frame(aucmtx[,c("ATF5","CREB3L2")])
info <- info[colnames(seurat),]
info[info>=0.17]
colnames(info) <- paste0("AUCell_",colnames(info))
seurat <- AddMetaData(seurat, metadata = info)
seurat$normalized_ATF5 <- seurat@assays$RNA@counts["ATF5",]
seurat$normalized_ATF5[seurat$normalized_ATF5>1] <- 1

pdf(paste0(path,"ATF5.pdf"),width = 5, height = 5)
FeaturePlot(seurat,features = "AUCell_ATF5",cols = c("grey","#DA2D31"), min.cutoff = 0, max.cutoff = 0.2)+
  scale_color_gradient2(low="grey",mid="grey",high="#DA2D31",midpoint = 0.1, limits = c(0,0.2))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+NoLegend()
FeaturePlot(seurat, features = "normalized_ATF5",cols = c("grey","#DA2D31"))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+NoLegend()
invisible(dev.off())