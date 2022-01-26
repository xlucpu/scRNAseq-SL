##### Figure plotting #####
## Figure S10A ## 
library(Seurat)
library(ggplot2)
library(dplyr)
setwd('D:/2021scRNA-seq')
load("seurat_Stromal.Rdat")
plot2<-DimPlot(seurat, reduction = "umap",   group.by ='seurat_clusters',label = T)
ggsave("Stromal cluster/UMAP  20211211.png", plot = plot2, width = 6/1.2, height = 5.5/1.2)

Idents(seurat) <- "group" 
seurat<- RenameIdents(seurat, 'NC'='NC','HP'='HP','SSL'='SSL','SSLD'='SSLD',
                      'TSA'='TSA')
seurat$groups<-Idents(seurat)
plot2<-DimPlot(seurat, reduction = "umap",
               group.by ='groups',label = T)
plot2
ggsave("Stromal cluster/UMAP 20220119 group.pdf", plot = plot2, width = 6, height = 5)



## Figure S10B ## 
plot3<-FeaturePlot(seurat, features = c('MMP3','MMP11'))
ggsave("Stromal cluster/MMP  20211211.png", plot = plot2, width = 6/1.2, height = 5.5/1.2)


## Figure S10C ## 
####BRAF-mut CRC间质
setwd('D:/2021scRNA-seq+3 BRAF-MUT')
load('seurat_fullCRC.Rdat')
#通过直接建模在单细胞数据固有的均方差关系改进了以前的版本，并在实现FindVariableFeatures功能。V3默认选择2000个差异基因。这些将用于下游分析，例如PCA。
seurat  <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
top10<-head(VariableFeatures(seurat),10)
top10
##如果内存足够最好对所有基因进行中心化
scale.genes <-  rownames(seurat)
seurat <- ScaleData(seurat, features = scale.genes)
###已经是scale后的数据；scale后的数据一般以0为中心，正负都有
seurat <- RunPCA(seurat , features = VariableFeatures(object = seurat ))
plot1 <- DimPlot(seurat, reduction = "pca", group.by="orig.ident") 
plot1
plot2 <- ElbowPlot(seurat, ndims=20, reduction="pca") 
plot2
pc.num=1:10

seurat <- FindNeighbors(seurat, dims = pc.num)
seurat <- FindClusters(seurat, resolution = 0.5)
table(seurat@meta.data$seurat_clusters)
metadata <- seurat@meta.data
#UMAP
seurat <- RunUMAP(seurat, n.neighbors = 50,n.epochs=50,dims = pc.num)
#seurat <- RunUMAP(seurat, dims = pc.num)
embed_umap <- Embeddings(seurat, 'umap')
setwd('D:/2021scRNA-seq+3 BRAF-MUT/20210803Stromal')
write.csv(embed_umap,'embed_umap n.neighbors = 50,n.epochs=50.csv') 
plot2 = DimPlot(seurat,label = T, reduction = "umap") 
plot2
setwd('D:/2021scRNA-seq+3 BRAF-MUT')
ggsave("full cluster/UMAP20211114.pdf", plot = plot2, width =5.5, height = 4)

seurat@meta.data$groups<-seurat@meta.data$group

seurat<-subset(seurat,idents='Stromal cells')
save(seurat,file='seurat_stromal3.Rdat')

table(seurat@meta.data$louvain)

#通过直接建模在单细胞数据固有的均方差关系改进了以前的版本，并在实现FindVariableFeatures功能。V3默认选择2000个差异基因。这些将用于下游分析，例如PCA???
seurat  <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
top10<-head(VariableFeatures(seurat),10)
top10
##如果内存足够最好对所有基因进行中心化
scale.genes <-  rownames(seurat)
seurat <- ScaleData(seurat, features = scale.genes)
###已经是scale后的数据；scale后的数据一般以0为中心，正负都有
seurat <- RunPCA(seurat , features = VariableFeatures(object = seurat ))
plot1 <- DimPlot(seurat, reduction = "pca", group.by="orig.ident") 
plot2 <- ElbowPlot(seurat, ndims=20, reduction="pca") 
plot2
pc.num=1:8

seurat <- FindNeighbors(seurat, dims = pc.num)
seurat <- FindClusters(seurat, resolution = 0.5)
table(seurat@meta.data$seurat_clusters)
metadata <- seurat@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
write.csv(cell_cluster,'cell_cluster resolution0.2.csv',row.names = F)
#非线性降
seurat = RunTSNE(seurat, dims = pc.num)
embed_tsne <- Embeddings(seurat, 'tsne')
write.csv(embed_tsne,'T cluster/embed_tsne.csv')
plot1 = DimPlot(seurat, reduction = "tsne") 
library(ggplot2)
#UMAP
seurat <- RunUMAP(seurat, n.neighbors = 100,n.epochs=100,dims = pc.num)
#seurat <- RunUMAP(seurat, dims = pc.num)
embed_umap <- Embeddings(seurat, 'umap')
write.csv(embed_umap,'embed_umap n.neighbors = 50,n.epochs=50.csv') 
plot2 = DimPlot(seurat,label = T, reduction = "umap") 
plot2
ggsave("20210803Stromal/UMAP_Stromal.pdf", plot = plot2, width = 5, height = 4)
plot3 =FeaturePlot(seurat, features = c('MMP11'), min.cutoff = 0, max.cutoff = 1)
ggsave("20210803Stromal/MMP11 brafCRC.pdf", plot = plot3, width = 4.6, height = 4)
plot4 =FeaturePlot(seurat, features = c('PDGFRA'))
ggsave("20210803Stromal/PDGFRA brafCRC.pdf", plot = plot4, width = 4.6, height = 4)
setwd('D:/2021scRNA-seq+3 BRAF-MUT')
load('seurat_stromal3.Rdat')
plot3 =FeaturePlot(seurat, features = c('POSTN'), min.cutoff = 0, max.cutoff = 4)
ggsave("20210803Stromal/POSTN brafCRC.pdf", plot = plot3, width = 4.6, height = 4)



## Figure S10D ## 
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
load("Stromal cell_identify\\NC与SL EC的差异基因\\y.Rdat")
MSigDB.hr <- read.gmt("D:\\2018R\\clusterProfiler\\reactome.hallmark.v7.3.symbols.gmt")
resData <- read.csv("Stromal cell_identify\\NC与SL EC的差异基因\\NC与SL EC的差异基因全0.025.csv",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
geneList <- resData$avg_log2FC; names(geneList) <- rownames(resData)
geneList <- sort(geneList,decreasing = T)
gsea.val2 <- GSEA(geneList = geneList,
                  TERM2GENE = MSigDB.hr,
                  pvalueCutoff = 0.25,
                  seed = T,
                  verbose = F)
gsea.val2.df <- as.data.frame(gsea.val2)
write.csv(gsea.val2.df,"gsea20211211.csv",row.names = F,quote = F)
yd <- as.data.frame(y)

jama <- c("#3B4E55","#D69044","#44A0D5","#A94747","#81AF96")


id <- c("HALLMARK_OXIDATIVE_PHOSPHORYLATION","REACTOME_THE_CITRIC_ACID_TCA_CYCLE_AND_RESPIRATORY_ELECTRON_TRANSPORT","REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION",
        'HALLMARK_INFLAMMATORY_RESPONSE')
gseaplot2(x = y,
          color = brewer.pal(n = 4,name = "Set2"),
          #color = clust.col,
          subplots = 1:3,
          geneSetID = id)
ggsave("gseaplot up20211211-4.pdf", width = 5,height = 4.5)



