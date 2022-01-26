##Figure S8A-C ## 
seurat <- RunUMAP(seurat, n.neighbors = 50,n.epochs=50,dims = pc.num)
save(seurat,file='seurat_Myeloid.Rdat')
Idents(seurat) <- "seurat_clusters" 
plot2 = DimPlot(seurat,label = T, reduction = "umap")
ggsave("Myeloid cluster/UMAP n.neighbors = 50,n.epochs=50 20210407.pdf", plot = plot2, width = 6*0.75, height = 5*0.75)

plot3<-FeaturePlot(seurat, features = c('CD14','FCN1','S100A8','S100A9','IL1B','MRC1','CD163',
                                        'APOE','SELENOP','MAF','CD1C','FCER1A'))#作图panel 
ggsave("Myeloid cluster/Featureplot 20211129.pdf", plot = plot3, width = 12, height = 9)



##Figure S8D ## 
# Feature Plot to visulize gene expression and regulon activity -----------

path <- "i:/genomicdata/jiaoda/SCENIC/output/";setwd(path)
raw.path <- "i:/genomicdata/jiaoda/raw"
load(file.path(raw.path, "seurat_Myeloid.Rdata"))

# add regulon activity into seurat object
aucmtx <- read.csv("auc_mtx.csv")
colnames(aucmtx) <- substr(colnames(aucmtx),1,nchar(colnames(aucmtx))-3)
rownames(aucmtx) <- aucmtx[,1]
rownames(aucmtx) <- gsub("\\.","-",rownames(aucmtx))
info <- as.data.frame(aucmtx[,c("ETV5","HSF4")])
info <- info[colnames(seurat),]
colnames(info) <- paste0("AUCell_",colnames(info))
seurat <- AddMetaData(seurat, metadata = info)

# Do plot
pdf(file.path(path,"ETV5+HSF4-1220.pdf"),width = 411/100, height = 675/100)
plot_grid(FeaturePlot(seurat,features = "ETV5",cols = c("grey","#DA2D31"))+
            theme(panel.border = element_rect(colour = "black", fill=NA, size=1)),
          FeaturePlot(seurat, features = "AUCell_ETV5",cols = c("grey","#DA2D31"), min.cutoff = 0, max.cutoff = 0.2)+
            scale_color_gradient2(low="grey",mid="grey",high="#DA2D31",midpoint = 0.12, limits = c(0,0.2))+
            theme(panel.border = element_rect(colour = "black", fill=NA, size=1)),
          ncol = 1)

plot_grid(FeaturePlot(seurat,features = "HSF4",cols = c("grey","#DA2D31"), min.cutoff = 0, max.cutoff = 0.1)+
            theme(panel.border = element_rect(colour = "black", fill=NA, size=1)),
          FeaturePlot(seurat, features = "AUCell_HSF4",cols = c("grey","#DA2D31"), min.cutoff = 0, max.cutoff = 0.2)+
            scale_color_gradient2(low="grey",mid="grey",high="#DA2D31",midpoint = 0.03, limits = c(0,0.2))+
            theme(panel.border = element_rect(colour = "black", fill=NA, size=1)),
          ncol = 1)
invisible(dev.off())