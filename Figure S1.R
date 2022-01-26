##### Figure plotting #####
## Figure S1E ## 
setwd('D:/2021scRNA-seq')
load("seurat_full.Rdat")
Idents(seurat) <- "new.patient" 
seurat<- RenameIdents(seurat, 'NC1'='NC1','NC2'='NC2','NC3'='NC3','HP1'='HP1',
                      'HP2'='HP2','HP3'='HP3','HP4'='HP4','SSL1'='SSL1',
                      'SSL2'='SSL2','SSL3'='SSL3','SSL4'='SSL4','SSL5'='SSL5',
                      'SSLD'='SSLD','TSA1'='TSA1','TSA2'='TSA2','TSA3'='TSA3','TSA4'='TSA4',
                      'TSA5'='TSA5')
seurat$new.patient<-Idents(seurat)
plot1<-VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA"),pt.size=0,group.by = 'new.patient')
  ggsave("serrated adenoma/SL cluster/QC 20211215.pdf", plot = plot1, width = 10, height = 5*0.8)


## Figure S1F ## 
plot2<-DimPlot(seurat, reduction = "umap", group.by ='louvain',label = T)
plot2
ggsave("serrated adenoma/SL cluster/UMAP 20211215.pdf", plot = plot2, width = 7.5*0.8, height = 5)


## Figure S1G ## 
col1<-brewer.pal(12,"Paired")
col2<-brewer.pal(6,"Set2")
#colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(22)
#scales::show_col(colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(22))
plot3<-DimPlot(seurat, reduction = "umap", group.by ='new.patient',label = F,cols=c(col2,col1))
plot3
ggsave("serrated adenoma/SL cluster/UMAPpatient 20211215.pdf", plot = plot3, width = 7.5*0.8, height = 5)


## Figure S1H ## 
Idents(seurat) <- "annotated" 
diff.wilcox = FindAllMarkers(seurat,min.pct = 0.8)
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(all.markers, "serrated adenoma/SL cell_identify/diff_genes_wilcox.csv", row.names = F)
write.csv(top10, "serrated adenoma/SL cell_identify/top10_diff_genes_wilcox.csv", row.names = F)

DoHeatmap(seurat,features=top10$gene)+NoLegend()

library(ComplexHeatmap)
standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
  outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata>halfwidth]=halfwidth
    outdata[outdata<(-halfwidth)]= -halfwidth
  }
  return(outdata)
}


df.gene <- read.csv(file.path(work.path, "top10_diff_genes_wilcox(1).csv"))
load(file.path(raw.path, "seurat_full.Rdata"))
# df.gene <- split(x = df.gene$gene, f = df.gene$cluster)
plot.gene <- df.gene$gene
mtx <- AverageExpression(seurat, features = plot.gene, group.by = "annotated")[[1]]
plot.data <- standarize.fun(indata = mtx, halfwidth = 2)
# plot.data <- standarize.fun(indata = mtx, halfwidth = 40, centerFlag = F, scaleFlag = F)

Epithelial.col <- "#FF8214"
Stromal.col <- "#AA40FC"
Myeloid.col <- "#DD4D4D"
Tcells.col <- "#925F54"
Bcells.col <- "#1E77B4"
Mast.col <- "#279E68"
annCol <- data.frame(celltypes = factor(x = colnames(plot.data), levels = colnames(plot.data)),
                     row.names = colnames(plot.data),
                     stringsAsFactors = F)
# annColors <- list()
# annColors[["celltypes"]] = c(
#   "T cells" = Tcells.col,
#   "B cells" = Bcells.col,
#   "Epithelial" = Epithelial.col,
#   "Stromal cells" = Stromal.col,
#   "Myeloids" = Myeloid.col,
#   "Mast cells" = Mast.col
# ) 
# 
# pdf(file.path(work.path, "HeatmapForAnotated-1216.pdf"), width = 565/100, height = 420/100,)
# pheatmap(plot.data,
#          cluster_cols = F, cluster_rows = F,
#          # border_color = "white",
#          # color = colorRampPalette(colors = c("#6353A3","white","#9A0038"))(100),
#          color = viridis::inferno(64),
#          annotation_col = annCol,
#          annotation_colors = annColors,
#          cellwidth = 12, cellheight = 4,
#          # gaps_row = seq(10, 50, 10),
#          fontsize = 3,
#          angle_col = "45",
#          treeheight_row = 20)
# invisible(dev.off())


ha = HeatmapAnnotation(
  celltypes = colnames(plot.data),
  col = list(celltypes = c(
    "T cells" = Tcells.col,
    "B cells" = Bcells.col,
    "Epithelial" = Epithelial.col,
    "Stromal cells" = Stromal.col,
    "Myeloids" = Myeloid.col,
    "Mast cells" = Mast.col
  )),
  simple_anno_size = unit(4, "points"),
  show_annotation_name = F
)

pdf(file.path(work.path, "HeatmapForAnotated-1216a.pdf"), width = 565/100, height = 420/100,)
Heatmap(matrix = as.matrix(plot.data),
        cluster_columns = F, cluster_rows = F,
        rect_gp = gpar(col = "#999999", lwd = 1),
        show_row_names = T, row_names_gp = gpar(fontsize = 3),
        show_column_names = T,column_names_gp = gpar(fontsize = 8),
        col = viridis::inferno(64),
        top_annotation = ha,
        width = unit(12, "points")*dim(plot.data)[2],
        height = unit(4, "points")*dim(plot.data)[1],
        name = "scaled log-exp")
invisible(dev.off())
  
  
sessionInfo()