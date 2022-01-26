##### Figure plotting #####
## Figure 6A ## 
setwd('D:/2021scRNA-seq')
load("seurat_Myeloid.Rdat")
seurat <- RunUMAP(seurat, n.neighbors = 50,n.epochs=50,dims = pc.num)
save(seurat,file='seurat_Myeloid.Rdat')
Idents(seurat) <- "seurat_clusters" 
plot2 = DimPlot(seurat,label = T, reduction = "umap") 
load("Myeloid cluster\\Myeloid.col.RData")
cell.col
seurat@meta.data$celltypes<-seurat@meta.data$seurat_clusters
Idents(seurat) <- "celltypes" 
seurat<- RenameIdents(seurat, '0'='Anti-inflammatory macrophages','2'='Anti-inflammatory macrophages',
                      '7'='FCN1+ inflammatory macrophages','10'='Anti-inflammatory macrophages','6'='Conventional DCs','3'='Conventional DCs',
                      '4'='Conventional DCs','5'='Plasmacytoid DCs','8'='Anti-inflammatory macrophages',
                      '9'='Anti-inflammatory macrophages','1'='Anti-inflammatory macrophages')
seurat$celltypes<-Idents(seurat)
plot2<-DimPlot(seurat, reduction = "umap",  group.by ='celltypes',cols=cell.col,label = F)
ggsave("Myeloid cluster/UMAP 20210912.pdf", plot = plot2, width = 8, height = 5)
ggsave("Myeloid cluster/UMAP n.neighbors = 50,n.epochs=50 20210407.png", plot = plot2, width = 8, height = 7)



load("Myeloid_metadata.Rdat")
dat <- Myeloid_metadata[,c("group","celltypes")]
dat$celltypes <- as.character(dat$celltypes)
table(dat$celltypes)
dat$group <- factor(dat$group,levels = c("NC","HP","SSL","SSLD","TSA"))
dat$celltypes <- factor(dat$celltypes,levels = c("Anti-inflammatory macrophages",
                                                 "FCN1+ inflammatory macrophages",
                                                 "Conventional DCs",
                                                 "Plasmacytoid DCs"))


dat.bar <- as.data.frame(table(dat$group,dat$celltypes))
dat.bar <- dat.bar%>% group_by(Var1) %>% 
  mutate(Pct = Freq / sum(Freq) * 100) %>% as.data.frame()
colnames(dat.bar) <- c("group","celltype","count","pct")

mycol <- brewer.pal(4,"Set2")
cell.col <- c("Anti-inflammatory macrophages" = mycol[1],
              "FCN1+ inflammatory macrophages" = mycol[2],
              "Conventional DCs" = mycol[3],
              "Plasmacytoid DCs" = mycol[4])
save(cell.col,file = "Myeloid.col.RData")

# left barplot
p.left <- ggplot(dat.bar, aes(x = group, y = count,fill = celltype)) +
  scale_fill_manual(values  = cell.col) +
  geom_bar(stat = 'identity') +
  scale_x_discrete(name = "",position = "bottom") +
  ylab("Myeloid cell number") + 
  theme_bw() +
  theme(axis.line.y = element_blank(),
        axis.line.x = element_line(size = 0.6),
        axis.ticks.y = element_line(size = 0.6),
        axis.title = element_text(size = 12,colour = "black"),
        axis.text.y = element_text(colour = "black",size = 10),
        axis.title.x = element_text(vjust = -0.3,size = 12),
        axis.text.x = element_text(size = 10, color = "black",angle = 45,hjust = 1),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.3, -1.7, 0.3, 0.3), "lines"),
        legend.title = element_blank())

p.left


p.right <- ggplot(dat.bar, aes(x = group, y = pct, fill=celltype)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values  = cell.col) +
  scale_x_discrete(name = "") +
  theme_bw() +
  ylab("Myeloid cell proportion (%)") + 
  theme(axis.line.y = element_blank(),
        axis.line.x = element_line(size = 0.6),
        axis.ticks.y = element_line(size = 0.6),
        axis.title = element_text(size = 12,colour = "black"),
        axis.title.x = element_text(vjust = -0.3,size = 12),
        axis.text.x = element_text(size = 10, color = "black",angle = 45,hjust = 1),
        axis.text.y = element_text(size = 10, color = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.3, -1.7, 0.3, 0.3), "lines"),
        legend.title = element_blank())
p.right

pal <- p.left + p.right +
  plot_layout(widths = c(6,6), guides = 'collect') & theme(legend.position = 'right',legend.key.size = unit(0.4, 'cm'))
pal
ggsave(filename = "distribution of Myeloid cells.pdf", width = 5.8,height = 3)


load("Myeloid_metadata.Rdat")
dat <- Myeloid_metadata[,c("patient","group","celltypes")]
dat$celltypes <- factor(dat$celltypes,levels = c("Anti-inflammatory macrophages",
                                                 "FCN1+ inflammatory macrophages",
                                                 "Conventional DCs",
                                                 "Plasmacytoid DCs"))
dat$new.patient <- as.character(dat$patient)
table(dat$new.patient)
dat[which(dat$new.patient == "SSL3"),"new.patient"] <- "HP2"
dat[which(dat$new.patient == "SSL6"),"new.patient"] <- "SSL3"
table(dat$new.patient)
dat$new.patient <- factor(dat$new.patient,levels = rev(c("NC1","NC2","NC3","HP1","HP2","HP3","HP4","SSL1","SSL2","SSL3","SSL4","SSL5","SSLD","TSA1","TSA2","TSA3","TSA4","TSA5")))
# mycol <- rev(viridis(10))
# mycol <- rev(inferno(10))
# mycol <- brewer.pal(10,"Spectral")

dat.bar <- as.data.frame(table(dat$new.patient,dat$celltypes))
dat.bar <- dat.bar%>% group_by(Var1) %>% 
  mutate(Pct = Freq / sum(Freq) * 100) %>% as.data.frame()
colnames(dat.bar) <- c("patient","celltype","count","pct")

# left barplot
p.left <- ggplot(dat.bar, aes(x = patient, y = count,fill = celltype)) +
  scale_fill_manual(values  = cell.col) +
  geom_bar(stat = 'identity') +
  scale_x_discrete(name = "",position = "top") +
  #theme_bw() +
  theme(axis.line.y = element_blank(),
        axis.line.x = element_line(size = 0.6),
        axis.ticks.y = element_line(size = 0.6),
        axis.text.y = element_text(colour = "black",size = 10),
        axis.title.x = element_text(vjust = -0.3,size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.3, -1.7, 0.3, 0.3), "lines"),
        legend.title = element_blank()) +
  coord_flip() +
  scale_y_reverse(expand = c(0.01,0),
                  name = "Myeloid cell number", position = "right")
p.left


p.right <- ggplot(dat.bar, aes(x = patient, y = pct, fill=celltype)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values  = cell.col) +
  scale_x_discrete(name = "") +
  #theme_bw() +
  theme(axis.line.y = element_blank(),
        axis.line.x = element_line(size = 0.6),
        axis.ticks.y = element_line(size = 0.6),
        axis.text.y = element_blank(),
        axis.title.x = element_text(vjust = -0.3,size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.3, -1.7, 0.3, 0.3), "lines"),
        legend.title = element_blank()) +
  coord_flip() +
  scale_y_continuous(expand = c(0.01,0),
                     name = "Myeloid cell proportion (%)", position = "right")

p.right

pal <- p.left + p.right +
  plot_layout(widths = c(7,7), guides = 'collect') & theme(legend.position = 'right',legend.key.size = unit(0.5, 'cm'))
pal
ggsave(filename = "distribution of Myeloid cells2.pdf", width = 6.8,height = 3)



## Figure 6B ## 
# set working path
workdir <- "F:/IGBMC/externalproject/ExternalProject/ZYJ/5"; setwd(workdir)

# set colors
blue   <- "#5bc0eb"
yellow <- "#fde74c"
green  <- "#9bc53d"
red    <- "#f25f5c"
purple <- "#531f7a"
grey   <- "#8693ab"
orange <- "#fa7921"
white  <- "#f2d7ee"
darkred   <- "#F2042C"
lightred  <- "#FF7FBF"
lightblue <- "#B2EBFF"
darkblue  <- "#1d00ff"
cherry    <- "#700353"
lightgrey <- "#dcddde"
nake <- "#F8C364"
gold <- "#ECE700"
cyan <- "#00B3D0"
sun  <- "#E53435"
peach  <- "#E43889"
violet <- "#89439B"
soil   <- "#EC7D21"
lightgreen <- "#54B642"
darkblue   <- "#21498D"
darkgreen  <- "#009047"
brown      <- "#874118"
seagreen   <- "#008B8A"
jco <- c("#2874C5","#EABF00","#868686","#C6524A","#80A7DE")
heatmap.BlWtRd <- c("#6699CC","white","#FF3C38")
heatmap.YlGnPe <- c("#440259","#345F8C","#228C8A","#78CE51","#FAE71F")
heatmap.GrWtRd <- c("#2b2d42","#8d99ae","#edf2f4","#ef233c","#d90429")
heatmap.L.BlYlRd <- c("#4281a4","#9cafb7","#ead2ac","#e6b89c","#fe938c")
heatmap.BlBkRd <- c("#54FEFF","#32ABAA","#125456","#000000","#510000","#A20000","#F30000")
heatmap.fancy <- c("#10040A", "#2A0B35", "#4D155B", "#73215B", "#9C3558", "#C34D44", "#E07038", "#F2981C", "#F2CA51", "#FAF6A3")

# customized function
standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
  outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata>halfwidth]=halfwidth
    outdata[outdata<(-halfwidth)]= -halfwidth
  }
  return(outdata)
}

gmt2list <- function(annofile){
  if (!file.exists(annofile)) {
    stop("There is no such gmt file.")
  }
  
  if (tools::file_ext(annofile) == "xz") {
    annofile <- xzfile(annofile)
    x <- scan(annofile, what="", sep="\n", quiet=TRUE)
    close(annofile)
  } else if (tools::file_ext(annofile) == "gmt") {
    x <- scan(annofile, what="", sep="\n", quiet=TRUE)
  } else {
    stop ("Only gmt and gmt.xz are accepted for gmt2list")
  }
  
  y <- strsplit(x, "\t")
  names(y) <- sapply(y, `[[`, 1)
  
  annoList <- lapply(y, `[`, c(-1,-2))
}

# single sample enrichment
load("seurat_Myeloid.Rdat")
MSigDB.hm.rt <- gmt2list("reactome.hallmark.v7.3.symbols.gmt")
MSigDB.hm.rt <- sapply(MSigDB.hm.rt,function(x) setdiff(x,""))
MSigDB.imm <- gmt2list("c7.all.v7.4.symbols.gmt")
expr <- seurat@assays$RNA@scale.data
hmrt.ssgsea <- gsva(as.matrix(expr),
                    MSigDB.hm.rt,
                    method = "ssgsea")
hmrt.gsva <- gsva(as.matrix(expr),
                  MSigDB.hm.rt,
                  method = "gsva")
imm.ssgsea <- gsva(as.matrix(expr),
                   MSigDB.imm,
                   method = "ssgsea")
imm.gsva <- gsva(as.matrix(expr),
                 MSigDB.imm,
                 method = "gsva")
write.table(hmrt.ssgsea,file = "hmrt.ssgsea.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(hmrt.gsva,file = "hmrt.gsva.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(imm.ssgsea,file = "imm.ssgsea.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(imm.gsva,file = "imm.gsva.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

pseudo.moic.res                 <- list("clust.res" = dat,
                                        "mo.method" = "Myeloid")

# make pseudo samID
pseudo.moic.res$clust.res$samID <- rownames(pseudo.moic.res$clust.res)
pseudo.moic.res$clust.res$clust <- sapply(pseudo.moic.res$clust.res$celltypes,
                                          switch,
                                          "Anti-inflammatory macrophages" = 1,
                                          "FCN1+ inflammatory macrophages" = 2, 
                                          "Conventional DCs" = 3, 
                                          "Plasmacytoid DCs" = 4) 
runDEA(dea.method = "limma",
       expr       = hmrt.ssgsea, # normalized expression data
       moic.res   = pseudo.moic.res,
       overwt     = TRUE,
       prefix     = "hmrt.ssgsea")
marker.up <- runMarker(moic.res      = pseudo.moic.res,
                       dea.method    = "limma", # name of DEA method
                       prefix        = "hmrt.ssgsea", # MUST be the same of argument in runDEA()
                       dat.path      = getwd(), # path of DEA files
                       res.path      = getwd(), # path to save marker files
                       p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                       p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                       dirct         = "up", # direction of dysregulation in expression
                       norm.expr     = hmrt.ssgsea,
                       n.marker      = 30, # number of biomarkers for each subtype
                       doplot        = TRUE) # no heatmap

runDEA(dea.method = "limma",
       expr       = hmrt.ssgsea, # normalized expression data
       moic.res   = pseudo.moic.res,
       overwt     = TRUE,
       prefix     = "hmrt.gsva")
marker.up <- runMarker(moic.res      = pseudo.moic.res,
                       dea.method    = "limma", # name of DEA method
                       prefix        = "hmrt.gsva", # MUST be the same of argument in runDEA()
                       dat.path      = getwd(), # path of DEA files
                       res.path      = getwd(), # path to save marker files
                       p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                       p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                       dirct         = "up", # direction of dysregulation in expression
                       norm.expr     = hmrt.gsva,
                       n.marker      = 30, # number of biomarkers for each subtype
                       doplot        = TRUE) # no heatmap

runDEA(dea.method = "limma",
       expr       = imm.ssgsea, # normalized expression data
       moic.res   = pseudo.moic.res,
       overwt     = TRUE,
       prefix     = "imm.ssgsea")
marker.up <- runMarker(moic.res      = pseudo.moic.res,
                       dea.method    = "limma", # name of DEA method
                       prefix        = "imm.ssgsea", # MUST be the same of argument in runDEA()
                       dat.path      = getwd(), # path of DEA files
                       res.path      = getwd(), # path to save marker files
                       p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                       p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                       dirct         = "up", # direction of dysregulation in expression
                       norm.expr     = imm.ssgsea,
                       n.marker      = 30, # number of biomarkers for each subtype
                       doplot        = TRUE) # no heatmap

runDEA(dea.method = "limma",
       expr       = imm.ssgsea, # normalized expression data
       moic.res   = pseudo.moic.res,
       overwt     = TRUE,
       prefix     = "imm.gsva")
marker.up <- runMarker(moic.res      = pseudo.moic.res,
                       dea.method    = "limma", # name of DEA method
                       prefix        = "imm.gsva", # MUST be the same of argument in runDEA()
                       dat.path      = getwd(), # path of DEA files
                       res.path      = getwd(), # path to save marker files
                       p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                       p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                       dirct         = "up", # direction of dysregulation in expression
                       norm.expr     = imm.gsva,
                       n.marker      = 30, # number of biomarkers for each subtype
                       doplot        = TRUE) # no heatmap

tmp <- cbind.data.frame(t(hmrt.ssgsea),celltypes = dat$celltypes)
tmp <- as.data.frame(t(apply(tmp[,setdiff(colnames(tmp), "celltypes")], 2, function(x) tapply(x, INDEX=factor(tmp$celltypes), FUN=mean, na.rm=TRUE))))
write.table(tmp,file = "hmrt.ssgsea.mean.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

tmp <- cbind.data.frame(t(hmrt.gsva),celltypes = dat$celltypes)
tmp <- as.data.frame(t(apply(tmp[,setdiff(colnames(tmp), "celltypes")], 2, function(x) tapply(x, INDEX=factor(tmp$celltypes), FUN=mean, na.rm=TRUE))))
write.table(tmp,file = "hmrt.gsva.mean.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

tmp <- cbind.data.frame(t(imm.ssgsea),celltypes = dat$celltypes)
tmp <- as.data.frame(t(apply(tmp[,setdiff(colnames(tmp), "celltypes")], 2, function(x) tapply(x, INDEX=factor(tmp$celltypes), FUN=mean, na.rm=TRUE))))
write.table(tmp,file = "imm.ssgsea.mean.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

tmp <- cbind.data.frame(t(imm.gsva),celltypes = dat$celltypes)
tmp <- as.data.frame(t(apply(tmp[,setdiff(colnames(tmp), "celltypes")], 2, function(x) tapply(x, INDEX=factor(tmp$celltypes), FUN=mean, na.rm=TRUE))))
write.table(tmp,file = "imm.gsva.mean.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

# generate heatmap
indata <- read.csv("hmrt.gsva.mean.heatmap.csv",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
indata <- standarize.fun(indata,halfwidth = 1)
rownames(indata) <- gsub("HALLMARK_","",rownames(indata))
rownames(indata) <- gsub("REACTOME_","",rownames(indata))

annCol <- data.frame(celltype = unique(dat$celltypes),
                     row.names = unique(dat$celltypes))
annColors <- list(celltype = cell.col)
hm <- pheatmap(as.matrix(indata[c(1,3,4,5,6,11,2,7,9,8,10),]),
               #color = viridis(64),
               color = NMF:::ccRamp(c("#01665E","white","#94580F"),64),
               annotation_col = annCol[colnames(indata),,drop = F],
               annotation_colors = annColors,
               cluster_rows = F,
               cluster_cols = F,
               show_rownames = T,
               show_colnames = T,
               border_color = "white",
               cellheight = 12,
               cellwidth = 30,
               treeheight_row = 10,
               treeheight_col = 10)
pdf("hmrt.gsva.mean.heatmap.pdf",width = 15,height = 6)
draw(hm, annotation_legend_side = "bottom",heatmap_legend_side = "left")
invisible(dev.off())


## Figure 6C ## 
# load R package
library(ggplot2)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(viridis)
library(tidyr)
library(Seurat)
library(cowplot)
library(ComplexHeatmap)
library(gplots)
library(GSVA)
library(MOVICS)



# generate another heatmap
genes <- c("CD14","CD68","FCN1","VCAN","S100A8","S100A9","IL1B",
           "MRC1","CD163","C1QA","APOE",
           "SELENOP","MAF",
           "ITGAX","CLEC9A","CD1C","FCER1A","MZB1")
subexpr <- expr[genes,rownames(Myeloid_metadata)]

indata <- as.data.frame(t(subexpr))
indata$cluster <- Myeloid_metadata$celltypes
indata <- apply(indata[,setdiff(colnames(indata), "cluster")], 2, function(x) tapply(x, INDEX=factor(indata$cluster), FUN=mean, na.rm=TRUE))
indata <- as.data.frame(t(indata))
colnames(indata)
indata <- indata[,c(1,3,2,4)]

plotdata <- standarize.fun(indata,halfwidth = 2)
range(plotdata)
annCol <- data.frame(celltype = c("Anti-inflammatory macrophages","FCN1+ inflammatory macrophages","Conventional DCs","Plasmacytoid DCs"),
                     row.names = c("Anti-inflammatory macrophages","FCN1+ inflammatory macrophages","Conventional DCs","Plasmacytoid DCs"))
annColors <- list("celltype" = cell.col)
hm <- pheatmap(plotdata,
               border_color = "black",
               color = NMF:::ccRamp(c("#E6EAF7","#B6D1E8","#498EB9","#204F8D"),64),
               cellwidth = 20,
               cellheight = 12,
               cluster_cols = F,
               cluster_rows = F,
               show_rownames = T,
               show_colnames = T,
               annotation_col = annCol,
               annotation_colors = annColors,
               name = "Myeloid",
               #gaps_col = c(6,8,12,13,16),
               gaps_row = c(7,13))
pdf("heatmap of Myeloid cells using specific markers.pdf", width = 7,height = 6)
draw(hm)
invisible(dev.off())

save.image(file = "Myeloid.RData")



## Figure 6D ## 

# SCENIC_heatmap_myeloid --------------------------------------------------

avg_celltype <- read.csv(file = "output/avg_celltype.csv", header = T, row.names = 1)
colnames(avg_celltype)[29]<-"delta Gamma T cells"
typelist<-table(cell_info[cell_info$type=="mye",]$celltypes)
typelist<-rownames(typelist)[typelist>0]
reg<-sapply(typelist,function(i){
  gettopreg(i,10)
})
reg<-as.vector(reg)
reg<-reg[!duplicated(reg)]
colnames(avg_celltype)[1:4]<-typelist
data<-avg_celltype[reg,typelist]
scaled.data<-t(scale(t(data)))

mye_order<-c("Anti-inflammatory macrophages","FCN1+ inflammatory macrophages","Conventional DCs","Plasmacytoid DCs")
order <- mye_order

scaled.data<-scaled.data[,order]
plot.data<-data.frame(scaled.data[,order],
                      "rank"=c(apply(scaled.data,1,function(data){match(max(data),data)})),
                      "index"=c(apply(scaled.data,1,function(data){match(max(data),data)})))
plot.data<-split(plot.data,plot.data$rank)
plot.data<-lapply(plot.data,function(x){x<-dplyr::arrange(x,-x[,unique(x$index)])})
plot.data<-do.call(rbind,plot.data)
plot.data<-plot.data[,-c(length(plot.data)-1,length(plot.data))]
rownames(plot.data)<-t(data.frame(strsplit(rownames(plot.data),"\\.")))[,2]
colnames(plot.data)<-order



anncolors[["celltype"]]<-c("Anti-inflammatory macrophages" = "#52BA9A",
                           "Conventional DCs" = "#7D92C4",
                           "FCN1+ inflammatory macrophages" = "#FC8254",
                           "Plasmacytoid DCs" = "#E581BE")

pdf("output/myeloid.pdf",width = 8, height = 10)
pheatmap(plot.data,
         border_color = "white",
         color = colorRampPalette(colors = c("#75A3D0","white","red"))(100),
         cluster_rows = F,
         cluster_cols = F,
         annotation_col = anncol[colnames(plot.data),,drop=F],
         annotation_colors = anncolors,
         cellwidth = 30, cellheight = 12)
invisible(dev.off())



## Figure 6E ## 
load("seurat_B.Rdat")
seurat  <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 4000)
top10<-head(VariableFeatures(seurat),10)
top10
scale.genes <-  rownames(seurat)
seurat <- ScaleData(seurat, features = scale.genes)
seurat <- RunPCA(seurat , features = VariableFeatures(object = seurat ))
plot2 <- ElbowPlot(seurat, ndims=20, reduction="pca") 
plot2
pc.num=1:6
seurat <- FindNeighbors(seurat, dims = pc.num)
seurat <- FindClusters(seurat, resolution = 0.5)
table(seurat@meta.data$group)
table(seurat@meta.data$seurat_clusters,seurat@meta.data$group)
#mat<-as.data.frame(table(seurat@meta.data$seurat_clusters,seurat@meta.data$group))
#write.csv(mat,'B cluster/各亚型比例.csv',row.names = F)
metadata <- seurat@meta.data
colnames(metadata)
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
write.csv(cell_cluster,'B cluster/cell_cluster nfeatures4000 20210327.csv',row.names = F)

seurat <- RunUMAP(seurat, n.neighbors = 50,n.epochs=50,dims = pc.num)
embed_umap <- Embeddings(seurat, 'umap')
write.csv(embed_umap,'B cluster/embed_umap n.neighbors = 50,n.epochs=50 20210327.csv') 
plot2 = DimPlot(seurat,label = T, group.by ='seurat_clusters', reduction = "umap") 
#DimPlot改颜色
load("B cluster\\B.col.Rdata")
cell.col
plot2<-DimPlot(seurat, reduction = "umap",cols=c( "#8DD3C7", "#FFFFB3","#BEBADA",  "#FB8072","#80B1D3"),
               group.by ='celltypes',label = T)
plot2
ggsave("B cluster/UMAP 20211129.pdf", plot = plot2, width = 6, height = 5)
ggsave("B cluster/UMAP n.neighbors = 50,n.epochs=50 20210327.png", plot = plot2, width = 8, height = 7)



dat <- B_metadata[,c("patient","group","celltypes")]
colnames(dat)[3] <- "celltypes"
dat$celltypes <- as.character(dat$celltypes)
dat[grep("ve",dat$celltypes),"celltypes"] <- "Naive mature B cells"
#dat[which(dat$celltypes == "CD8+ Memory"),"celltypes"] <- "CD8+ Effecor-2"
table(dat$celltypes)
dat$celltypes <- factor(dat$celltypes,levels = c("IgA+ plasma cells",
                                                 "Naive mature B cells",
                                                 #"CD8+ Memory",
                                                 "Memory B cells",
                                                 "Plasmablasts",
                                                 "Doublets (T cell)"))
dat$new.patient <- as.character(dat$patient)
table(dat$new.patient)
dat[which(dat$new.patient == "SSL3"),"new.patient"] <- "HP2"
dat[which(dat$new.patient == "SSL6"),"new.patient"] <- "SSL3"
table(dat$new.patient)
dat$new.patient <- factor(dat$new.patient,levels = rev(c("NC1","NC2","NC3","HP1","HP2","HP3","HP4","SSL1","SSL2","SSL3","SSL4","SSL5","SSLD","TSA1","TSA2","TSA3","TSA4","TSA5")))
# mycol <- rev(viridis(10))
# mycol <- rev(inferno(10))
# mycol <- brewer.pal(10,"Spectral")

dat.bar <- as.data.frame(table(dat$new.patient,dat$celltypes))
dat.bar <- dat.bar%>% group_by(Var1) %>% 
  mutate(Pct = Freq / sum(Freq) * 100) %>% as.data.frame()
colnames(dat.bar) <- c("patient","celltype","count","pct")

# left barplot
p.left <- ggplot(dat.bar, aes(x = patient, y = count,fill = celltype)) +
  scale_fill_manual(values  = cell.col) +
  geom_bar(stat = 'identity') +
  scale_x_discrete(name = "",position = "top") +
  #theme_bw() +
  theme(axis.line.y = element_blank(),
        axis.line.x = element_line(size = 0.6),
        axis.ticks.y = element_line(size = 0.6),
        axis.text.y = element_text(colour = "black",size = 10),
        axis.title.x = element_text(vjust = -0.3,size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.3, -1.7, 0.3, 0.3), "lines"),
        legend.title = element_blank()) +
  coord_flip() +
  scale_y_reverse(expand = c(0.01,0),
                  name = "B cell number", position = "right")
p.left


p.right <- ggplot(dat.bar, aes(x = patient, y = pct, fill=celltype)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values  = cell.col) +
  scale_x_discrete(name = "") +
  #theme_bw() +
  theme(axis.line.y = element_blank(),
        axis.line.x = element_line(size = 0.6),
        axis.ticks.y = element_line(size = 0.6),
        axis.text.y = element_blank(),
        axis.title.x = element_text(vjust = -0.3,size = 12),
        axis.text.x = element_text(size = 10, color = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.3, -1.7, 0.3, 0.3), "lines"),
        legend.title = element_blank()) +
  coord_flip() +
  scale_y_continuous(expand = c(0.01,0),
                     name = "B cell proportion (%)", position = "right")

p.right

# middle label
# pp <- ggplot() +
#   geom_text(data = dat.bar,
#              aes(label = patient, x = patient),
#              y = 0.5,
#              color = "black",
#              size = 0.8*11/.pt, # match font size to theme
#              hjust = 0.5, vjust = 0.5) +
#   theme_minimal()+
#   theme(axis.line.y =element_blank(),
#         axis.ticks.y =element_blank(),
#         axis.text.y =element_blank(),
#         axis.title.y =element_blank(),
#         axis.title.x =element_blank(),
#         plot.margin = unit(c(0.3, 0, 0.3, 0), "lines")
#   ) +
#   #guides(fill = FALSE) +
#   coord_flip() +
#   scale_y_reverse()
# pp

# combine figures
# pal <- p.left + pp + p.right +
#   plot_layout(widths = c(7,1,7), guides = 'collect') & theme(legend.position = 'right',legend.key.size = unit(0.5, 'cm'))
# pal

pal <- p.left + p.right +
  plot_layout(widths = c(7,7), guides = 'collect') & theme(legend.position = 'right',legend.key.size = unit(0.5, 'cm'))
pal
ggsave(filename = "B cluster//distribution of B cells2 20210829.pdf", width = 6,height = 3)



## Figure 6F ## 
plot4<-VlnPlot(Bcell, features =c('CXCR5','TNF','AIM2'), pt.size=0, group.by="groups", ncol=3)
ggsave("B cell_identify/VlnplotB 20210912.pdf", plot = plot4, width = 8, height = 3)
plot3<-FeaturePlot(seurat, features = c('CD79A','IGHA1','CD19','CD27','IGHD','CD38','SDC1','MS4A1'),cols=c('#D3D3D3','#1F7146'))#作图panel B
ggsave("B cell_identify/Featureplot 20210713.pdf", plot = plot3, width = 10, height = 8)
plot3<-FeaturePlot(seurat, features = c('CD79A','IGHA1','CD19','CD27','IGHD','CD38','SDC1','MS4A1'),cols=c('#D3D3D3','#1F7146'), min.cutoff = 0, max.cutoff = 2)#作图panel B
ggsave("B cell_identify/Featureplot 20211118.pdf", plot = plot3, width = 10, height = 8)


## Figure 6G ## 
## read output of CellChat
B_mye <- readRDS(paste0(path, "B_mye.rds"))
B_CD8T <- readRDS(paste0(path, "B_CD8T.rds"))
B_CD8T_data <- lapply(B_CD8T, subsetCommunication)
df.net <- subset(B_CD8T_data$NC, interaction_name%in%interlist & source == "IgA+ plasma cells")

interlist <- c("MIF_CD74_CXCR4","MDK_ITGA4_ITGB1","MDK_NCL",
               "GAS6_MERTK","CLEC2B_KLRB1","CD69_KLRB1",
               "ICAM2_ITGAL_ITGB2","NAMPT_ITGA5_ITGB1",
               "SEMA4D_PLXNB2","SEMA4A_PLXNB2")

## set function to extract value from cellchat object
GetInteractionScore <- function(interlist, group1, group2, group){
  gene <- lapply(interlist, FUN = function(x){unlist(strsplit(x,"_"))})
  score <- lapply(gene, FUN = function(g){
    sapply(group1, FUN = function(x){
      sapply(group2, FUN = function(y){
        source_expr <- expr[g[1],cell_info$celltypes == x & cell_info$group %in% group]
        target_expr <- expr[g[-1],cell_info$celltypes == y & cell_info$group %in% group]
        # value <- mean(c(as.vector(source_expr),as.vector(target_expr)))
        value <- (mean(source_expr)+mean(target_expr)*(length(g)-1))/length(g)
      })
    })
  })
  score <- do.call(cbind,score)
  colnames(score) =  interlist
  return(score)
}

turn <- function(mtx, name){
  record <- lapply(rownames(mtx), FUN = function(x){
    lapply(colnames(mtx), FUN = function(y){
      c(x,y,mtx[x,y])
    })
  })
  record <- lapply(record, FUN = function(x){do.call(rbind,x)})
  record <- do.call(rbind, record)
  record <- as.data.frame(record)
  colnames(record) <- name
  return(record)
}
FetchPvalue <- function(data,array){
  index <- dimnames(array)
  Pvalue <- array[match("IgA+ plasma cells", index[[1]]),
                  match(data$target, index[[2]]),
                  match(data$interaction_name1, index[[3]])]
  return(Pvalue)
}

## write function preparing for plotdata
DataToPlot <- function(group,division){
  score <- GetInteractionScore(interlist, "IgA+ plasma cells", c(CD8T,mye), group)
  xorder <- rownames(score)
  yorder <- interlist
  plotdata <- turn(score, c("target","interaction_name","score"))
  plotdata$interaction_name1 <- plotdata$interaction_name
  plotdata$interaction_name1[plotdata$interaction_name1=="CD69_KLRB1"] = "CLEC2C_KLRB1"
  plotdata$target <- factor(plotdata$target, levels = xorder)
  plotdata$interaction_name <- factor(plotdata$interaction_name, levels = yorder)
  plotdata$score <- as.numeric(plotdata$score)
  plotdata$Pvalue1 <- sapply(rownames(plotdata), FUN = function(x){
    FetchPvalue(plotdata[x,],B_mye[[division]]@net$pval)
  })
  plotdata$Pvalue2 <- sapply(rownames(plotdata), FUN = function(x){
    FetchPvalue(plotdata[x,],B_CD8T[[division]]@net$pval)
  })
  plotdata$Pvalue <- paste0(plotdata$Pvalue1,plotdata$Pvalue2)
  plotdata$Pvalue <- as.numeric(gsub("NA","",plotdata$Pvalue))
  plotdata$P <- cut(plotdata$Pvalue,c(1,0.05,0.01,0),include.lowest = T)
  plotdata$P <- factor(plotdata$P,levels = c("(0.05,1]","(0.01,0.05]","[0,0.01]"))
  plotdata$score[plotdata$score>5]=5
  plotdata <- plotdata[!is.na(plotdata$P),]
}
TMdata <- DataToPlot(c("HP","SSL","SSLD","TSA"),"TM")
NCdata <- DataToPlot("NC","NC")


## output plot
NC <- ggplot(NCdata, aes(x = target, y = interaction_name, col = score, size = P))+
  geom_point()+
  scale_colour_gradientn(colors = colorRampPalette(c("darkblue","yellow","red"))(99), na.value = "white") +
  scale_y_discrete(limits = rev)+
  theme_linedraw() + theme(panel.grid.major = element_blank())+
  theme(axis.text.y = element_text(color = "#70B7E6", size = 10), axis.title = element_blank(),
        axis.text.x = element_text(color = "#70B7E6", size = 10, hjust = 1, vjust = 0.5, angle = 90))+
  geom_vline(xintercept = 6.5)

TM <- ggplot(TMdata, aes(x = target, y = interaction_name, col = score, size = P))+
  geom_point()+
  scale_colour_gradientn(colors = colorRampPalette(c("darkblue","yellow","red"))(99), na.value = "white") +
  scale_y_discrete(limits = rev)+
  theme_linedraw() + theme(panel.grid.major = element_blank())+
  theme(axis.text.y = element_text(color = "#70B7E6", size = 10), axis.title = element_blank(),
        axis.text.x = element_text(color = "#70B7E6", size = 10, hjust = 1, vjust = 0.5, angle = 90))+
  geom_vline(xintercept = 6.5)
ggsave(paste0(path,"output/lgA_NC.pdf"), NC, width = 580/100, height = 687/100, dpi = 100)
ggsave(paste0(path,"output/lgA_TM.pdf"), TM, width = 580/100, height = 687/100, dpi = 100)



## Figure 6I ## 
# SCENIC heatmap for B
avg_celltype <- read.csv(file = "output/avg_celltype.csv", header = T, row.names = 1)
colnames(avg_celltype)[29]<-"delta Gamma T cells"
typelist<-table(cell_info[cell_info$type=="B",]$celltypes)
typelist<-rownames(typelist)[typelist>0]
reg<-sapply(typelist,function(i){
  gettopreg(i,10)
})
reg<-as.vector(reg)
reg<-reg[!duplicated(reg)]
colnames(avg_celltype)[1:4]<-typelist
data<-avg_celltype[reg,typelist]
scaled.data<-t(scale(t(data)))

B_order <- c("Naïve mature B cells", "Memory B cells", "Plasmablasts", "IgA+ plasma cells")
order <- B_order
scaled.data<-scaled.data[,order]
plot.data<-data.frame(scaled.data[,order],
                      "rank"=c(apply(scaled.data,1,function(data){match(max(data),data)})),
                      "index"=c(apply(scaled.data,1,function(data){match(max(data),data)})))
plot.data<-split(plot.data,plot.data$rank)
plot.data<-lapply(plot.data,function(x){x<-dplyr::arrange(x,-x[,unique(x$index)])})
plot.data<-do.call(rbind,plot.data)
plot.data<-plot.data[,-c(length(plot.data)-1,length(plot.data))]
rownames(plot.data)<-t(data.frame(strsplit(rownames(plot.data),"\\.")))[,2]
colnames(plot.data)<-order

library(ComplexHeatmap)
# downstream.Rmd
anncol<-data.frame(celltype = factor(colnames(plot.data),levels = colnames(plot.data)),
                   row.names = colnames(plot.data),
                   stringsAsFactors = F)
anncolors<-list()
anncolors[["celltype"]]<-c("IgA+ plasma cells" = "#8ED3C7",
                           "Naive mature B cells" = "#FFFFB3",
                           "Memory B cells" = "#B5B0D4",
                           "Plasmablasts" = "#FA7A6C")

plot.data[plot.data>2]=2
plot.data[plot.data<(-2)]=-2
plot.data <- plot.data[-c(13, 14, 17),]
pdf("output/B.pdf",width = 8, height = 10)
pheatmap(plot.data,
         border_color = "white",
         color = colorRampPalette(colors = c("#75A3D0","white","red"))(100),
         cluster_rows = F,
         cluster_cols = F,
         annotation_col = anncol[colnames(plot.data),,drop=F],
         annotation_colors = anncolors,
         cellwidth = 30, cellheight = 12,
         # gaps_row = c(2,8,13,19,22,26,32)
         # gaps_row = c(14,26,27)
         # gaps_row = c(7,12,16,22,30)
         gaps_row = c(8,12,17)
)
invisible(dev.off())




