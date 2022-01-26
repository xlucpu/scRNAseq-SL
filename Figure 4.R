##### Figure plotting #####
## Figure 4A ## 
# set working path
workdir <- "~/data/ZYJ/9"; setwd(workdir)

# load package
library(ggplot2)

# load data
go <- read.table(file = "go.select2.txt",sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T)
go <- go[c(6:2,12:8),]
go$Term <- factor(go$Description, levels = go$Description)
#go$Term <- go$Description

go$GeneRatio <- as.numeric(sapply(strsplit(go$GeneRatio,"/",fixed = T),"[",1))/as.numeric(sapply(strsplit(go$GeneRatio,"/",fixed = T),"[",2))
ggplot(go, # you can replace the numbers to the row number of pathway of your interest
       aes(x = GeneRatio, y = Term)) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 10) +
  theme(axis.text = element_text(size = 10, colour = "black")) + 
  scale_colour_gradient(limits=c(0, 0.05), low="red") +
  ylab(NULL) + facet_grid(ONTOLOGY~.,scales = "free") + 
  ggtitle("GO pathway enrichment")
ggsave(filename = "go enrichment.pdf", width = 7,height = 4)



## Figure 4B ## 
library(ComplexHeatmap)
library(Seurat)
load("/raw/seurat_epi.Rdata")
Idents(seurat) <- seurat$celltypes
mtx <- AverageExpression(seurat)[[1]]
plotgenes <- c("GSTP1", "NQO1", "NQO2", "GPX2", "MGST1", "TXN", "PRDX5", "AQP5")
plotgenes%in%rownames(mtx)
plotcelltypes <- names(table(seurat$celltypes))[1:8]
plotcelltypes <- plotcelltypes[-c(4,6)]
plotcelltypes <- plotcelltypes[c(1,3,2,4:6)]
all(plotgenes%in%rownames(mtx))
plot.data <- mtx[plotgenes,plotcelltypes]
plot.data <- t(scale(t(plot.data)))
range(plot.data)
plot.data[plot.data>1]=1
plot.data[plot.data<(-1)]=-1

anncol<-data.frame(celltype = factor(colnames(plot.data),levels = colnames(plot.data)),
                   row.names = colnames(plot.data),
                   stringsAsFactors = F)
anncolors <- list()
anncolors[["celltype"]]<-c("Epi-2" = "#0A6BAD",
                           "Epi-3" = "#E10C0D",
                           "Epi-1" = "#33A02C",
                           "Transit amplifying cells" = "#C9B1D5",
                           "Goblet cells" = "#FB9998",
                           "Intermediate" = "#FDBE6E",
                           "BEST4+ colonocytes" = "#A6CEE3",
                           "Mature colonocytes" = "#B2DF8A") 
pdf("plot/Heatmap/Epi1113.pdf", width = 402/100, height = 327/100,)
pheatmap(plot.data,
         cluster_cols = F, cluster_rows = F,
         border_color = "white",
         # color = colorRampPalette(colors = c("#6353A3","white","#9A0038"))(100),
         color = colorRampPalette(colors = c("#75A3D0","white","red"))(100),
         annotation_col = anncol[colnames(plot.data),,drop=F],
         annotation_colors = anncolors,
         cellwidth = 12, cellheight = 12,
         treeheight_row = 20)
# gaps_row = 2)
invisible(dev.off())



## Figure 4D ## 
library(ComplexHeatmap)
auc.seu<-readRDS("auc.seu.rds") # AUCscore
rssMat <- readRDS("rssMat.rds") # RSSscore
avg_celltype <- read.csv(file = "output/avg_celltype.csv", header = T, row.names = 1) # Average AUCscore of each cluster
cell_info <- FetchData(auc.seu, vars = c("group","type","division","celltypes")) # cell info
gettopreg<-function(celltypes,k){
  data<-rssMat[celltypes]
  data<-dplyr::arrange(data,-data[,1])
  top_reg<-rownames(data)[1:k]
  return(top_reg)
}
typelist<-table(cell_info[cell_info$type=="epithelic",]$celltypes)
typelist<-rownames(typelist)[typelist>0]
reg<-sapply(typelist,function(i){
  gettopreg(i,10)
})
reg<-as.vector(reg)
reg<-reg[!duplicated(reg)]
colnames(avg_celltype)[11:20]<-typelist
data<-avg_celltype[reg,typelist]
scaled.data<-t(scale(t(data)))
epi_order<- c("Epi-2","Epi-3","Epi-1","Transit amplifying cells","Goblet cells","Intermediate",
              "BEST4+ colonocytes","Mature colonocytes","Tuft cells","Enteroendocrine cells")
order<-epi_order
scaled.data<-scaled.data[,order]
plot.data<-data.frame(scaled.data[,order],
                      "rank"=c(apply(scaled.data,1,function(data){match(max(data),data)})),
                      "index"=c(apply(scaled.data,1,function(data){match(max(data),data)})))
plot.data<-split(plot.data,plot.data$rank)
plot.data<-lapply(plot.data,function(x){x<-dplyr::arrange(x,-x[,unique(x$index)])})
plot.data<-do.call(rbind,plot.data)
plot.data<-do.call(rbind,plot.data)
plot.data<-plot.data[,-c(length(plot.data)-1,length(plot.data))]
rownames(plot.data)<-t(data.frame(strsplit(rownames(plot.data),"\\.")))[,2]
colnames(plot.data)<-order

anncol<-data.frame(celltype = factor(colnames(plot.data),levels = colnames(plot.data)),
                   row.names = colnames(plot.data),
                   stringsAsFactors = F)
anncolors<-list()
#epithelial
anncolors[["celltype"]]<-c("Epi-2" = "#0A6BAD",
                           "Epi-3" = "#E10C0D",
                           "Epi-1" = "#33A02C",
                           "Transit amplifying cells" = "#C9B1D5",
                           "Goblet cells" = "#FB9998",
                           "Intermediate" = "#FDBE6E",
                           "BEST4+ colonocytes" = "#A6CEE3",
                           "Mature colonocytes" = "#B2DF8A",
                           "Tuft cells" = "#FF7F00",
                           "Enteroendocrine cells" = "#6A3D9A")
plot.data[plot.data>2]=2
plot.data[plot.data<(-2)]=-2
plot.data <- plot.data[c(1:2, 6, 31, 29, 30, 32, 3:5, 7:11, 12:28, 33:40),]
pdf("output/Epi1027.pdf",width = 8, height = 10)
pheatmap(plot.data,
         border_color = "white",
         color = colorRampPalette(colors = c("#75A3D0","white","red"))(100),
         cluster_rows = F,
         cluster_cols = F,
         annotation_col = anncol[colnames(plot.data),,drop=F],
         annotation_colors = anncolors,
         cellwidth = 25, cellheight = 12,
         gaps_row = c(7,11,17,23,26,30,32)
         # gaps_row = c(14,26,27)
         # gaps_row = c(7,12,16,22,30)
         # gaps_row = c(8,12,17)
)
invisible(dev.off())



## Figure 4E ## 
Idents(seurat) <- "group" 
Epicell<-subset(seurat,idents=c('NC','SSL','SSLD','TSA'))
Epicell@meta.data$groups<-Epicell@meta.data$group
Idents(Epicell) <- "groups" 
Epicell<- RenameIdents(Epicell, 'NC'='NC','SSL'='SSL','SSLD'='SSLD',
                       'TSA'='TSA')
Epicell$groups<-Idents(Epicell)
plot4<-VlnPlot(Epicell, features =c('MKI67','FOXQ1','NQO2','AQP5'), pt.size=0.5, group.by="groups", ncol=3)
ggsave("Epi cell_identify/氧化应激 TA heatmap/VlnplotB 20211113.pdf", plot = plot4, width = 8, height = 6)



## Figure 4F ## 
path <- "~/data/ZYJ/11";setwd(path)

heatmap.BlWtRd <- c("#6699CC","white","#FF3C38")
heatmap.YlGnPe <- c("#440259","#345F8C","#228C8A","#78CE51","#FAE71F")
heatmap.GrWtRd <- c("#2b2d42","#8d99ae","#edf2f4","#ef233c","#d90429")
heatmap.L.BlYlRd <- c("#4281a4","#9cafb7","#ead2ac","#e6b89c","#fe938c")
heatmap.BlBkRd <- c("#54FEFF","#32ABAA","#125456","#000000","#510000","#A20000","#F30000")
heatmap.fancy <- c("#10040A", "#2A0B35", "#4D155B", "#73215B", "#9C3558", "#C34D44", "#E07038", "#F2981C", "#F2CA51", "#FAF6A3")



expr <- openxlsx::read.xlsx("GSE76987_RightColonProcessed(1).xlsx")
plot.gene <- c("GSTP1", "NQO1", "NQO2", "DUOX1", "DUOX2", "DUOXA1", "DUOXA2", "AQP5")
plot.sample <- c(paste0("FPKM.CR-", seq(1:10)), 
                 paste0("FPKM.AP-", seq(1:10)), 
                 paste0("FPKM.SSA/P-", seq(1:21)))
plot.sample%in%colnames(expr)
plotdata <- expr[match(plot.gene, expr$gene), plot.sample]
rownames(plotdata) <- plot.gene
plotdata <- as.matrix(plotdata)
plotdata <- log(plotdata+1)
range(plotdata)
plotdata <- t(scale(t(plotdata)))
range(plotdata)
plotdata[plotdata<(-2)] <- -2
plotdata[plotdata>2] <- 2
range(plotdata)
plot(density(plotdata))
Group <- c(rep("NC", 10), rep("AP", 10), rep("SSL", 21))
anncol <- data.frame("Group" = factor(Group, levels = unique(Group)),
                     row.names = colnames(plotdata))
anncolors <- list()
anncolors[["Group"]] <- c("NC" = "#DFC27D",
                          "AP" = "#689E45",
                          "SSL" = "#5A8CA8")
# plot(x = as.factor(Group), y = plotdata[4,])

library(ComplexHeatmap)
pdf("GSE76987_FPKM_Heatmap.pdf", width = 1011/100, height = 443/100)
pheatmap(mat = plotdata, 
         cluster_rows = F, cluster_cols = F,
         border_color = "NA", gaps_col = c(10,20),
         # color = colorRampPalette(colors = c("#75A3D0","white","red"))(100),
         color = colorRampPalette(colors = heatmap.fancy)(100),
         cellwidth = 10, cellheight = 20,
         angle_col = "45",
         fontsize_col = 12, fontsize_row = 10,
         annotation_col = anncol, annotation_colors = anncolors,
         fontsize = 15,
         show_colnames = F,
         name = "Scaled log(FPKM+1)")
invisible(dev.off())