## Figure S2A ## 
plot2<-DimPlot(seurat, reduction = "umap",
               group.by ='seurat_clusters',label = T)
plot2
ggsave("Epi cluster/UMAP 20211204 seurat_clusters.pdf", plot = plot2, width = 6, height = 5)

## Figure S2C ## 
Idents(seurat) <- "group" 
seurat<- RenameIdents(seurat, 'NC'='NC','HP'='HP','SSL'='SSL','SSLD'='SSLD',
                      'TSA'='TSA')
seurat$groups<-Idents(seurat)
plot2<-DimPlot(seurat, reduction = "umap",
               group.by ='groups',label = T)
plot2
ggsave("Epi cluster/UMAP 20220119 group.pdf", plot = plot2, width = 6, height = 5)


## Figure S2D ## 
plot1<-VlnPlot(seurat, features = c('ENO1'), pt.size=0, group.by="groups", ncol=3)+geom_boxplot(width=0.2,col="black",fill="#E1E3E6")
ggsave("备用/Figure S2/Vlnplot ENO1 .pdf", plot = plot1, width = 10, height = 4)
plot2<-VlnPlot(seurat, features = c('PKM'), pt.size=0, group.by="groups", ncol=3)+geom_boxplot(width=0.2,col="black",fill="#E1E3E6")
plot2
ggsave("备用/Figure S2/Vlnplot PKM .pdf", plot = plot2, width = 10, height = 4)

plot3<-VlnPlot(seurat, features = c('NQO1'), pt.size=0, group.by="groups", ncol=3)+geom_boxplot(width=0.2,col="black",fill="#E1E3E6")
ggsave("备用/Figure S2/Vlnplot NQO1 .pdf", plot = plot3, width = 10, height = 4)
plot4<-VlnPlot(seurat, features = c('LYZ'), pt.size=0, group.by="groups", ncol=3)+geom_boxplot(width=0.2,col="black",fill="#E1E3E6")
ggsave("备用/Figure S2/Vlnplot LYZ .pdf", plot = plot4, width = 10, height = 4)



## Figure S2E ## 
library(ggplot2)
expr <- openxlsx::read.xlsx("GSE76987_RightColonProcessed(1).xlsx")
expr$gene <- substr(paste0(expr$Ensembl_ID, expr$gene), 16, 30)
expr <- subset(expr, !gene %in% c("NA", ""))
expr <- expr[!duplicated(expr$gene), ]
rownames(expr) <- expr$gene
expr <- expr[, -c(1:3)]
plot.sample <- c(paste0("FPKM.CR-", seq(1:10)), 
                 paste0("FPKM.AP-", seq(1:10)), 
                 paste0("FPKM.SSA/P-", seq(1:21)))
plot.sample%in%colnames(expr)
expr <- expr[, plot.sample]

Group <- c(rep("NC", 10), rep("AP", 10), rep("SSL", 21))

cell.col <- c("NC" = "#DFC27D",
              "AP" = "#689E45",
              "SSL" = "#5A8CA8")

geneset <- msigdbr::msigdbr(species = "Homo sapiens")
geneset <- split(geneset$gene_symbol, geneset$gs_name)
names(geneset)[grep("OXIDATIVE_STRESS", names(geneset))]


res <- GSVA::gsva(expr = as.matrix(expr),
                  gset.idx.list = geneset,
                  kcdf = "Poisson",
                  min.sz = 15, max.sz = 500,
                  parallel.sz = 40)

pathway = rownames(res)[3]
pdf(paste0("GSE76987_", pathway, "_boxplot.pdf"), width = 350/100, height = 500/100)
plot.data <- data.frame(
  "Sample" = colnames(expr),
  "Group" = Group,
  "GSVA" = res[pathway, ]
)
plot.data$Group <- factor(plot.data$Group, levels = c("NC", "AP", "SSL"))
ggplot(data = plot.data,aes(x = Group, #分组列名
                            y = GSVA , #连续变量列名
                            fill = Group))+ #按分组填充颜色
  scale_fill_manual(values = cell.col) + #用自定义颜色填充
  scale_color_manual(values = cell.col) + #用自定义颜色填充
  geom_violin(alpha = 0.4, position = position_dodge(width = .75), 
              size = 0.8, color="white") +
  
  geom_point(shape = 21, size=1.2, # 点的性状和大小
             position = position_jitterdodge(), # 让点散开
             aes(color=Group), alpha = 0.6) +
  geom_boxplot(notch = TRUE, outlier.size = -1, 
               color="black", lwd = 0.6, alpha = 0.7) +
  theme_bw()+
  ggtitle(paste0(pathway, " GSVA in GSE76987")) +
  theme(axis.text.x = element_text(hjust = 1, angle = 45,color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "top",
        legend.title = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12)) +
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("AP", "NC"), c("NC", "SSL")))
invisible(dev.off())



## Figure S2F ## 
data.path <- file.path(work.path, "input");dir.create(data.path, showWarnings = F)
fig.path <- file.path(work.path, "figure");dir.create(data.path, showWarnings = F)
library(TCGAbiolinks)
library(ggplot2)
Counts <- read.table(file.path(data.path, "TCGA_CRC_Counts.txt"), sep = "\t", header = T, row.names = 1)
mut <- read.table(file.path(data.path, "TCGA_CRC_MAF.txt"), sep = "\t", header = T, row.names = 1)
mut <- t(mut["BRAF", ])
geneset <- msigdbr::msigdbr(species = "Homo sapiens")
geneset <- split(geneset$gene_symbol, geneset$gs_name)
names(geneset)[grep("OXIDATIVE_STRESS", names(geneset))]
geneset <- list(
  "RESPONSE_TO_OXIDATIVE_STRESS" = geneset[["GOBP_RESPONSE_TO_OXIDATIVE_STRESS"]]
)
geneset <- lapply(geneset, function(x){
  ensembldb::select(EnsDb.Hsapiens.v86, 
                    keys = x, 
                    keytype = "SYMBOL", columns = "GENEID")$GENEID
})

res <- GSVA::gsva(expr = as.matrix(Counts),
                  gset.idx.list = geneset,
                  kcdf = "Poisson",
                  min.sz = 15, max.sz = 500,
                  parallel.sz = 40)
all(colnames(Counts) == colnames(res))

cell.col <- RColorBrewer::brewer.pal(name = "Set1", n = 8)
cell.col <- c("WT" = cell.col[3], "M" = cell.col[2])


pdf("GSVAaboutBRAFmut0118.pdf", width = 330/100, height = 590/100, onefile = T)
lapply(rownames(res), function(pathway){
  plot.data <- data.frame(
    row.names = colnames(Counts),
    "GSVA" = res[pathway,],
    "BRAF_mut" = mut[match(paste0(substr(colnames(Counts), 1, 4), substr(colnames(Counts), 8, 16)), rownames(mut)), ]
  )
  plot.data$BRAF_mut[is.na(plot.data$BRAF_mut)] <- 0
  plot.data$BRAF_mut = ifelse(test = plot.data$BRAF_mut == 0, yes = "WT", no = "M")
  plot.data$BRAF_mut = factor(plot.data$BRAF_mut, levels = c("WT", "M"))
  
  ggplot(data = plot.data,aes(x = BRAF_mut, #分组列名
                              y = GSVA , #连续变量列名
                              fill = BRAF_mut))+ #按分组填充颜色
    scale_fill_manual(values = cell.col) + #用自定义颜色填充
    scale_color_manual(values = cell.col) + #用自定义颜色填充
    geom_violin(alpha = 0.4, position = position_dodge(width = .75), 
                size = 0.8, color="white") +
    
    geom_point(shape = 21, size=1.2, # 点的性状和大小
               position = position_jitterdodge(), # 让点散开
               aes(color=BRAF_mut), alpha = 0.6) +
    geom_boxplot(notch = TRUE, outlier.size = -1, 
                 color="black", lwd = 0.6, alpha = 0.7) +
    theme_bw()+
    xlab("BRAF_mut") +
    ggtitle(paste0(pathway, " GSVA")) +
    theme(axis.text.x = element_text(hjust = 1, size = 10, angle = 45,color = "black"),
          axis.text.y = element_text(size = 10, color = "black"),
          axis.ticks = element_line(size=0.2, color="black"),
          axis.ticks.length = unit(0.2, "cm"),
          legend.position = "top",
          legend.title = element_blank(),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 10)) +
    stat_compare_means(method = "wilcox.test", label.y = max(plot.data$GSVA)+0.2, label.x = 1.3)
})
invisible(dev.off())
