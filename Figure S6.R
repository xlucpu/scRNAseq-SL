##### Figure plotting #####
## Figure S6A ## 
setwd('D:/2021scRNA-seq')
load("seurat_CD8T加naive.Rdat")
#####CD8+T细胞再分群
CD8T<-subset(seurat,idents='CD8+ T cells')
#新CD8T 包括naive C11, C14
Idents(seurat) <- "celltypes"
CD8T<-subset(seurat,idents=c('CD8+ T cells','Naïve CD8+ T cells'))
CD8T <- FindVariableFeatures(CD8T, selection.method = "vst", nfeatures = 4000)
top10<-head(VariableFeatures(CD8T),10)
top10
scale.genes <-  rownames(CD8T)
CD8T <- ScaleData(CD8T, features = scale.genes)
CD8T<- RunPCA(CD8T, features = VariableFeatures(object = CD8T))
plot2 <- ElbowPlot(CD8T, ndims=20, reduction="pca") 
plot2
pc.num=1:11
CD8T<- FindNeighbors(CD8T, dims = pc.num)
CD8T <- FindClusters(CD8T, resolution = 0.5)
table(CD8T@meta.data$seurat_clusters)
metadata <- CD8T@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)

CD8T <- RunUMAP(CD8T, n.neighbors = 10,n.epochs=10,dims = pc.num)
embed_umap <- Embeddings(CD8T, 'umap')
write.csv(embed_umap,'T cluster/CD8T cluster/embed_umap n.neighbors = 50,n.epochs=50 20210427.csv') 
plot2 = DimPlot(CD8T,label = T, reduction = "umap") 
plot2
ggsave("T cluster/CD8T cluster/UMAP 20210427加naive.pdf", plot = plot2, width = 8, height = 7)
diff.wilcox = FindAllMarkers(CD8T)
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write.csv(all.markers, "T cell_identify/CD8T cell_identify/diff_genes_wilcox加naive.csv", row.names = F)
write.csv(top10, "T cell_identify/CD8T cell_identify/top10_diff_genes_wilcox加naive.csv", row.names = F)

Idents(CD8Teff) <- "group" 
diff_dat <- FindMarkers(CD8Teff,ident.1=c('SSL','SSLD'),ident.2 = 'NC',
                        group.by='group')#SSL+SSLD T的差异基因
diff_dat2 <- diff_dat[diff_dat$p_val<0.05,]
head(diff_dat2)
write.csv(diff_dat2, "T cell_identify/CD8T cell_identify/CD8Teff SSL+SSLD T的差异基因.csv", row.names = T)
plot4<-DimPlot(CD8T, reduction = "umap",group.by = "seurat_clusters",label = T)
ggsave("T cluster/CD8T cluster/UMAP 20211216.pdf", plot = plot4, width = 6, height = 5)



## Figure S6B ## 
CD8T@meta.data$groups<-CD8T@meta.data$group
Idents(CD8T) <- "groups" 
CD8T<- RenameIdents(CD8T, 'NC'='NC','HP'='HP','SSL'='SSL','SSLD'='SSLD',
                    'TSA'='TSA')
CD8T$groups<-Idents(CD8T)
plot2<-VlnPlot(CD8T, features =c('ITGAE'), pt.size=0, group.by="groups", ncol=3)
ggsave("T cell_identify\\CD8T cell_identify\\20210823 TRM\\CD103 groups.pdf", plot = plot2, width = 8, height = 3.2)



## Figure S6F ## 
work.path <- "i:/genomicdata/external/ZYJ/16";setwd(work.path)
data.path <- file.path(work.path, "input");dir.create(data.path, showWarnings = F)
fig.path <- file.path(work.path, "figure");dir.create(data.path, showWarnings = F)

library(EnsDb.Hsapiens.v86)
library(SummarizedExperiment)

fpkmToTpm <- function(fpkm){
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

gse39582.expr <- read.table("gse39582.expr.txt", sep = "\t", row.names = 1, header = T)

projects = c("TCGA-COAD", "TCGA-READ")
Expr_info <- GDCquery(project = projects,
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "HTSeq - FPKM"
)
expr <- GDCprepare(query = Expr_info, directory = "i:/genomicdata/external/ZYJ/15")
mut <- read.table("i:/genomicdata/external/ZYJ/15/input/TCGA_CRC_MAF.txt", sep = "\t", header = T, row.names = 1)
mut <- t(mut["BRAF", ])

tcga.FPKM <- assay(expr)
tcga.expr <- log(fpkmToTpm(tcga.FPKM)+1)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, 
                              keys= rownames(tcga.expr),
                              keytype = "GENEID", 
                              columns = c("SYMBOL","GENEID"))
rownames(tcga.expr) <- geneIDs1$SYMBOL

length(intersect(rownames(gse39582.expr), rownames(tcga.expr)))
gene <- intersect(rownames(gse39582.expr), rownames(tcga.expr))

b <- readRDS("i:/genomicdata/jiaoda/gsva/TRMgeneset.rds")
geneset <- list(
  "TRM" = b$TRM
)

merge.expr <- sva::ComBat(dat = cbind(tcga.expr[gene, ], gse39582.expr[gene, ]),
                          batch = c(rep("tcga", ncol(tcga.expr)), rep("gse", ncol(gse39582.expr))),
                          par.prior=TRUE, prior.plots=TRUE)

res <- GSVA::gsva(expr = as.matrix(merge.expr),
                  gset.idx.list = geneset,
                  kcdf = "Gaussian",
                  min.sz = 15, max.sz = 500,
                  parallel.sz = 40)

info <- read.table("GSE39582 MSI(1).csv", sep = ",", header = T)
info <- subset(info, new != "")
all(info$GEOID %in% colnames(res))
info$GSVA <- res[, match(info$GEOID, colnames(res))]

msi_results <- GDCprepare_clinic(GDCquery(project = projects, data.category = "Other", legacy = TRUE,
                                          access = "open", data.type = "Auxiliary test"),
                                 clinical.info = "msi", 
                                 directory = "i:/genomicdata/external/ZYJ/15/GDCdata")
data <- data.frame(
  "Sample" = msi_results$bcr_patient_barcode,
  "MSI" = msi_results$mononucleotide_and_dinucleotide_marker_panel_analysis_status,
  "BRAF" = mut[match(gsub("(TCGA-)(\\w+\\-)(\\w+)", "\\1\\3", msi_results$bcr_patient_barcode), 
                     substr(gsub("\\.", "-", rownames(mut)), 1, 9)),],
  "GSVA" = res[, match(msi_results$bcr_patient_barcode, substr(colnames(res), 1, 12))]
)
data$BRAF[is.na(data$BRAF)] = 0
data$BRAF <- ifelse(data$BRAF>0, yes = "BRAF", no = "WT")
data <- subset(data, MSI %in% c("MSI-L", "MSI-H", "MSS"))
data$MSI <- ifelse(test = data$MSI == "MSI-H", yes = "MSI", no = "MSS")
data$Group <- paste0(data$BRAF, " ", data$MSI)

plot.data <- data.frame(
  "Sample" = c(info$GEOID, data$Sample),
  "Group" = c(info$new, data$Group),
  "GSVA" = c(info$GSVA, data$GSVA)
)


pathway = "TRM"

cell.col <- c("BRAF MSI" = "#FF928B", "BRAF MSS" = "#99C231",
              "WT MSI" = "#22C7CC", "WT MSS" = "#D8A4FF")

plot.data$Group <- factor(plot.data$Group, levels = rev(c("WT MSS", "BRAF MSS", "WT MSI", "BRAF MSI")))

p.val <- kruskal.test(GSVA ~ Group,
                      data = plot.data)
p.lab <- paste0("P",
                ifelse(p.val$p.value < 0.001, " < 0.001",
                       paste0(" = ",round(p.val$p.value, 3)))) 
p.lab

p_top <- ggplot(plot.data, aes(x = GSVA, color = Group, fill = Group)) +
  geom_density() +
  # 让箱子的所有位置都颜色统一，如例文所示
  scale_color_manual(values = alpha(cell.col, 0.7)) + # 设置透明色
  scale_fill_manual(values = alpha(cell.col, 0.7)) +
  theme_classic() + # 如果显示采用这一行
  
  # 这里提取输入文件的第一列药物名称，写入x轴标签
  xlab(paste0("Estimated GSVA of ", pathway)) + 
  # 第一列非必需，可以像下面这样直接写xlab
  #xlab("Estimated AUC of Cisplatin") +
  ylim(-0.14, 2.5) +
  ylab(NULL) + 
  theme(legend.position = "none", 
        legend.title = element_blank(),
        axis.text.x = element_text(size = 12,color = "black"),
        axis.text.y = element_blank(), # 原文不显示纵轴的密度
        #axis.text.y = element_text(size = 12,color = "black"), # 如果要显示采用这一行
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_rug(length = unit(0.08, "npc"))

p_bot <- ggplot(plot.data, aes(x = Group, y = GSVA, fill = Group)) + 
  geom_boxplot(aes(col = Group), notch = TRUE, outlier.size = -1, 
               color="black", lwd = 0.3, alpha = 0.7) + 
  scale_fill_manual(values = cell.col) + 
  scale_color_manual(values = cell.col) + 
  xlab(NULL) + ylab(paste0(pathway, "GSVA score")) + 
  theme_void() +
  theme(legend.position = "right",
        legend.title = element_blank(),
        axis.text.x = element_blank(), # 原文不显示箱线图的x轴
        #axis.text.x = element_text(size = 12,color = "black"), # 如要显示箱线图x轴采用这一行
        axis.text.y = element_text(size = 11,color = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  annotate(geom="text",
           x = 1.5,
           hjust = 1,
           y = 0.54,#max(plot.data$GSVA)
           size = 4, angle = 270, fontface = "bold",
           label = p.lab) +
  coord_flip()

p <- p_top %>% aplot::insert_bottom(p_bot, height = 0.4)
p
ggsave(paste0("BoxDensityPlot", pathway, ".pdf"), width = 650/100, height = 400/100)
ggsave(paste0("BoxDensityPlot", pathway, ".pdf"), width = 650/100, height = 400/100)
