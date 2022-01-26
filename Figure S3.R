## Figure S3A ## 
plot5<-VlnPlot(seurat, features = c('SERPINB6'), pt.size=0, group.by="groups", ncol=3)+geom_boxplot(width=0.2,col="black",fill="#E1E3E6")
ggsave("备用/Figure S3 SERPINB6/Vlnplot B6.pdf", plot = plot5, width = 10, height = 4)


## Figure S3B ##
plot1<-VlnPlot(seurat, features = c('SERPINB6'), pt.size=0, group.by='celltypes2', ncol=3)+geom_boxplot(width=0.2,col="black",fill="#E1E3E6", outlier.size = 0.5)
plot1
ggsave("备用/Figure S3 SERPINB6/Vlnplot B6 EpiSL.pdf", plot = plot1, width = 15, height = 5)


## Figure S3G ##
## transcriptome profile
library(TCGAbiolinks)
library(SummarizedExperiment)
TCGAbiolinks::getGDCprojects()

project = c("TCGA-COAD", "TCGA-READ")
Expr_info <- GDCquery(project = project,
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "HTSeq - FPKM"
)
GDCdownload(query = Expr_info, directory = work.path)
expr <- GDCprepare(query = Expr_info, directory = work.path)
FPKM <- assay(expr)
clinical <- GDCquery_clinic(project = project[2],type = "clinical")
write.table(x = FPKM, file = file.path(data.path, "TCGA_CRC_FPKM.txt"), sep = "\t", row.names = T, col.names = T)


## mutation
library(TCGAmutations)
library(reshape2)
TCGAmutations::tcga_available()
maf <- TCGAmutations::tcga_load(study = c("COAD", "READ"))
if (class(maf) == "list"){
  maf.tab <- lapply(maf, function(x){x@data})
  maf.tab <- do.call(rbind, maf.tab)
}else{
  maf.tab <- as.data.frame(maf@data)
}

mut <- as.data.frame.matrix(table(maf.tab$Hugo_Symbol, maf.tab$Tumor_Sample_Barcode))
length(unique(mut))==length(unique(paste0(substr(colnames(mut),1,4),substr(colnames(mut),8,16))))
colnames(mut) <- paste0(substr(colnames(mut),1,4),substr(colnames(mut),8,16))
write.table(x = mut, file = file.path(data.path, "TCGA_CRC_MAF.txt"), sep = "\t", row.names = T, col.names = T)
mut <- t(mut["BRAF", ])


# Plot --------------------------------------------------------------------
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(ggpubr)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(FPKM), keytype = "GENEID", columns = c("SYMBOL","GENEID"))
rownames(FPKM) <- geneIDs1$SYMBOL

patient <- paste0(substr(colnames(FPKM),1,4),substr(colnames(FPKM),8,16))
sum(is.na(match(patient, colnames(mut))))
patient.mut <- mut[match(patient, rownames(mut)),]
plot.gene <- c("SERPINB6", "SDR16C5")
plot.data <- as.data.frame(t(FPKM[plot.gene,]))
plot.data <- log10(plot.data+1)
plot.data$BRAF <- patient.mut
# plot.data <- subset(plot.data, !is.na(BRAF))
plot.data$BRAF[is.na(plot.data$BRAF)] <- 0
plot.data$BRAF_mut <- ifelse(test = plot.data$BRAF == 0, yes = "WT", no = "M")


cell.col <- RColorBrewer::brewer.pal(name = "Set1", n = 8)
cell.col <- c("WT" = cell.col[2], "M" = cell.col[3])

# plot1
ggplot(data = plot.data,aes(x = BRAF_mut, #分组列名
                            y = SERPINB6 , #连续变量列名
                            fill = BRAF_mut))+ #按分组填充颜色
  scale_fill_manual(values = cell.col) + #用自定义颜色填充
  geom_violin(alpha = 0.4, position = position_dodge(width = .75), 
              size = 0.8, color="white") +
  
  geom_point(shape = 21, size=1.2, # 点的性状和大小
             position = position_jitterdodge(), # 让点散开
             aes(color=BRAF_mut), alpha = 0.6) +
  scale_color_manual(values = cell.col) + #用自定义颜色填充
  geom_boxplot(notch = TRUE, outlier.size = -1, 
               color="black", lwd = 0.6, alpha = 0.7) +
  theme_bw()+
  # ylab("SERPINB6 expression") +
  xlab("BRAF_mut") +
  ggtitle("SERPINB6 expression") +
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
  stat_compare_means(method = "wilcox.test", label.y = 2.3, label.x = 1.5)
ggsave(filename = file.path(work.path, "TCGA_CRC_SERPINB6_withBRAFmutation1223.pdf"), width = 311/100, height = 511/100)
