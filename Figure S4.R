## Figure S4A ## 
####细胞周期
seurat<- CellCycleScoring(seurat,
                          g2m.features = cc.genes.updated.2019$g2m.genes,
                          s.features = cc.genes.updated.2019$s.genes)
VlnPlot(seurat,features = c("nCount_RNA","nFeature_RNA"),group.by = "Phase",pt.size = 0)
plot2<-DimPlot(seurat,group.by = "Phase",label = TRUE,repel = T)  
ggsave("Epi cluster/EpiSL cluster/UMAP CellCycle20211201.pdf", plot = plot2, width = 6, height = 5)



## Figure S4B-C ## 
path <- "/share/genomics/cwx/jiaoda/inferCNV";setwd(path)
rawData.path <- "/share/genomics/cwx/jiaoda/raw"


library(Seurat)
library(infercnv)
library(ggplot2)
library(ggpubr)
library(circlize)

load(file.path(rawData.path, "seurat_epi.Rdata"))
load(file.path(rawData.path, "Epi_metadata20211125V2.Rdat"))
# expr <- readRDS(file.path(rawData.path, "expr0617.rds"))
# colnames(expr) <- gsub("\\.", "-", colnames(expr))
expr <- seurat@assays$RNA@data
seurat@meta.data <- Epi_metadata
seurat$group <- plyr::mapvalues(x = seurat$group,
                                from = c("NC", "HP", "SSL", "SSLD", "TSA"),
                                to = c("NC", "TM", "TM", "TM", "TM"))
seurat$ident <- paste0(seurat$group, "_", seurat$celltypes4)
seurat <- subset(seurat, celltypes4 %in% c("Epi-1", "Epi-2", "Mature colonocytes", "BEST4+ colonocytes"))
cell_info <- FetchData(seurat, vars = c("group", "celltypes4", "ident"))
all(colnames(seurat)%in%colnames(expr))
expr <- expr[, colnames(seurat)]

GeneToENSG <- read.table("/share/genomics/cwx/jiaoda/raw/genes.tsv",
                         header = T, sep = "\t")
rownames(expr) <- GeneToENSG$gene_ids[match(x = rownames(expr), table = GeneToENSG$gene)]

gene_order_file <- read.table(file.path(path,"gencode.v38.txt"), 
                              sep = "\t", header = F,
                              row.names = 1,
                              check.names = F)
gene38<-strsplit(rownames(gene_order_file),"\\.")
gene38<-do.call(rbind,gene38)[,1]
gene_order_file <- gene_order_file[!duplicated(gene38),]
rownames(gene_order_file) <-gene38[!duplicated(gene38)]


infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = expr,
                                     annotations_file = cell_info,
                                     gene_order_file = gene_order_file,
                                     ref_group_names = names(table(cell_info))[substr(names(table(cell_info)),0,2)=="NC"])

infercnv_obj <- infercnv::run(infercnv_obj = infercnv_obj,
                              cutoff = 0.1,
                              out_dir = paste0(path,"epi_by_celltype_nosub1"),
                              cluster_by_groups = T,
                              denoise = T,
                              # noise_filter = 0.2,
                              # sd_amplifier = 2,
                              HMM = T,
                              # analysis_mode = "subclusters",
                              analysis_mode = "sample",
                              num_threads = 20)

infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = expr,
                                     annotations_file = cell_info,
                                     gene_order_file = gene_order_file,
                                     ref_group_names = names(table(cell_info))[substr(names(table(cell_info)),0,2)=="NC"])

infercnv_obj <- infercnv::run(infercnv_obj = infercnv_obj,
                              cutoff = 0.1,
                              out_dir = paste0(path,"epi_by_celltype_sub"),
                              cluster_by_groups = T,
                              denoise = T,
                              # noise_filter = 0.2,
                              # sd_amplifier = 2,
                              HMM = T,
                              analysis_mode = "subclusters",
                              # analysis_mode = "sample",
                              num_threads = 20)


## boxplot
cnv <- readRDS(file.path(path,"epi_by_celltype_sub/run.final.infercnv_obj"))
cnv_state <- readRDS(file.path(path,"epi_by_celltype_sub/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.infercnv_obj"))
all(rownames(cell_info)==colnames(cnv@expr.data))
cell_info <- subset(cell_info, substr(ident, 1, 2)!="NC")
data <- cnv@expr.data[, rownames(cell_info)]-1
data <- apply(data,2,function(x){sum(x^2)})
sum(names(data)==rownames(cell_info))
plotdata <- data.frame("CNV_score" = data,
                       "celltype" = cell_info$celltypes4,
                       "group" = as.factor(as.character(cell_info$group)),
                       "CNV_signal" = data/dim(cnv@expr.data)[1])
compare_type <- list(c("Epi-2","Mature colonocytes"))
compare_group <- list(c("NC","SSL"))
# pdf(paste0(path,"pan_epi_by_celltype/boxplot.pdf"))
ggplot(data = plotdata, mapping = aes(x = celltype, y = CNV_score, fill = celltype))+
  geom_boxplot(width=0.4)+
  theme(legend.position = "top")+
  theme_classic()+
  scale_color_manual(values = c("Epi-2" = "#E31A1C",
                                "Epi-1" = "#33A02C",
                                "BEST4+ colonocytes" = "#A6CEE3",
                                "Mature colonocytes" = "#B2DF8A"),
                     aesthetics = "fill")

cell.col = c("Epi-2" = "#E31A1C",
             "Epi-1" = "#33A02C",
             "BEST4+ colonocytes" = "#A6CEE3",
             "Mature colonocytes" = "#B2DF8A")

pdf("CNV_boxplot1207.pdf", width = 359/100, height = 514/100)
ggplot(data = plotdata, aes(x = celltype, 
                            y = log2(CNV_score + 1) , 
                            fill = celltype))+ 
  scale_fill_manual(values = cell.col) + 
  geom_violin(alpha = 0.4, position = position_dodge(width = .75), 
              size = 0.8, color="white") +
  
  geom_point(shape = 21, size=1.2, 
             position = position_jitterdodge(), 
             aes(color=celltype), alpha = 0.6) +
  scale_color_manual(values = cell.col) + 
  geom_boxplot(notch = TRUE, outlier.size = -1, 
               color="black", lwd = 0.6, alpha = 0.7) +
  theme_bw() + 
  ylab("log2(CNV score + 1)") +
  xlab("") + ylim(0,5.5) +
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
  stat_compare_means(method = "kruskal.test", label.y = 5.3, label.x = 1.5)
invisible(dev.off())


## heatmap
# cnv <- readRDS(file.path(path,"epi_by_celltype_sub/run.final.infercnv_obj"))
# ref_data <- read.csv(file.path(path,"epi_by_celltype_sub/infercnv.references.txt"), sep = " ")
# obs_data <- read.csv(file.path(path,"epi_by_celltype_sub/infercnv.observations.txt"), sep = " ")
cnv <- readRDS(file.path(path,"epi_by_celltype_sub/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.infercnv_obj"))
ref_data <- read.csv(file.path(path,"epi_by_celltype_sub/infercnv.17_HMM_predHMMi6.leiden.hmm_mode-subclusters.references.txt"), sep = " ")
obs_data <- read.csv(file.path(path,"epi_by_celltype_sub/infercnv.17_HMM_predHMMi6.leiden.hmm_mode-subclusters.observations.txt"), sep = " ")

# rownames(cell_info) <- gsub(pattern = "-", replacement = "\\.", x = rownames(cell_info))
ref_type <- data.frame("celltypes" = cell_info$ident[match(colnames(ref_data), rownames(cell_info))],
                       row.names = colnames(ref_data))
obs_type <- data.frame("celltypes" = cell_info$ident[match(colnames(obs_data), rownames(cell_info))],
                       row.names = colnames(obs_data))
gene <- data.frame("chr" = gene_order_file$V2[match(rownames(obs_data),rownames(gene_order_file))],
                   row.names = rownames(obs_data))

gene_annotation <- cnv@gene_order$chr
names(gene_annotation) <- rownames(cnv@gene_order)
genecol <- colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(24)
names(genecol) <- names(table(gene_annotation))
anncol <- data.frame("chr" = factor(gene$chr, levels = unique(gene$chr)),
                     row.names = rownames(gene))

obs_annrow <- data.frame("celltypes" = factor(obs_type$celltypes, levels = unique(obs_type$celltypes)),
                         row.names = rownames(obs_type))
ref_annrow <- data.frame("celltypes" = factor(ref_type$celltypes, levels = unique(ref_type$celltypes)),
                         row.names = rownames(ref_type))

anncolors <- list()
anncolors[["celltypes"]] <- c("NC_Epi-2" = "#0A6BAD", "TM_Epi-2" = "#0A6BAD",
                              "NC_Epi-1" = "#33A02C", "TM_Epi-1" = "#33A02C",
                              "NC_Epi-3" = "#E10C0D", "TM_Epi-3" = "#E10C0D",
                              "Transit amplifying cells" = "#C9B1D5",
                              "Goblet cells" = "#FB9998",
                              "Intermediate" = "#FDBE6E",
                              "NC_BEST4+ colonocytes" = "#A6CEE3",
                              "NC_Mature colonocytes" = "#B2DF8A", "TM_Mature colonocytes" = "#B2DF8A",
                              "Tuft cells" = "#FF7F00",
                              "Enteroendocrine cells" = "#6A3D9A")

pheatmap(t(ref_data),
         border_color = NA,
         color = colorRampPalette(c("darkblue","white","darkred"))(100),
         clustering_method = "ward.D",
         cluster_rows = T, cluster_cols = F,
         show_rownames = F,show_colnames = F,
         annotation_row  = ref_annrow,
         annotation_col = anncol,
         annotation_colors = anncolors,
)


## a heatmap using Heatmap function
#############################
cnv1 <- readRDS(file.path(path,"epi_by_celltype_sub/run.final.infercnv_obj"))
cnv <- readRDS(file.path(path,"epi_by_celltype_sub/20_HMM_pred.repr_intensitiesHMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.infercnv_obj"))
# cnv <- readRDS(file.path(path,"epi_by_celltype_sub/19_HMM_pred.Bayes_NetHMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.infercnv_obj"))
# cnv <- readRDS(file.path(path,"epi_by_celltype_sub/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.infercnv_obj"))

gene_annotation <- cnv@gene_order$chr
names(gene_annotation) <- rownames(cnv@gene_order)
genecol <- colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(22)
names(genecol) <- names(table(gene_annotation))
anncol <- HeatmapAnnotation(gene = gene_annotation,
                            col = list(gene = genecol))
ha <- HeatmapAnnotation(gene = anno_block(gp = gpar(fill = colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(22)),
                                          labels = paste0("chr",seq(1:22))))

ref <- c("NC_Epi-1", "NC_Epi-2", "NC_Mature colonocytes", "NC_BEST4+ colonocytes")
ref_plotdata <- lapply(ref, function(x){
  expr = cnv@expr.data[,cnv@reference_grouped_cell_indices[[x]]]
  # x = expr[,order.dendrogram(as.dendrogram(cnv@tumor_subclusters$hc[[x]]))]
  hcl <- hclust(d = dist(t(expr)), method = "ward.D")
  # a <- order.dendrogram(as.dendrogram(hcl))
  expr = expr[,order.dendrogram(as.dendrogram(hcl))]
})
ref_plotdata <- do.call(cbind,ref_plotdata)
# ref_annotation <- epi_info1[colnames(ref_plotdata),];names(ref_annotation) = colnames(ref_annotation)
ref_annotation <- factor(cell_info$ident[match(colnames(ref_plotdata), rownames(cell_info))], 
                         levels = unique(cell_info$ident[match(colnames(ref_plotdata), rownames(cell_info))]))
ref_annrow <- rowAnnotation(celltype = ref_annotation,
                            col = list(
                              celltype = c(
                                "NC_Epi-1" = "#33A02C",
                                "NC_Epi-2" = "#E31A1C",
                                "NC_Mature colonocytes" = "#B2DF8A",
                                "NC_BEST4+ colonocytes" = "#A6CEE3"
                              )
                            ))
tp <- Heatmap(t(ref_plotdata),
              width = ncol(t(ref_plotdata))*unit(0.04, "mm"),
              height = nrow(t(ref_plotdata))*unit(0.03, "mm"),
              # col = colorRamp2(c(0.9, 1 ,1.1),c("darkblue","white","darkred")),
              col = colorRamp2(c(0, 1 ,3),c("darkblue","white","darkred")),
              row_split = ref_annotation, row_gap = unit(0, "mm"), row_title = "NC",
              column_split = gene_annotation, column_gap = unit(0, "mm"), column_title = NULL,
              border = T, 
              cluster_rows = F, cluster_columns = F,
              show_row_names = F, show_column_names = F,
              left_annotation = ref_annrow
              # bottom_annotation =  ha
)

obs <- c("TM_Epi-1", "TM_Epi-2", "TM_Mature colonocytes", "TM_BEST4+ colonocytes")
obs_plotdata <- lapply(obs, function(x){
  expr = cnv@expr.data[,cnv@observation_grouped_cell_indices[[x]]]
  # x = expr[,order.dendrogram(as.dendrogram(cnv@tumor_subclusters$hc[[x]]))]
  hcl <- hclust(d = dist(t(expr)), method = "ward.D")
  # d <- order.dendrogram(as.dendrogram(hcl))
  expr = expr[, order.dendrogram(as.dendrogram(hcl))]
})
obs_plotdata <- do.call(cbind,obs_plotdata)
# obs_annotation <- epi_info1[colnames(obs_plotdata),];names(obs_annotation) = colnames(obs_annotation)
obs_annotation <- factor(cell_info$ident[match(colnames(obs_plotdata), rownames(cell_info))], 
                         levels = unique(cell_info$ident[match(colnames(obs_plotdata), rownames(cell_info))]))
obs_annrow <- rowAnnotation(celltype = obs_annotation,
                            col = list(
                              celltype = c(
                                "TM_Epi-1" = "#33A02C",
                                "TM_Epi-2" = "#E31A1C",
                                "TM_Mature colonocytes" = "#B2DF8A",
                                "TM_BEST4+ colonocytes" = "#A6CEE3"
                              )
                            ))
bp <- Heatmap(t(obs_plotdata),
              width = ncol(t(obs_plotdata))*unit(0.04, "mm"),
              height = nrow(t(obs_plotdata))*unit(0.03, "mm"),
              # col = colorRamp2(c(0.9, 1 ,1.1), c("darkblue","white","darkred")),
              col = colorRamp2(c(0, 1, 3), c("darkblue","white","darkred")),
              row_split = obs_annotation, row_gap = unit(0, "mm"), row_title = "TM",
              column_split = gene_annotation, column_gap = unit(0, "mm"), column_title = NULL,
              border = T, 
              cluster_rows = F, cluster_columns = F,
              show_row_names = F, show_column_names = F,
              left_annotation = obs_annrow,
              top_annotation =  ha
)
pdf(file.path(path,"inferCNV_predicted_state_subcluster_heatmap1208.pdf"), width = 38*0.39, height = 26*0.39)
draw(tp%v%bp, column_title = "state predicted by inferCNV")
invisible(dev.off())

plot(density(as.vector(cnv@expr.data)))
plot(density(ref_plotdata))
plot(density(obs_plotdata))

plot.data <- read.table(file.path(path,"epi_by_celltype_sub/infercnv.20_HMM_predHMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt"))
p1 <- pheatmap(t(plot.data), 
               color = colorRamp2(c(0, 1, 3), c("darkblue","white","darkred")), 
               cluster_cols = F, cluster_rows = T,
               show_rownames = F, show_colnames = F)


cl_cb <- function(hcl, mat){
  # Recalculate manhattan distances for reorder method
  dists <- dist(mat)
  
  # Perform reordering according to OLO method
  hclust_olo <- reorder(hcl, dists)
  return(hclust_olo)
}
expr = cnv@expr.data[,cnv@observation_grouped_cell_indices[[x]]]
# x = expr[,order.dendrogram(as.dendrogram(cnv@tumor_subclusters$hc[[x]]))]
hcl <- hclust(d = dist(t(expr)), method = "ward.D")
a <- cl_cb(hcl = hcl, mat = t(expr))
a <- order.dendrogram(as.dendrogram(hcl))
expr = expr[,a]
p1 <- pheatmap(t(plot.data[a,]), 
               color = colorRamp2(c(0, 1, 3), c("darkblue","white","darkred")), 
               cluster_cols = F, cluster_rows = T,
               show_rownames = F, show_colnames = F)



## Figure S4D ## 
plot6<-VlnPlot(seurat, features = c('MKI67'), pt.size=0, group.by="groups", ncol=3)+geom_boxplot(width=0.2,col="black",fill="#E1E3E6")
ggsave("备用/Figure S4/Vlnplot MKI67.pdf", plot = plot6, width = 10, height = 4)
plot7<-VlnPlot(seurat, features = c('HES1'), pt.size=0, group.by="groups", ncol=3)+geom_boxplot(width=0.2,col="black",fill="#E1E3E6")
ggsave("备用/Figure S4/Vlnplot HES1.pdf", plot = plot7, width = 10, height = 4)
plot8<-VlnPlot(seurat, features = c('SOX4'), pt.size=0, group.by="groups", ncol=3)+geom_boxplot(width=0.2,col="black",fill="#E1E3E6")
ggsave("备用/Figure S4/Vlnplot SOX4.pdf", plot = plot8, width = 10, height = 4)
plot9<-VlnPlot(seurat, features = c('ATP5ME'), pt.size=0, group.by="groups", ncol=3)+geom_boxplot(width=0.2,col="black",fill="#E1E3E6")
ggsave("备用/Figure S4/Vlnplot ATP5ME.pdf", plot = plot9, width = 10, height = 4)


## Figure S4E ## 
#####Epi-2与Epi-1比做GO富集
EpiSLG <- read.csv(file = "Epi cell_identify/GO/Epi2 vs 1.diff_genes_wilcox20211125V2.csv", header = T)
genes <- EpiSLG[[1]]

#GO富集分析
library(org.Hs.eg.db)
library(clusterProfiler)
#对于加载的注释库的使用，以上述为例，就直接在 OrgDb 中指定人（org.Hs.eg.db）或绵羊（sheep）
enrich.go <- enrichGO(gene = genes,  #基因列表文件中的基因名称
                      OrgDb = 'org.Hs.eg.db',  #指定物种的基因数据库，示例物种是绵羊（sheep）
                      keyType = 'SYMBOL',  #指定给定的基因名称类型，例如这里以 entrze id 为例
                      ont = 'ALL',  #可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                      pAdjustMethod = 'fdr',  #指定 p 值校正方法
                      pvalueCutoff = 0.05,  #指定 p 值阈值，不显著的值将不显示在结果中
                      qvalueCutoff = 0.2,  #指定 q 值阈值，不显著的值将不显示在结果中
                      readable = FALSE)

#例如上述指定 ALL 同时计算 BP、MF、CC，这里将它们作个拆分后输出
BP <- enrich.go[enrich.go$ONTOLOGY=='BP', ]
CC <- enrich.go[enrich.go$ONTOLOGY=='CC', ]
MF <- enrich.go[enrich.go$ONTOLOGY=='MF', ]

write.table(as.data.frame(BP), 'Epi cell_identify\\EPISL(V2)DEG\\go.BP Epi-2与1比.txt', sep = '\t', row.names = FALSE, quote = FALSE)
write.table(as.data.frame(CC), 'Epi cell_identify\\EPISL(V2)DEG\\go.CC Epi-2与1比.txt', sep = '\t', row.names = FALSE, quote = FALSE)
write.table(as.data.frame(MF), 'Epi cell_identify\\EPISL(V2)DEG\\go.MF Epi-2与1比.txt', sep = '\t', row.names = FALSE, quote = FALSE)

gene = bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(gene)

ego <- enrichKEGG(
  gene          = gene$ENTREZID,
  keyType     = "kegg",
  organism   = 'hsa',
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.1
)
write.table(as.data.frame(ego), 'Epi cell_identify\\EPISL(V2)DEG\\KEGG Epi-2与1比.txt', sep = '\t', row.names = FALSE, quote = FALSE)

# load package
library(ggplot2)
library(stringr)
# load data
go <- read.table(file = "Epi cell_identify\\EPI23DEG\\go.select5.txt",sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T)
go <- go[c(3:1,8:4,12:9),]
go$Term <- factor(go$Description, levels = go$Description)
#go$Term <- go$Description

go$GeneRatio <- as.numeric(sapply(strsplit(go$GeneRatio,"/",fixed = T),"[",1))/as.numeric(sapply(strsplit(go$GeneRatio,"/",fixed = T),"[",2))
ggplot(go, # you can replace the numbers to the row number of pathway of your interest
       aes(x = GeneRatio, y = Term)) + 
  geom_point(aes(size = Count, color = p.adjust)) +
  theme_bw(base_size = 10) +
  theme(axis.text = element_text(size = 10, colour = "black")) + 
  scale_colour_gradient(limits=c(0, 0.01), low="red") +
  ylab(NULL) + facet_grid(ONTOLOGY~.,scales = "free") + 
  ggtitle("GO pathway enrichment") + scale_y_discrete(labels=function(x) str_wrap(x, width=40))
ggsave(filename = "go enrichment20211209 Epi-2.pdf", width = 6,height = 4)



## Figure S4F ## 
standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
  outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata>halfwidth]=halfwidth
    outdata[outdata<(-halfwidth)]= -halfwidth
  }
  return(outdata)
}

seurat@meta.data <- Epi_metadata
mtx <- AverageExpression(object = seurat, group.by = "celltypes4")[[1]]
plot.gene <- c(paste0("ANXA", seq(1:10)), "S100P", "S100A2", "S100A4", "S100A6", "S100A8", "S100A11")
plot.celltype <- c("Epi-1", "Epi-2", "Transit amplifying cells", "Goblet cells", "Intermediate", "BEST4+ colonocytes", "Mature colonocytes")
plot.data <- mtx[plot.gene, plot.celltype]
plot.data <- standarize.fun(indata = plot.data, halfwidth = 2)
annCol <- data.frame(celltypes = factor(x = colnames(plot.data), levels = colnames(plot.data)),
                     row.names = colnames(plot.data),
                     stringsAsFactors = F)
annColors <- list()
annColors[["celltypes"]] = c("Epi-2" = "#E31A1C",
                             # "Epi-3" = "#E10C0D",
                             "Epi-1" = "#33A02C",
                             "Transit amplifying cells" = "#C9B1D5",
                             "Goblet cells" = "#FB9998",
                             "Intermediate" = "#FDBE6E",
                             "BEST4+ colonocytes" = "#A6CEE3",
                             "Mature colonocytes" = "#B2DF8A",
                             "Tuft cells" = "#FF7F00",
                             "Enteroendocrine cells" = "#6A3D9A") 

pdf("i:/genomicdata/jiaoda/plot/Heatmap/Epi1201.pdf", width = 796/100, height = 690/100,)
pheatmap(plot.data,
         cluster_cols = F, cluster_rows = F,
         # border_color = "white",
         # color = colorRampPalette(colors = c("#6353A3","white","#9A0038"))(100),
         color = viridis::inferno(64),
         annotation_col = annCol,
         annotation_colors = annColors,
         cellwidth = 15, cellheight = 15,
         angle_col = "45",
         treeheight_row = 20)
invisible(dev.off())

