#####Figure S5
library(igraph)

# read regulon.csv --------------------------------------------------------


regulon.tab <- read.csv(file.path(data.path, "reg.csv"))
colnames(regulon.tab) <- paste0(regulon.tab[1,], regulon.tab[2,])
regulon.tab <- regulon.tab[-c(1,2), ]
regulon.weight <- lapply(1:nrow(regulon.tab), function(x){
  tmp = regulon.tab[x, c("TF", "TargetGenes")]
  target = gsub(pattern = "\\[|\\]|\\(|\\)|'", replacement = "", x = tmp$TargetGenes)
  target = unlist(strsplit(x = target, split = ","))
  target = matrix(target, ncol = 2, byrow = T)
  target = data.frame(
    row.names = NULL,
    "TF" = tmp["TF"],
    "target" = target[, 1],
    "weight" = as.numeric(target[, 2])
  )
})

regulon.weight <- do.call(rbind, regulon.weight)

### check whether the weight in different regulon is the same
# test <- split(x = regulon.weight$weight, f = paste0(regulon.weight$TF, regulon.weight$target))
# test <- lapply(test, function(x){length(unique(x))})
# summary(unlist(test));rm(test) 

regulon.weight <- regulon.weight[!duplicated(paste0(regulon.weight$target, "_", regulon.weight$TF)), ]
regulon.weight <- reshape2::recast(data = regulon.weight, formula = target ~ TF)
rownames(regulon.weight) <- regulon.weight$target
regulon.weight <- regulon.weight[, -1]
sum(!is.na(regulon.weight))
regulon.weight[is.na(regulon.weight)] <- 0

# adj.mat <- as.matrix(parallelDist::parallelDist(x = t(as.matrix(regulon.weight))))
adj.mat <- regulon.weight[!duplicated(paste0(regulon.weight$target, "_", regulon.weight$TF)), ]
KNN.mat <- FNN::get.knn(data = t(regulon.weight), k = floor(sqrt(ncol(regulon.weight))))
adj.mat.index <- as.data.frame(KNN.mat$nn.index)
adj.mat.index <- lapply(1:nrow(adj.mat.index), function(i){
  sapply(1:nrow(adj.mat.index), function(j){
    if (j %in% adj.mat.index[i, ]){
      adj.mat.index = 1
    }else{adj.mat.index = 0}
  })
})
adj.mat.index <- do.call(rbind, adj.mat.index)
adj.mat[adj.mat.index == 0] <- 0

# G <- igraph::graph.adjacency(adjmatrix = adj.mat,
#                              mode = "undirected", weighted = T, diag = F)
# 
# vcount(G);ecount(G)
# cluster <- cluster_louvain(G)
# plot(G, vertex.color=rainbow(3, alpha=0.6)[cluster$membership])

G <- leiden::leiden(object = adj.mat.index,
                    resolution_parameter = 2)

cluster.res <- data.frame(
  "TF" = colnames(regulon.weight),
  "cluster" = G
)
cluster.res <- dplyr::arrange(cluster.res, cluster.res$cluster)
openxlsx::write.xlsx(cluster.res, "superregulon.xlsx")

# regulon with different regulon ------------------------------------------

auc.mtx <- read.csv(file.path(data.path, "auc_mtx.csv"))
rownames(auc.mtx) <- auc.mtx$Cell
auc.mtx <- auc.mtx[, -1]
colnames(auc.mtx) <- gsub("\\...", "", colnames(auc.mtx))
rownames(auc.mtx) <- gsub("\\.", "-", rownames(auc.mtx))

load("i:/genomicdata/jiaoda/raw/Epi_metadata20211125V2.Rdat")
all(rownames(Epi_metadata)%in%rownames(auc.mtx))
auc.mtx <- auc.mtx[rownames(Epi_metadata), ]
auc.mtx <- scale(auc.mtx)
Epi1.AUC <- auc.mtx[rownames(Epi_metadata)[Epi_metadata$celltypes4 == "Epi-1"], ]
Epi2.AUC <- auc.mtx[rownames(Epi_metadata)[Epi_metadata$celltypes4 == "Epi-2"], ]
Epi_metadata <- subset(Epi_metadata, celltypes4 %in% c("Epi-1", "Epi-2"))
cell_info <- data.frame(
  row.names = rownames(Epi_metadata),
  ident = droplevels(Epi_metadata$celltypes4)
)

DE.regulon <- pbapply::pblapply(colnames(auc.mtx), function(TF){
  c(wilcox.test(Epi1.AUC[, TF], Epi2.AUC[, TF], conf.int = TRUE)$estimate,
    wilcox.test(Epi1.AUC[, TF], Epi2.AUC[, TF], conf.int = TRUE)$p.value)
})
DE.regulon <- do.call(rbind, DE.regulon)
DE.regulon



# network plot ------------------------------------------------------------

SuperRegulon <- data.frame(
  row.names = cluster.res$TF,
  "TF" = cluster.res$TF,
  "cluster" = cluster.res$cluster,
  "Epi-1" = round(colMeans(Epi1.AUC), 2),
  "Epi-2" = round(colMeans(Epi2.AUC), 2)
)

vertex.color <- colorRampPalette(c("#75A3D0","white","red"))(201)
names(vertex.color) <- round(seq(-1, 1, 0.01), digits = 2)

pdf("test3.pdf", width = 2300/100, height = 1200/100)
lapply(2:3, function(index){
  TF <- SuperRegulon$TF[SuperRegulon$cluster == index]
  TF.adj <- regulon.weight[TF, TF]
  rownames(TF.adj) <- TF
  TF.adj[is.na(TF.adj)] <- 0
  plotdata <- graph_from_adjacency_matrix(as.matrix(TF.adj), weighted = TRUE, diag = F, mode = "directed")
  
  Epi1.col = vertex.color[match(SuperRegulon[TF,]$Epi.1, names(vertex.color))]
  Epi2.col = vertex.color[match(SuperRegulon[TF,]$Epi.2, names(vertex.color))]
  par(mfrow=c(1,2), cex = 1.1)
  # layout(matrix(c(1:3), ncol = 3), widths = c(4,1,4))
  plot(plotdata, 
       layout = layout.circle, 
       vertex.color = Epi1.col, 
       edge.width = sapply((E(plotdata)$weight/2.5), function(x){min(x, 10)}), 
       edge.arrow.size = sapply((E(plotdata)$weight/20), function(x){min(x, 3)}),
       edge.arrow.size = 0.7,
       vertex.shape="rectangle",
       vertex.size = setNames(object = nchar(TF)*3, nm = TF),
       vertex.size2 = 6,
       main = "Epi-1")
  plotrix::color.legend(1.3, -0.3, 1.35, 0.3, legend = c(-1, 0, 1), rect.col = vertex.color, gradient = "y")
  plot(plotdata, 
       label.font = "2",
       layout = layout.circle, 
       vertex.color = Epi2.col, 
       edge.width = sapply((E(plotdata)$weight/2.5), function(x){min(x, 10)}), 
       edge.arrow.size = sapply((E(plotdata)$weight/20), function(x){min(x, 3)}),
       # edge.arrow.size = 0.7,
       vertex.shape="rectangle",
       vertex.size = setNames(object = nchar(TF)*3, nm = TF),
       vertex.size2 = 6,
       main = "Epi-2")
})
invisible(dev.off())

adj <- read.table(file.path(data.path, "adj.tsv"), sep = "\t", header = T)
adj1 <- reshape2::recast(data = adj, formula = target ~ TF)
rownames(adj1) <- adj1$target; adj1 <- adj1[, -1]


# # Trash -------------------------------------------------------------------
# 
# 
# adj <- TF.adj
# adj$target <- rownames(adj)
# adj <- reshape2::melt(data = adj, id.vars = "target", variable.name = "TF", value.name = "adjacency")
# adj <- data.frame("TF" = adj$TF, "target" = adj$target, "adjacency" = adj$adjacency)
# write.csv(adj, "RegulonNetwork.csv")
# a <- matrix(
#   c(0, 5, 10, 8,
#     5, 0, 6, 7,
#     10, 6, 0, 4,
#     8, 7, 4, 0),
#   ncol = 4
# )
# leiden::leiden(object = a, resolution_parameter = 1.5)
# 
