##### Figure plotting #####
## Figure 5E ## 
sinfo <- read.csv("GSE39582 BRAF-mut20220121.csv",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
sinfo <- as.data.frame(sinfo[!is.na(sinfo$BRAF.Mutation),])
sinfo <- as.data.frame(sinfo[which(sinfo$BRAF.Mutation != "#N/A"),])

trm <- read.table("TRM signature.txt",sep = "\t",row.names = NULL,stringsAsFactors = F,header = T)
trm <- list(trm = trm$gene)
trm.gsva <- gsva(expr = as.matrix(expr[,rownames(sinfo)]),
                 trm,
                 method = "gsva")

sinfo$trm <- as.numeric(trm.gsva[1,])
kruskal.test(sinfo$notch~sinfo$BRAF.Mutation)

boxplot(sinfo$trm~sinfo$BRAF.Mutation)

sinfo$group <- ifelse(sinfo$BRAF.Mutation == "M","BRAF-mut","WT")
sinfo$group <- factor(sinfo$group, levels = c("WT","BRAF-mut"))
ggplot(data = sinfo,aes(x = group, #分组列名
                        y = trm, #连续变量列名
                        fill = group))+ #按分组填充颜???
  scale_fill_manual(values = jco[2:1]) + #用自定义颜色填充
  geom_violin(alpha = 0.4, position = position_dodge(width = .75),
              size = 0.8, color="white") +
  
  geom_point(shape = 21, size=2, # 点的性状和大???
             position = position_jitterdodge(), # 让点散开
             aes(color=group), alpha = 1) +
  scale_color_manual(values = jco[2:1]) + #用自定义颜色填充
  geom_boxplot(notch = TRUE, outlier.size = -1,
               color="black", lwd=0.8, alpha = 0.7) +
  # scale_y_continuous(
  #   breaks = c(2,4,6,8,10), 
  #   labels = c(2,4,6,8,10),
  #   limits = c(2,10) 
  # ) +
  #scale_x_discrete(labels= c("MT (n=55)","WT (n=431)")) + 
  theme_bw() + 
  ylab("TRM signature\nGSVA") +
  xlab("") +
  theme(axis.text.x = element_text(hjust = 0.5, size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.ticks = element_line(size=0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "none",
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10)) +
  stat_compare_means(method = "wilcox.test")
ggsave("boxplot for TRM signature in gse39582 regarding BRAF.pdf",width = 3,height = 4)



## Figure 5F ## 
gse39582.expr <- read.table("gse39582.expr.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
gse39582.surv1 <- read.table("GSE39582 BRAF-mut.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
gse39582.surv2 <- read.table("GSE39582 KRAS-mut.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
gse39582.surv <- rbind.data.frame(gse39582.surv1,gse39582.surv2)
rownames(gse39582.surv) <- gse39582.surv$GEOID
gse39582.expr <- gse39582.expr[,rownames(gse39582.surv)]
colnames(gse39582.surv)[8:11] <- c("RFS","RFS.time","OS","OS.time")

sam <- rownames(gse39582.surv[which(gse39582.surv$BRAF.Mutation == "M" | gse39582.surv$KRAS.Mutation == "M"),])
df <- data.frame(lyz = as.numeric(gse39582.expr["LYZ",sam]),
                 muc2 = as.numeric(gse39582.expr["MUC2",sam]),
                 OS = as.numeric(gse39582.surv[sam,"OS"]),
                 OS.time = as.numeric(gse39582.surv[sam,"OS.time"]),
                 PFI = as.numeric(gse39582.surv[sam,"RFS"]),
                 PFI.time = as.numeric(gse39582.surv[sam,"RFS.time"]))
coxph(Surv(OS.time,OS)~lyz, df)
coxph(Surv(PFI.time,PFI)~lyz, df)


signature <- rbind.data.frame(T.signature,B.signature,CD8.signature,Myeloid.signature)
cell.type <- unique(signature$cluster)
scsignature <- list()
for (i in cell.type) {
  scsignature[[i]] <- intersect(signature[which(signature$cluster == i),"gene"],rownames(gse39582.expr))
}

scsignature.gsva <- gsva(as.matrix(gse39582.expr),
                         scsignature,
                         method = "ssgsea")

tmp <- cbind.data.frame(gse39582.surv[,c("OS.time","OS","RFS.time","RFS")],as.data.frame(t(scsignature.gsva[,rownames(gse39582.surv)])))
Coxoutput.os <- Coxoutput.pfs <- NULL
for(k in 5:ncol(tmp)){
  gene <- colnames(tmp)[k]
  cox <- coxph(Surv(OS.time,OS) ~ scale(tmp[,k]), data = tmp)
  coxSummary = summary(cox)
  Coxoutput.os <- rbind.data.frame(Coxoutput.os,
                                   data.frame(gene = gene,
                                              HR = as.numeric(coxSummary$coefficients[,"exp(coef)"])[1],
                                              z = as.numeric(coxSummary$coefficients[,"z"])[1],
                                              pval = as.numeric(coxSummary$coefficients[,"Pr(>|z|)"])[1],
                                              lower = as.numeric(coxSummary$conf.int[,3][1]),
                                              upper = as.numeric(coxSummary$conf.int[,4][1]),
                                              row.names = gene,
                                              stringsAsFactors = F),
                                   stringsAsFactors = F)
  
  cox <- coxph(Surv(RFS.time,RFS) ~ scale(tmp[,k]), data = tmp)
  coxSummary = summary(cox)
  Coxoutput.pfs <- rbind.data.frame(Coxoutput.pfs,
                                    data.frame(gene = gene,
                                               HR = as.numeric(coxSummary$coefficients[,"exp(coef)"])[1],
                                               z = as.numeric(coxSummary$coefficients[,"z"])[1],
                                               pval = as.numeric(coxSummary$coefficients[,"Pr(>|z|)"])[1],
                                               lower = as.numeric(coxSummary$conf.int[,3][1]),
                                               upper = as.numeric(coxSummary$conf.int[,4][1]),
                                               row.names = gene,
                                               stringsAsFactors = F),
                                    stringsAsFactors = F)
}
write.table(Coxoutput.os, "new univariate cox results of os using cluster marker in all gse39582.txt",sep = "\t",row.names = F,col.names = T,quote = F)
write.table(Coxoutput.pfs, "new univariate cox results of rfs using cluster marker in all gse39582.txt",sep = "\t",row.names = F,col.names = T,quote = F)

# forest plot
red <- "#C3281D"
blue <- "#899ECE"
dat1 <- read.table("forestinput1.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
dat1[which(dat1$HR > 1 & dat1$pval < 0.05),"col"] <- red
dat1[which(dat1$HR < 1 & dat1$pval < 0.05),"col"] <- blue
dat1[which(dat1$pval >= 0.05),"col"] <- "black"
dat1$label <- paste0(round(dat1$HR,2)," (",round(dat1$lower,2),"-",round(dat1$upper,2),")")
dat1$label <- paste0(round(dat1$HR,2)," (",round(dat1$lower,2),"-",round(dat1$upper,2),")    ",format(round(dat1$pval,3),nsmall = 3))

dat2 <- read.table("forestinput2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
dat2[which(dat2$HR > 1 & dat2$pval < 0.05),"col"] <- red
dat2[which(dat2$HR < 1 & dat2$pval < 0.05),"col"] <- blue
dat2[which(dat2$pval >= 0.05),"col"] <- "black"
dat2$label <- paste0(format(round(dat2$HR,2),nsmall = 2)," (",format(round(dat2$lower,2),nsmall = 2),"-",format(round(dat2$upper,2),nsmall = 2),")    ",format(round(dat2$pval,3),nsmall = 3))

pdf("forestplot1.pdf", width = 8,height = 6)
par(bty="l", mgp = c(1.5,.33,0), mar=c(3,8,3,12), las=1, tcl=-.25, xpd = F)
plot(NULL,NULL,
     xlim = range(pretty(c(dat1$lower,dat1$upper))),
     ylim = c(1,nrow(dat1)),
     xlab = "Hazard ratio",
     ylab = "",
     xaxt = "n",
     yaxt = "n")
abline(v = 1, lwd = 2, lty = 2)
axis(side = 1, at = pretty(c(dat1$lower,dat1$upper)),labels = pretty(c(dat1$lower,dat1$upper)))
axis(side = 2, at = 1:nrow(dat1), labels = rev(rownames(dat1)))
axis(side = 4, at = 1:nrow(dat1), labels = rev(dat1$label),line = NA, col = NA)
for (i in 1:nrow(dat1)) {
  tmp <- dat1[i,]
  lines(c(tmp$lower,tmp$upper),c(nrow(dat1)-i + 1,nrow(dat1)-i + 1),
        col = tmp$col,
        lwd = 1.2)
  points(x = tmp$HR,
         y = nrow(dat1)-i + 1,
         col = tmp$col,
         pch = 19,
         cex = 1.5)
  lines(c(tmp$lower,tmp$lower),c(nrow(dat1)-i + 1 + 0.1,nrow(dat1)-i + 1-0.1),
        col = tmp$col,
        lwd = 1.2)
  lines(c(tmp$upper,tmp$upper),c(nrow(dat1)-i + 1 + 0.1,nrow(dat1)-i + 1-0.1),
        col = tmp$col,
        lwd = 1.2)
}
par("usr")
# [1]  0.468  1.332  0.600 11.400
par(xpd = T)
text(1.25,13,bquote("Hazard ratio (95% CI)    "~italic(p)~"-value"),adj = 0)
#Arrows(x0 = 1, y0 = 12,x1 = 0.9,y1 = 12, col = blue,)
invisible(dev.off())

pdf("forestplot2.pdf", width = 8,height = 6)
par(bty="l", mgp = c(1.5,.33,0), mar=c(3,8,3,12), las=1, tcl=-.25, xpd = F)
plot(NULL,NULL,
     xlim = range(pretty(c(dat2$lower,dat2$upper))),
     ylim = c(1,nrow(dat2)),
     xlab = "Hazard ratio",
     ylab = "",
     xaxt = "n",
     yaxt = "n")
abline(v = 1, lwd = 2, lty = 2)
axis(side = 1, at = pretty(c(dat2$lower,dat2$upper)),labels = pretty(c(dat2$lower,dat2$upper)))
axis(side = 2, at = 1:nrow(dat2), labels = rev(rownames(dat2)))
axis(side = 4, at = 1:nrow(dat2), labels = rev(dat2$label),line = NA, col = NA)
for (i in 1:nrow(dat2)) {
  tmp <- dat2[i,]
  lines(c(tmp$lower,tmp$upper),c(nrow(dat2)-i + 1,nrow(dat2)-i + 1),
        col = tmp$col,
        lwd = 1.2)
  points(x = tmp$HR,
         y = nrow(dat2)-i + 1,
         col = tmp$col,
         pch = 19,
         cex = 1.5)
  lines(c(tmp$lower,tmp$lower),c(nrow(dat2)-i + 1 + 0.1,nrow(dat2)-i + 1-0.1),
        col = tmp$col,
        lwd = 1.2)
  lines(c(tmp$upper,tmp$upper),c(nrow(dat2)-i + 1 + 0.1,nrow(dat2)-i + 1-0.1),
        col = tmp$col,
        lwd = 1.2)
}
par("usr")
# [1]  0.468  1.332  0.600 11.400
par(xpd = T)
text(1.35,11,bquote("Hazard ratio (95% CI)    "~italic(p)~"-value"),adj = 0)
#Arrows(x0 = 1, y0 = 12,x1 = 0.9,y1 = 12, col = blue,)
invisible(dev.off())



## Figure 5G ##
#scMetabolism
library(Seurat)
path <- "~/data/cwx/jiaoda/scMetabolism";setwd(path)

library(scMetabolism)
count <- readRDS("~/data/cwx/jiaoda/raw/countStro.rds")
KEGG<-sc.metabolism(countexp = as.matrix(count),
                    method = "VISION", 
                    imputation = F, 
                    ncores = 10, 
                    metabolism.type = "KEGG")
REACTOME<-sc.metabolism(countexp = as.matrix(count),
                        method = "VISION", 
                        imputation = F, 
                        ncores = 10, 
                        metabolism.type = "REACTOME")
saveRDS(list(KEGG = KEGG, REACTOME = REACTOME), "res.rds")