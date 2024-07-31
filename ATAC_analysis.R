setwd("/nfs4/chaozhang/proj/Embryo/Yota/BRG1_ATAC_0901/script")
library(edgeR)
library(ggplot2)
library(pheatmap)
library(reshape2)
library(ggvenn)
library(pheatmap)
library(patchwork)

expr = read.table("../peaks/ATAC.peaks.read.bed", header = F)

rownames(expr) = paste(expr$V1, expr$V2, expr$V3)

expr = expr[, -c(1:10)]
#colnames(expr) = gsub(".mm.sort.rmdup.bam", "", list.files("../mapping", ".mm.sort.rmdup.bam$"))
colnames(expr) = c("Ctrl_rep1", "Ctrl_rep2", "100n_rep1","100n_rep2","10n_rep1","10n_rep2", "1n_rep1","1n_rep2","3n_rep1","3n_rep2","0.3n_rep1","0.3n_rep2", "WT_rep1", "WT_rep2")

grp = stringr::str_replace(colnames(expr), "_rep[1|2]", "")

y <- DGEList(counts=expr, group=grp)
y$samples

keep <- rowMeans(expr) > 10
y <- y[keep, , keep.lib.sizes=FALSE]

y <- calcNormFactors(y)
y$samples

normalized_counts = cpm(y, normalized.lib.sizes=T, log=F, prior.count = 1)
normalized_counts_group = cpmByGroup(y, normalized.lib.sizes=T, log=F, prior.count = 1)

write.csv(normalized_counts, "ATAC_normalized_counts.csv", quote = F)
write.csv(normalized_counts_group, "ATAC_normalized_counts_group.csv", quote = F)

## 
DM = read.table("../BW/size.factor.txt")
DM$size = DM$V2/DM$V3

library(sva)
batch = stringr::str_replace(colnames(expr), ".*_rep", "")
normalized_counts_removeBatch = ComBat(dat=normalized_counts, batch=batch)


normalized_counts_log = log2(normalized_counts+1)
plot(normalized_counts_log[,1], normalized_counts_log[,2])
cor(normalized_counts_log[,1], normalized_counts_log[,2], method = "spearman")
smoothScatter(normalized_counts_log[,1], normalized_counts_log[,2], xlab = "Ctrl_rep1", ylab = "Ctrl_rep2", main = "cor: 0.87")

cor(normalized_counts_log[,3], normalized_counts_log[,4], method = "spearman")
smoothScatter(normalized_counts_log[,3], normalized_counts_log[,4], xlab = "100n_rep1", ylab = "100n_rep2", main = "cor: 0.85")


df.cor = cor(normalized_counts, method = "spearman")
diag(df.cor) <- .88
pdf("../figure/correlation_peaks.pdf", 6, 5)
pheatmap::pheatmap(df.cor, border_color = NA)
dev.off()


pca = prcomp(t(normalized_counts), scale = T)
pca = as.data.frame(pca$x)

pca = cbind(grp, pca)
pca$grp = factor(pca$grp, levels = c("WT", "Ctrl", "0.3n", "1n", "3n", "10n", "100n"))
ggplot(pca, aes(PC1, PC2)) +
  geom_point(aes(color = grp), size=4)+
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(colour = "black")) 
ggsave("../figure/Samples_PCA.pdf", width = 6, height = 5)



#

# DE Peaks
expr1 = expr[,c(3,4,1,2)]
cond = matrix(c("100n", "rep1",
                "100n", "rep2",
                "Ctrl", "rep1",
                "Ctrl", "rep2"), ncol = 2, byrow = T)
rownames(cond) = colnames(expr1)
colnames(cond) = c("condition", "rep")
cond = as.data.frame(cond)

y <- DGEList(counts=expr1, group=cond$condition)
y$samples

keep <- rowMeans(expr) > 10
y <- y[keep, , keep.lib.sizes=FALSE]

y <- calcNormFactors(y)
y$samples
cond$condition = factor(cond$condition, levels = c("Ctrl", "100n"))
design <- model.matrix(~cond$condition + cond$rep)

y <- estimateDisp(y, design)

fit <- glmFit(y,design)
lrt <- glmLRT(fit, coef=2)
DEGs = as.data.frame(topTags(lrt, n = nrow(y$counts)))
DEGs.sig = DEGs[DEGs$FDR < 0.05 & abs(DEGs$logFC) > 1, ]

plotMD(lrt)

DEGs$sig = "none"
DEGs$sig[DEGs$logFC > 0 & rownames(DEGs) %in% rownames(DEGs.sig)] = "Up"
DEGs$sig[DEGs$logFC < 0 & rownames(DEGs) %in% rownames(DEGs.sig)] = "Down"

ggplot(DEGs, aes(logFC, -log10(PValue), col = sig)) +
  geom_point() +
  theme_classic(base_size = 20) + xlab("logFC: dTag_100n/Ctrl") +
  scale_color_manual(values = c("blue", "gray90", "red"))
ggsave("../figure/DEG_dTag100n_Ctrl.scatterplot.pdf", width = 8, height = 10)

table(DEGs$sig)


write.table(DEGs, "DEG_dTag100n_Ctrl.txt", quote = F, sep = "\t")
DEGs = read.table("DEG_dTag100n_Ctrl.txt", header = T, sep = "\t", row.names = 1)


normalized_counts1 = normalized_counts[rownames(DEGs.sig), ]
breaksList = seq(-2, 2, by = .01)
pheatmap::pheatmap(normalized_counts1[, ], show_rownames = F, 
                   scale = "row",
                   #cluster_cols = F,
                   #clustering_method = "ward.D",
                   color = colorRampPalette(c("blue", "white", "red"))(length(breaksList)),
                   breaks = breaksList)


normalized_counts_group1 = normalized_counts_group[rownames(DEGs.sig), ]
pheatmap::pheatmap(normalized_counts_group1[, c("WT", "Ctrl", "0.3n", "1n", "3n", "10n", "100n")],
                   show_rownames = F, scale = "row",
                   cluster_cols = F,
                   clustering_method = "ward.D",
                   color = colorRampPalette(c("blue", "white", "red"))(length(breaksList)),
                   breaks = breaksList)



###
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

files = list(None = "./DEG_dTag100n_Ctrl.non.bed",
             Up = "./DEG_dTag100n_Ctrl.up.bed",
             Down = "./DEG_dTag100n_Ctrl.down.bed",
             BRG1 = "../../BRG1_CutRun_Nov2023/peaks/Ctr/Ctr_peaks.narrowPeak")

peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb,
                       annoDb="org.Mm.eg.db",
                       tssRegion=c(-1000, 1000), verbose=FALSE)

plotAnnoBar(peakAnnoList)



genes= lapply(peakAnnoList, function(i) as.data.frame(i)$SYMBOL)
vennplot(genes)



# total peak number


pn = data.frame(sample = c("WT", "Ctrl",  "1n", "3n", "10n", "100n"),
                peak_n = c(40995, 43212, 40673, 40164, 35207,37441))
pn$sample = factor(pn$sample, level = c("WT", "Ctrl", "1n", "3n", "10n", "100n"))
ggplot(pn, aes(sample, peak_n, fill = sample)) +
  geom_bar(stat = "identity") +
  theme_classic()+
  scale_fill_manual(values = c("red4", "red3", "red2", "red", "orange", "yellow2"))
ggsave("../figure/Peak_number.pdf", width = 8, height = 7)


pn = data.frame(sample = c("None", "Up", "Down"),
                peak_n = c(7413/75941, 50/553, 2621/5863))

ggplot(pn, aes(sample, peak_n, fill = sample)) +
  geom_bar(stat = "identity") +
  theme_classic()+ 
  ylab("Percent of ATAC Peaks Overlap with Brg1 Peaks")+ xlab("ATAC-seq Peaks")+
  scale_fill_manual(values = c( "blue","gray60", "red"))




## base on Brg1 peaks

#Brg1_peak = read.table("../peaks/Brg1.ATAC_reads.bed", header = F)
Brg1_peak = read.table("../peaks/Brg_new.ATAC_reads.bed", header = F)
Brg1_peak = Brg1_peak[!duplicated(Brg1_peak$V4),]
rownames(Brg1_peak) = paste(Brg1_peak$V1, Brg1_peak$V2, Brg1_peak$V3)
Brg1_peak = Brg1_peak[, -c(1:9)]

colnames(Brg1_peak) = c("Ctrl_rep1", "Ctrl_rep2", "100n_rep1","100n_rep2","10n_rep1","10n_rep2", "1n_rep1","1n_rep2","3n_rep1","3n_rep2", "0.3n_rep1","0.3n_rep2","WT_rep1", "WT_rep2")

grp = stringr::str_replace(colnames(Brg1_peak), "_rep[1|2]", "")

y1 <- DGEList(counts=Brg1_peak, group=grp)
y1$samples


y1 <- calcNormFactors(y1, method = "upperquartile", p = 0.95)


##
out_Brg1_peak = read.table("../peaks/Out_Brg1_ATAC.ATAC_reads.bed", header = F)
out_Brg1_peak = out_Brg1_peak[, -c(1:10)]
colnames(out_Brg1_peak) = c("Ctrl_rep1", "Ctrl_rep2", "100n_rep1","100n_rep2","10n_rep1","10n_rep2", "1n_rep1","1n_rep2","3n_rep1","3n_rep2", "0.3n_rep1","0.3n_rep2","WT_rep1", "WT_rep2")
y0 <- DGEList(counts=out_Brg1_peak, group=grp)
y0 <- calcNormFactors(y0, method = "upperquartile", p = 0.9)
y0$samples$norm.factors
##

#y1$samples$norm.factors = y0$samples$norm.factors
y1$samples

y1_normalized_counts = cpm(y1, normalized.lib.sizes=T, log=F, prior.count = 1)
write.csv(y1_normalized_counts, "ATAC_at_brg1_normalized_counts_unlog.csv", quote = F)
y1_normalized_counts_group = cpmByGroup(y1, normalized.lib.sizes=T, log=F, prior.count = 1,
                                        lib.size = y$samples$lib.size)
write.csv(y1_normalized_counts_group, "ATAC_at_brg1_normalized_counts_group_unlog.csv", quote = F)


y1_normalized_counts = cpm(y1, normalized.lib.sizes=T, log=T, prior.count = 1)
y1_normalized_counts_group = cpmByGroup(y1, normalized.lib.sizes=T, log=T, prior.count = 1,
                                        lib.size = y$samples$lib.size)

write.csv(y1_normalized_counts, "ATAC_at_brg1_normalized_counts.csv", quote = F)


ggplot(as.data.frame(y1_normalized_counts), aes(Ctrl_rep1, `100n_rep1`)) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") +
  #geom_density_2d_filled(alpha = 1, colour = "white")+
  geom_point(size = .1, alpha = .3) +
  geom_abline(slope = 1) + xlim(1,10) + ylim(1,10) +
  scale_fill_viridis_c() + theme_classic()

ggplot(as.data.frame(y1_normalized_counts), aes(Ctrl_rep2, `100n_rep2`)) +
  #stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") +
  #geom_density_2d_filled(alpha = 1, colour = "white")+
  geom_point(size = .1, alpha = .3) +
  geom_abline(slope = 1) + xlim(0,10) + ylim(0,10) +
  scale_fill_viridis_c() + theme_classic()


p0 = ggplot(as.data.frame(y1_normalized_counts_group), aes(Ctrl, `0.3n`)) +
  #stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") +
  geom_point(size = .1, alpha = .3) +
  xlim(1,10) + ylim(1,10) +geom_abline(slope = 1) +
  scale_fill_viridis_c() + theme_classic()
p1 = ggplot(as.data.frame(y1_normalized_counts_group), aes(Ctrl, `1n`)) +
  #stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") +
  geom_point(size = .1, alpha = .3) +
  xlim(1,10) + ylim(1,10) +geom_abline(slope = 1) +
  scale_fill_viridis_c() + theme_classic()
p2 = ggplot(as.data.frame(y1_normalized_counts_group), aes(Ctrl, `3n`)) +
  #stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") +
  geom_point(size = .1, alpha = .3) +
  xlim(1,10) + ylim(1,10)+geom_abline(slope = 1) +
  scale_fill_viridis_c() + theme_classic()
p3 = ggplot(as.data.frame(y1_normalized_counts_group), aes(Ctrl, `10n`)) +
  #stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") +
  geom_point(size = .1, alpha = .3) +
  xlim(1,10) + ylim(1,10) +geom_abline(slope = 1) +
  scale_fill_viridis_c() + theme_classic()
p4 = ggplot(as.data.frame(y1_normalized_counts_group), aes(Ctrl, `100n`)) +
  #stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") +
  geom_point(size = .1, alpha = .3) +
  xlim(1,10) + ylim(1,10) +geom_abline(slope = 1) +
  scale_fill_viridis_c() + theme_classic()
p0|p1|p2|p3|p4
ggsave("../figure/ATAC_at_BRG1_peaks.pair-wise_scatterplot.png", width = 25, height = 5)


#


df.cor = cor(c, method = "spearman")
diag(df.cor) <- .95
pheatmap::pheatmap(df.cor, border_color = NA)



pca = prcomp(t(y1_normalized_counts), scale = T)
pca = as.data.frame(pca$x)

pca = cbind(grp, pca)
pca$grp = factor(pca$grp, levels = c("WT", "Ctrl", "0.3n", "1n", "3n", "10n", "100n"))

pca$sample = rownames(pca)
ggplot(pca, aes(PC1, PC2, label = sample)) +
  geom_point(aes(fill = grp),shape = 21, colour = "white", size=6)+
  cowplot::theme_cowplot() +
  scale_fill_brewer(palette="Dark2") +
  geom_text_repel(min.segment.length = Inf, box.padding = 0.5)

ggsave("../figure/ATAC_on_Brg1_PCA.pdf", width = 6, height = 5)






# DE Peaks in Brg1 loc
DAR_call = function(c1, c2){
  
  
  cond = matrix(c(c1, "rep1",
                  c1, "rep2",
                  c2, "rep1",
                  c2, "rep2"), ncol = 2, byrow = T)
  
  rownames(cond) = paste(cond[,1], cond[,2], sep = "_")
  colnames(cond) = c("condition", "rep")
  cond = as.data.frame(cond)
  
  Brg1_peak1 = Brg1_peak[,rownames(cond)]
  
  y2 <- DGEList(counts=Brg1_peak1, group=cond$condition)
  y2$samples
  
  #keep <- rowMeans(Brg1_peak1) > 10
  #y2 <- y2[keep, , keep.lib.sizes=FALSE]
  
  y2 <- calcNormFactors(y2, method = "upperquartile", p = 0.95)

  cond$condition = factor(cond$condition, levels = c(c2, c1))
  design <- model.matrix(~ cond$condition)
  
  y2 <- estimateDisp(y2, design)
  
  fit <- glmFit(y2,design)
  lrt <- glmLRT(fit)
  Brg1_DEGs = as.data.frame(topTags(lrt, n = nrow(y2$counts)))
  Brg1_DEGs.sig = Brg1_DEGs[Brg1_DEGs$FDR < 0.05, ]
  
  Brg1_DEGs$sig = "none"
  Brg1_DEGs$sig[Brg1_DEGs$logFC > 0 & rownames(Brg1_DEGs) %in% rownames(Brg1_DEGs.sig)] = "Up"
  Brg1_DEGs$sig[Brg1_DEGs$logFC < 0 & rownames(Brg1_DEGs) %in% rownames(Brg1_DEGs.sig)] = "Down"
  
  
  write.table(Brg1_DEGs, paste0("Brg1_DEGs_", c1, "_", c2, ".txt"), quote = F, sep = "\t")
  
  
   ggplot(as.data.frame(y1_normalized_counts_group), aes(get(c2), get(c1))) +
    geom_point(size = .1, alpha = .3) +
    geom_point(data = as.data.frame(y1_normalized_counts_group[rownames(Brg1_DEGs)[Brg1_DEGs$sig == "Down"],]), size = .5, color = "blue") +
    geom_point(data = as.data.frame(y1_normalized_counts_group[rownames(Brg1_DEGs)[Brg1_DEGs$sig == "Up"],]), size = .5, color = "red") +
    geom_abline(slope = 1) + xlim(1,10) + ylim(1,10) + xlab(c2) + ylab(c1)+
    scale_fill_viridis_c() + theme_classic(base_size = 15) 
}

p1 = DAR_call("0.3n", "Ctrl")
p2 = DAR_call("1n", "Ctrl")
p3 = DAR_call("3n", "Ctrl")
p4 = DAR_call("10n", "Ctrl")
p5 = DAR_call("100n", "Ctrl")
p1|p2|p3|p4|p5
ggsave("../figure/ATAC_at_BRG1_peaks.DAR.pair-wise_scatterplot.pdf", width = 20, height = 4)
ggsave("../figure/ATAC_at_BRG1_peaks.DAR.pair-wise_scatterplot.png", width = 25, height = 5)


# 
DAR_plot = function(c1, c2){
  
  
  cond = matrix(c(c1, "rep1",
                  c1, "rep2",
                  c2, "rep1",
                  c2, "rep2"), ncol = 2, byrow = T)
  
  rownames(cond) = paste(cond[,1], cond[,2], sep = "_")
  colnames(cond) = c("condition", "rep")
  cond = as.data.frame(cond)
  
  Brg1_peak1 = Brg1_peak[,rownames(cond)]
  
  y2 <- DGEList(counts=Brg1_peak1, group=cond$condition)
  y2$samples
  
  #keep <- rowMeans(Brg1_peak1) > 10
  #y2 <- y2[keep, , keep.lib.sizes=FALSE]
  
  y2 <- calcNormFactors(y2, method = "upperquartile", p = 0.95)
  
  cond$condition = factor(cond$condition, levels = c(c2, c1))
  design <- model.matrix(~ cond$condition)
  
  y2 <- estimateDisp(y2, design)
  
  fit <- glmFit(y2,design)
  lrt <- glmLRT(fit)
  Brg1_DEGs = as.data.frame(topTags(lrt, n = nrow(y2$counts)))
  Brg1_DEGs.sig = Brg1_DEGs[Brg1_DEGs$FDR < 0.05, ]
  
  Brg1_DEGs$sig = "none"
  Brg1_DEGs$sig[Brg1_DEGs$logFC > 0 & rownames(Brg1_DEGs) %in% rownames(Brg1_DEGs.sig)] = "Up"
  Brg1_DEGs$sig[Brg1_DEGs$logFC < 0 & rownames(Brg1_DEGs) %in% rownames(Brg1_DEGs.sig)] = "Down"
  
  df = as.data.frame(y1_normalized_counts_group)
  df$name = rownames(df)
  
  df$ann = Brg1_ann[match(df$name, Brg1_ann$name), "V4"]

  ggplot(na.omit(df), aes(get(c2), get(c1))) +
    geom_point(size = .1, alpha = .3) +
    geom_point(data = df[rownames(Brg1_DEGs)[Brg1_DEGs$sig == "Down"],], size = .5, color = "blue") +
    geom_point(data = df[rownames(Brg1_DEGs)[Brg1_DEGs$sig == "Up"],], size = .5, color = "red") +
    geom_abline(slope = 1) + xlim(1,10) + ylim(1,10) + xlab(c2) + ylab(c1)+
    scale_fill_viridis_c() + theme_classic(base_size = 15) +
    facet_wrap(.~ann, ncol = 1)
}

Brg1_ann = read.table("Brg1_peak_ann.all.bed")
Brg1_ann = unique(Brg1_ann)
Brg1_ann$name = paste(Brg1_ann$V1, Brg1_ann$V2, Brg1_ann$V3)

p1 = DAR_plot("0.3n", "Ctrl")
p2 = DAR_plot("1n", "Ctrl")
p3 = DAR_plot("3n", "Ctrl")
p4 = DAR_plot("10n", "Ctrl")
p5 = DAR_plot("100n", "Ctrl")
p1|p2|p3|p4|p5

ggsave("../figure/ATAC_at_BRG1_peaks.DAR.3elements.pair-wise_scatterplot.pdf", width = 20, height = 12)
ggsave("../figure/ATAC_at_BRG1_peaks.DAR.3elements.pair-wise_scatterplot.png", width = 25, height = 15)



for (i in list.files(".", "Brg1_DEGs_.*.txt")){
  f = read.table(i, sep = "\t")
  write.table(stringr::str_replace_all(rownames(f), " ", "\t"), stringr::str_replace(i, ".txt", ".bed"),
              quote = F, row.names = F, col.names = F)
  
  ff = f[f$sig == "Down", ]
  write.table(stringr::str_replace_all(rownames(ff), " ", "\t"), stringr::str_replace(i, ".txt", ".Down.bed"),
              quote = F, row.names = F, col.names = F)
  
  f = f[f$sig == "Up", ]
  write.table(stringr::str_replace_all(rownames(f), " ", "\t"), stringr::str_replace(i, ".txt", ".Up.bed"),
              quote = F, row.names = F, col.names = F)
}

pn = data.frame(group = rep(c( "up", "down"), each=5),
  sample = rep(c("0.3n", "1n", "3n", "10n", "100n"), 2),
  number = -c(2,18,107,211,269,
              -94,-2640,-5391,-5043,-6989))

pn$sample = factor(pn$sample, levels = c("0.3n", "1n", "3n", "10n", "100n"))
ggplot(pn, aes(sample, number, fill = group)) +
  geom_bar(stat = "identity", position = "identity") +
  theme_classic(base_size = 15) +
  scale_fill_manual(values = c("blue","red")) +
  ylab("number of DARs")
ggsave("../figure/DARs_numbers.pdf", width = 6, height = 5)



# Brg1 binding independent/dependent regions
DAR_100n = read.table("Brg1_DEGs_100n_Ctrl.txt", sep = "\t")

DAR_100n_feature = list(dependent = makeGRangesFromDataFrame(as.data.frame(stringr::str_split_fixed(rownames(DAR_100n[DAR_100n$sig == "Down",]), " ", 3)),
                                     seqnames.field = "V1", start.field = "V2", end.field = "V3"),
     independent =  makeGRangesFromDataFrame(as.data.frame(stringr::str_split_fixed(rownames(DAR_100n[DAR_100n$sig == "none",]), " ", 3)),
                              seqnames.field = "V1", start.field = "V2", end.field = "V3"))

peakAnnoList <- lapply(DAR_100n_feature, annotatePeak, TxDb=txdb,
                       annoDb="org.Mm.eg.db",
                       tssRegion=c(-1000, 1000), verbose=FALSE)

plotAnnoBar(peakAnnoList)

a = as.data.frame(peakAnnoList[[1]])

plotTSSbar = function(a){
  a = a[abs(a$distanceToTSS) < 50000, ]
  a$distanceToTSS = round( a$distanceToTSS / 1000)
  ggplot(a, aes(distanceToTSS))+
    geom_bar(aes(y = (..count..)/sum(..count..))) + ylim(0,0.3)
}

# gene expr
normalized_rpkm = read.table("/nfs4/chaozhang/proj/Embryo/Yota/RNAseq_Apr2024/script/normalized_rpkm.tab")

genes = lapply(peakAnnoList, function(i){
  g = as.data.frame(i)
  g = g[abs(g$distanceToTSS) < 1000, ]
  g = g$SYMBOL
  unique(g)
})

genes= lapply(peakAnnoList, function(i) {
  i = as.data.frame(i)
  i = i[i$annotation == "Promoter (<=1kb)", ]
  unique(i$SYMBOL)
})

genes[[1]] = setdiff(genes[[1]], genes[[2]])
genes[[2]] = setdiff(genes[[2]], genes[[1]])


breaksList = seq(0, 2, by = .1)
pheatmap(log(normalized_rpkm[as.character(unlist(genes)), c(11,12,1,2,7,8,9,10,5,6,3,4)]),
         scale = "row", 
         gaps_row = length(genes[[1]]),
         show_rownames = F, cluster_cols = F, cluster_rows = F,
         color = pals::coolwarm(length(breaksList)),
         breaks = breaksList)

boxplot(log(c[genes[[1]], c(11,12,1,2,7,8,9,10,5,6,3,4)]+1), ylim=c(0,3))
boxplot(log(normalized_rpkm[genes[[2]], c(11,12,1,2,7,8,9,10,5,6,3,4)]+1), ylim=c(0,3))

colMeans((na.omit(normalized_rpkm[genes[[1]], c(11,12,1,2,7,8,9,10,5,6,3,4)]+1)))
colMeans((na.omit(normalized_rpkm[genes[[2]], c(11,12,1,2,7,8,9,10,5,6,3,4)]+1)))


gene_dependent = data.frame(gene = unlist(genes),
                            group = c(rep("dependent", length(genes[[1]])), rep("independent", length(genes[[2]]))) )
gene_dependent = cbind(gene_dependent, normalized_rpkm[gene_dependent$gene, ])
gene_dependent = na.omit(gene_dependent)

gene_dependent.m = melt(gene_dependent)
gene_dependent.m$variable = sub("_total_rep[1|2]", "", gene_dependent.m$variable)
gene_dependent.m$variable = factor(gene_dependent.m$variable, 
                                   levels = c("Brg1_Ctrl", "Brg1_0.3n", "Brg1_1n", "Brg1_3n", "Brg1_10n", "Brg1_100n"))


ggplot(gene_dependent.m, aes(variable, log(value+1), fill = variable))+
  geom_boxplot(outlier.size = -1) + ylim(0,5)+
  facet_wrap(.~group, ncol = 1) +
  scale_fill_brewer(palette="YlOrRd", direction=-1) +
  cowplot::theme_cowplot()
ggsave("../figure/Gene_expr_BRG1_dependent.pdf", width = 12, height = 8)



lapply(c("Brg1_0.3n", "Brg1_1n", "Brg1_3n", "Brg1_10n", "Brg1_100n"), function(i){
  t.test(gene_dependent.m[gene_dependent.m$group == "dependent" & gene_dependent.m$variable == "Brg1_Ctrl", "value"],
         gene_dependent.m[gene_dependent.m$group == "dependent" & gene_dependent.m$variable == i, "value"], paired = T)
})


gene_dependent.m = melt(gene_dependent)
gene_dependent.m$variable = factor(gene_dependent.m$variable, 
                                   levels = colnames(normalized_rpkm)[c(11,12,1,2,7,8,9,10,5,6,3,4)])

ggplot(gene_dependent.m[gene_dependent.m$variable %in% c("Brg1_Ctrl_total_rep1", "Brg1_Ctrl_total_rep2", "Brg1_100n_total_rep1", "Brg1_100n_total_rep2"), ],
       aes(variable, log(value+1), fill = variable))+
  geom_boxplot(outlier.size = -1) + ylim(0,5)+
  facet_wrap(.~group, ncol = 1) +
  scale_fill_manual(values = c("#BD0026","#BD0026", "#FFFFB2", "#FFFFB2")) +
  cowplot::theme_cowplot()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("") + ylab("log2(expr+1)")
ggsave("../figure/Gene_expr_BRG1_dependent_100n.pdf", width = 8, height = 8)



gene_dependent.in = gene_dependent[gene_dependent$group == "dependent", ]
t.test(gene_dependent.in[,"Brg1_Ctrl_total_rep1"], gene_dependent.in[,"Brg1_100n_total_rep1"], paired = T)
t.test(gene_dependent.in[,"Brg1_Ctrl_total_rep2"], gene_dependent.in[,"Brg1_100n_total_rep2"], paired = T)

# GO
library(clusterProfiler)
compGO <- compareCluster(geneCluster   = genes,
                         fun           = "enrichGO",
                         keyType = "SYMBOL",
                         ont = "BP",
                         OrgDb = org.Mm.eg.db,
                         pAdjustMethod = "BH")
dotplot(compGO, showCategory = 30, title = "GO Enrichment Analysis")
ggsave("../figure/Gene_expr_BRG1_dependent_GOenrich.pdf", width = 8, height = 6)


# chrHMM
HMM = lapply(c("Down", "none", "shuff"), function(i){
  f = read.table(paste0("./chromHMM/Brg1_DEGs_100n_Ctrl_",i,".chrHMM.bed"))
  n = nrow(read.table(paste0("Brg1_DEGs_100n_Ctrl.",i,".bed")))
  table(f$V7)/n
})

df = as.data.frame(do.call(cbind, HMM))
colnames(df) = c("dependent", "independent", "random")

pdf("../figure/Brg1_independent_chrHMM.pdf", 5, 5)
pheatmap::pheatmap(df[c(1,4:10,2,3),],cluster_rows = F, cluster_cols = F,
                   display_numbers = T,
                   border_color = "white",
                   number_format = "%.2f", 
                   fontsize = 15,
                   number_color = "black",
                   color = colorRampPalette(c("white", "yellow", "red"))(100))
dev.off()

df_enrich = df
df_enrich$dependent = (df$dependent/df$random)
df_enrich$independent = (df$independent/df$random)
df_enrich$Elements = rownames(df_enrich)
df_enrich$Elements = stringr::str_split_fixed(df_enrich$Elements, "_", 2)[,2]
df_enrich = df_enrich[c(3,5,8,9),]
df_enrich.m = melt(df_enrich[,-3])
ggplot(df_enrich.m, aes(Elements, value, fill = variable)) +
  geom_bar(stat = "identity",position = "dodge") +
  cowplot::theme_cowplot()+ ylab("Enrichment: obs/random")+
  scale_fill_manual(values = c("purple", "orange")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
ggsave("../figure/Brg1_independent_enrich_chrHMM.pdf", width = 5, height = 4)

# width
peak_width = data.frame(group =  c(rep("dependent", length(DAR_100n_feature[[1]])), rep("independent", length(DAR_100n_feature[[2]]))),
                        width = c(DAR_100n_feature[[1]]@ranges@width, DAR_100n_feature[[2]]@ranges@width))


ggplot(peak_width, aes(width, color = group))+
  stat_ecdf(geom = "step") + xlim(100, 1000)

ggplot(peak_width, aes(group, log10(width), color = group))+
  geom_boxplot() + 
  cowplot::theme_cowplot()


#motif
top10 = lapply(list.files("/nfs4/chaozhang/proj/Embryo/Yota/BRG1_ATAC_0901/script/motif_dependent", "knownResults.txt", recursive =T, full.names = T), function(i){
  f = read.table(i,sep="\t", stringsAsFactors = F, comment.char = "*", header = T)
  f$Motif.Name[1:20]
})
top10 = unique(unlist(top10))

top10_d = lapply(list.files("/nfs4/chaozhang/proj/Embryo/Yota/BRG1_ATAC_0901/script/motif_dependent", "knownResults.txt", recursive =T, full.names = T), function(i){
  f = read.table(i,sep="\t", stringsAsFactors = F, comment.char = "*", header = T)
  f = f[!duplicated(f$Motif.Name), ]
  rownames(f) = f$Motif.Name
  f = f[top10, c("P.value", "X..of.Target.Sequences.with.Motif")]
  f$X..of.Target.Sequences.with.Motif = as.numeric(sub("%", "", f$X..of.Target.Sequences.with.Motif))
  f
})

top10_p = do.call(rbind, top10_d)
colnames(top10_p) = c("p", "percentage")
top10_p$group = rep(c("dependent", "independent"), each = length(top10))
top10_p$Motif = stringr::str_split_fixed(top10, "/", 3)[,1]
top10_p$p[top10_p$p < 1e-200] <- 1e-200
ggplot(top10_p, aes(group, Motif, size = percentage, color = -log10(p)))+
  geom_point()+
  cowplot::theme_cowplot()+
  scale_color_gradientn(colours = c("blue4", "gray90", "yellow", "red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
ggsave("../figure//Brg1_independent_Motif_enrich.pdf", width = 5, height = 5)

#
top10_p = do.call(cbind, top10_d)
top10_p$diff = abs(log(top10_p[,2] / top10_p[,4]))
top10_p = top10_p[order(top10_p$diff, decreasing = T), ]
top10_p = rownames(top10_p[1:10,])

top10_d = lapply(list.files("/nfs4/chaozhang/proj/Embryo/Yota/BRG1_ATAC_0901/script/motif_dependent", "knownResults.txt", recursive =T, full.names = T), function(i){
  f = read.table(i,sep="\t", stringsAsFactors = F, comment.char = "*", header = T)
  f = f[!duplicated(f$Motif.Name), ]
  rownames(f) = f$Motif.Name
  f = f[top10_p, c("P.value", "X..of.Target.Sequences.with.Motif")]
  f$X..of.Target.Sequences.with.Motif = as.numeric(sub("%", "", f$X..of.Target.Sequences.with.Motif))
  f
})

top10_p = do.call(rbind, top10_d)
colnames(top10_p) = c("p", "percentage")
top10_p$group = rep(c("dependent", "independent"), each = 10)
top10_p$Motif = stringr::str_split_fixed(rownames(top10_p), "/", 3)[,1]
top10_p$p[top10_p$p==0] <- 1e-300
ggplot(top10_p, aes(group, Motif, size = percentage, color = -log10(p)))+
  geom_point()+
  cowplot::theme_cowplot()+
  scale_color_gradientn(colours = c("blue4", "blue", "gray90", "yellow", "red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
ggsave("../figure//Brg1_independent_Motif_enrich.pdf", width = 5, height = 5)


#motif promoter distal
top10 = lapply(list.files("/nfs4/chaozhang/proj/Embryo/Yota/BRG1_ATAC_0901/script/motif_dependent_pro/", "knownResults.txt", recursive =T, full.names = T), function(i){
  f = read.table(i,sep="\t", stringsAsFactors = F, comment.char = "*", header = T)
  f$Motif.Name[1:10]
})
top10 = unique(unlist(top10))

top10_d = lapply(list.files("/nfs4/chaozhang/proj/Embryo/Yota/BRG1_ATAC_0901/script/motif_dependent_pro/", "knownResults.txt", recursive =T, full.names = T), function(i){
  f = read.table(i,sep="\t", stringsAsFactors = F, comment.char = "*", header = T)
  f = f[!duplicated(f$Motif.Name), ]
  rownames(f) = f$Motif.Name
  f = f[top10, c("P.value", "X..of.Target.Sequences.with.Motif")]
  f$X..of.Target.Sequences.with.Motif = as.numeric(sub("%", "", f$X..of.Target.Sequences.with.Motif))
  f
})

top10_p = do.call(rbind, top10_d)
colnames(top10_p) = c("p", "percentage")
top10_p$group = rep(c("dependent_distal", "dependent_promoter", "independent_distal", "independent_promoter"), each = length(top10))
top10_p$Motif = stringr::str_split_fixed(top10, "/", 3)[,1]
top10_p$p[top10_p$p < 1e-200] <- 1e-200
ggplot(top10_p, aes(group, Motif, size = percentage, color = -log10(p)))+
  geom_point()+
  cowplot::theme_cowplot()+
  scale_color_gradientn(colours = c("blue4", "gray90", "yellow", "red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
ggsave("../figure//Brg1_independent_pro_Motif_enrich.pdf", width = 5, height = 5)

top10_d = lapply(list.files("/nfs4/chaozhang/proj/Embryo/Yota/BRG1_ATAC_0901/script/motif_dependent_pro/", "knownResults.txt", recursive =T, full.names = T), function(i){
  f = read.table(i,sep="\t", stringsAsFactors = F, comment.char = "*", header = T)
  f = f[!duplicated(f$Motif.Name), ]
  rownames(f) = f$Motif.Name
  f = f[top10_p, c("P.value", "X..of.Target.Sequences.with.Motif")]
  f$X..of.Target.Sequences.with.Motif = as.numeric(sub("%", "", f$X..of.Target.Sequences.with.Motif))
  f
})

plotMotifPie = function(f, m){
  pct = sapply(list.files(f, "knownResults.txt", recursive =T, full.names = T), function(i){
    f = read.table(i,sep="\t", stringsAsFactors = F, comment.char = "*", header = T)
    f$Motif = stringr::str_split_fixed(f$Motif.Name, "/", 3)[,1]
    f$X..of.Target.Sequences.with.Motif = as.numeric(sub("%", "", f$X..of.Target.Sequences.with.Motif))
    f[f$Motif == m, "X..of.Target.Sequences.with.Motif"]
  })
  pct = as.data.frame(pct)
  pct$other = 100-pct$pct
  
  par(mfrow=c(1, nrow(pct)))
  plot_list = lapply(1:nrow(pct), function(i){
    pie(as.numeric(pct[i,]), labels = c("with motif", ""), col = c("pink", "gray90" ), border = "white")
  })
}

plotMotifPie("/nfs4/chaozhang/proj/Embryo/Yota/BRG1_ATAC_0901/script/motif_dependent_pro/", "CTCF(Zf)")
plotMotifPie("/nfs4/chaozhang/proj/Embryo/Yota/BRG1_ATAC_0901/script/motif_dependent_pro/", "Oct4(POU,Homeobox)")
plotMotifPie("/nfs4/chaozhang/proj/Embryo/Yota/BRG1_ATAC_0901/script/motif_dependent_pro/", "NFY(CCAAT)")

plotMotifPie("/nfs4/chaozhang/proj/Embryo/Yota/BRG1_ATAC_0901/script/motif_dependent/", "CTCF(Zf)")
plotMotifPie("/nfs4/chaozhang/proj/Embryo/Yota/BRG1_ATAC_0901/script/motif_dependent/", "NFY(CCAAT)")
plotMotifPie("/nfs4/chaozhang/proj/Embryo/Yota/BRG1_ATAC_0901/script/motif_dependent/", "Oct4(POU,Homeobox)")
plotMotifPie("/nfs4/chaozhang/proj/Embryo/Yota/BRG1_ATAC_0901/script/motif_dependent/", "Brn1(POU,Homeobox)")

#
files = list(Down_0.3n = "./Brg1_DEGs_0.3n_Ctrl.Down.bed",
              Down_1n = "./Brg1_DEGs_1n_Ctrl.Down.bed",
             Down_3n = "./Brg1_DEGs_3n_Ctrl.Down.bed",
             Down_10n = "./Brg1_DEGs_10n_Ctrl.Down.bed",
             Down_100n = "./Brg1_DEGs_100n_Ctrl.Down.bed",
             All_peak = "/nfs4/chaozhang/proj/Embryo/Yota/BRG1_CutRun_Nov2023/peaks/Ctr/Ctr_peaks.narrowPeak")

peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb,
                       annoDb="org.Mm.eg.db",
                       tssRegion=c(-1000, 1000), verbose=FALSE)

plotAnnoBar(peakAnnoList)
ggsave("../figure/DAR_Down_peaks_Features.pdf", width = 6, height = 5)



plotAnnoPie(peakAnnoList[[2]])
plotAnnoPie(peakAnnoList[[3]])
plotAnnoPie(peakAnnoList[[4]])
plotAnnoPie(peakAnnoList[[5]])
plotAnnoPie(peakAnnoList[[6]])



#

peak_ann = annotatePeak("../peaks/Brg_new.ATAC_reads.bed",
             annoDb="org.Mm.eg.db",
             TxDb=txdb,
             tssRegion=c(-1000, 1000))
plotAnnoPie(peak_ann)
peak_ann = as.data.frame(peak_ann)
write.table(peak_ann, "Brg1_peak_ann.bed", quote = F, col.names = F, row.names = F, sep = "\t")
peak_ann$ann = stringr::str_replace_all(peak_ann$annotation, " .*", "")
peak_ann_promoter = peak_ann[peak_ann$annotation == "Promoter", ]
peak_ann_promoter$start = peak_ann_promoter$start - 1
peak_ann_promoter = paste(peak_ann_promoter$seqnames, peak_ann_promoter$start, peak_ann_promoter$end)
#
DAR_100n = read.table("Brg1_DEGs_100n_Ctrl.txt", sep = "\t")
write.table(stringr::str_split_fixed(rownames(DAR_100n[DAR_100n$sig == "Down" & rownames(DAR_100n) %in%peak_ann_promoter, ]), " ",3),
            "Brg1_DEGs_100n_Ctrl_Down_Promoter.bed", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(stringr::str_split_fixed(rownames(DAR_100n[DAR_100n$sig == "Down" & !rownames(DAR_100n) %in%peak_ann_promoter, ]), " ",3),
            "Brg1_DEGs_100n_Ctrl_Down_Distal.bed", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(stringr::str_split_fixed(rownames(DAR_100n[DAR_100n$sig == "none" & rownames(DAR_100n) %in%peak_ann_promoter, ]), " ",3),
            "Brg1_DEGs_100n_Ctrl_none_Promoter.bed", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(stringr::str_split_fixed(rownames(DAR_100n[DAR_100n$sig == "none" & !rownames(DAR_100n) %in%peak_ann_promoter, ]), " ",3),
            "Brg1_DEGs_100n_Ctrl_none_Distal.bed", sep = "\t", quote = F, row.names = F, col.names = F)


# Promoter
DAR_promoter = function(c1, c2){
  Brg1_DEGs = read.table(paste0("Brg1_DEGs_", c1, "_", c2, ".txt"),sep = "\t")
  Brg1_DEGs = Brg1_DEGs[rownames(c) %in% peak_ann_promoter, ]
  print(table(Brg1_DEGs$sig))
  
  y1_normalized_counts_group1 = y1_normalized_counts_group[rownames(y1_normalized_counts_group)%in% peak_ann_promoter, ]
  y1_normalized_counts_group1 = as.data.frame(y1_normalized_counts_group1)
  
  ggplot(y1_normalized_counts_group1, aes(get(c2), get(c1))) +
    geom_point(size = .1, alpha = .3) +
    geom_point(data = (y1_normalized_counts_group1[rownames(Brg1_DEGs)[Brg1_DEGs$sig == "Down"],]), size = .5, color = "blue") +
    geom_point(data = (y1_normalized_counts_group1[rownames(Brg1_DEGs)[Brg1_DEGs$sig == "Up"],]), size = .5, color = "red") +
    geom_abline(slope = 1) + xlim(1,10) + ylim(1,10) + xlab(c2) + ylab(c1)+
    scale_fill_viridis_c() + theme_classic(base_size = 15) 
  
}

p1 = DAR_promoter("0.3n", "Ctrl")
p2 = DAR_promoter("1n", "Ctrl")
p3 = DAR_promoter("3n", "Ctrl")
p4 = DAR_promoter("10n", "Ctrl")
p5 = DAR_promoter("100n", "Ctrl")
p1|p2|p3|p4|p5
ggsave("../figure/ATAC_at_BRG1_peaks.Promoter_DAR.pair-wise_scatterplot.png", width = 25, height = 5)

# enhancer
DAR_enhancer = function(c1, c2){
  
  Brg1_DEGs = read.table(paste0("Brg1_DEGs_", c1, "_", c2, ".txt"),sep = "\t")
  Brg1_DEGs = Brg1_DEGs[!rownames(Brg1_DEGs) %in% peak_ann_promoter, ]
  print(table(Brg1_DEGs$sig))
  
  y1_normalized_counts_group1 = y1_normalized_counts_group[!rownames(y1_normalized_counts_group)%in% peak_ann_promoter, ]
  y1_normalized_counts_group1 = as.data.frame(y1_normalized_counts_group1)
  
  ggplot(y1_normalized_counts_group1, aes(get(c2), get(c1))) +
    geom_point(size = .1, alpha = .3) +
    geom_point(data = (y1_normalized_counts_group1[rownames(Brg1_DEGs)[Brg1_DEGs$sig == "Down"],]), size = .5, color = "blue") +
    geom_point(data = (y1_normalized_counts_group1[rownames(Brg1_DEGs)[Brg1_DEGs$sig == "Up"],]), size = .5, color = "red") +
    geom_abline(slope = 1) + xlim(1,10) + ylim(1,10) + xlab(c2) + ylab(c1)+
    scale_fill_viridis_c() + theme_classic(base_size = 15) 
  
}

p1 = DAR_enhancer("0.3n", "Ctrl")
p2 = DAR_enhancer("1n", "Ctrl")
p3 = DAR_enhancer("3n", "Ctrl")
p4 = DAR_enhancer("10n", "Ctrl")
p5 = DAR_enhancer("100n", "Ctrl")
p1|p2|p3|p4|p5
ggsave("../figure/ATAC_at_BRG1_peaks.Enhancer_DAR.pair-wise_scatterplot.png", width = 25, height = 5)


#
DEGs_100n = read.table(paste0("Brg1_DEGs_100n_Ctrl.txt"),sep = "\t")
dim(DEGs_100n[DEGs_100n$sig == "Down" & rownames(DEGs_100n) %in% peak_ann_promoter, ])
dim(DEGs_100n[DEGs_100n$sig == "Down" & !rownames(DEGs_100n) %in% peak_ann_promoter, ])
write.table(stringr::str_replace_all(rownames(DEGs_100n[DEGs_100n$sig == "Down" & rownames(DEGs_100n) %in% peak_ann_promoter, ]), " ", "\t"),
            "Brg1_peak_ann.Promoter.Down.bed", quote = F, row.names = F, col.names = F)
write.table(stringr::str_replace_all(rownames(DEGs_100n[DEGs_100n$sig == "none" & rownames(DEGs_100n) %in% peak_ann_promoter, ]), " ", "\t"),
            "Brg1_peak_ann.Promoter.none.bed", quote = F, row.names = F, col.names = F)
write.table(stringr::str_replace_all(rownames(DEGs_100n[DEGs_100n$sig == "Down" & !rownames(DEGs_100n) %in% peak_ann_promoter, ]), " ", "\t"),
            "Brg1_peak_ann.Enhancer.Down.bed", quote = F, row.names = F, col.names = F)
write.table(stringr::str_replace_all(rownames(DEGs_100n[DEGs_100n$sig == "none" & !rownames(DEGs_100n) %in% peak_ann_promoter, ]), " ", "\t"),
            "Brg1_peak_ann.Enhancer.none.bed", quote = F, row.names = F, col.names = F)


pheatmap::pheatmap(y1_normalized_counts_group[rownames(y1_normalized_counts_group)%in% peak_ann_promoter, c("WT", "Ctrl", "0.3n", "1n", "3n", "10n", "100n")],
                   show_rownames = F, scale = "none",
                   cluster_cols = F,
                   clustering_method = "ward.D",
                   color = colorRampPalette(c("blue", "white", "red"))(length(breaksList)),
                   breaks = breaksList)


breaksList = seq(0, 4, by = .1)
pheatmap::pheatmap(log(y1_normalized_counts_group+1)[rownames(y1_normalized_counts_group)%in% peak_ann_promoter, c("WT", "Ctrl", "0.3n", "1n", "3n", "10n", "100n")], 
                   scale = "row", show_rownames = F,
                   border_color = NA,
                   cluster_cols = F, cluster_rows = T,
                   color = pals::gnuplot(length(breaksList)),
                   breaks = breaksList)
dev.off()


