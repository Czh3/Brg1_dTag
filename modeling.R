# modeling
# sensitive & buffered
setwd("/nfs4/chaozhang/proj/Embryo/Yota/BRG1_ATAC_0901/script")
library(qlcMatrix)

normalized_counts = read.csv("ATAC_at_brg1_normalized_counts_unlog.csv", header = T, row.names = 1)
DAR = read.table("DEG_dTag100n_Ctrl.txt", header = T, row.names = 1, sep = "\t")
DAR_down = DAR[DAR$sig == "Down", ]

Down_normalized_counts = normalized_counts[rownames(DAR_down), ]

Brg1_level = c(.78, .78, .03, .03, .09, .09, .42, .42, .22, .22, .61, .61, 1, 1)

Down_normalized_counts = Down_normalized_counts/rowMeans(Down_normalized_counts[,c("WT_rep1", "WT_rep2")])
Zscores = apply(Down_normalized_counts, 1, function(x) (x-mean(x))/sd(x))
summary(rowMax(t(Zscores)))
Down_normalized_counts = Down_normalized_counts[rowMax(t(Zscores)) <= 2.5, ]
DAR_down = DAR_down[rowMax(t(Zscores)) <= 2.5, ]

res = apply(Down_normalized_counts, 1, function(i){
  dis1(1-Brg1_level, i) - dis2(1-Brg1_level, i)
})

DAR_down$delta = res
DAR_down$group = "Group_1"
DAR_down$group[DAR_down$delta > -0.5 & DAR_down$delta < 0] = "Group_2"
DAR_down$group[DAR_down$delta >= 0] = "Group_3"

write.csv(DAR_down, "DAR_down_delta.csv", quote = F)

ggplot(DAR_down, aes(logCPM, delta))+
  geom_point()

ggplot(DAR_down, aes(x=delta)) + 
  geom_histogram(aes(y=..count..), colour="black", fill="gray80")+
  #geom_density(alpha=.2, fill="blue") +
  geom_vline(xintercept = 0, color = "red3") +
  cowplot::theme_cowplot() 

ggsave("../figure/Model_bias.barplot.pdf", width = 5, height = 5)

library(fitdistrplus)
descdist(DAR_down$delta, discrete = FALSE)


pheatmap::pheatmap(Down_normalized_counts[rownames(DAR_down[DAR_down$delta > 0, ]), c(13,14,1,2,11,12,7,8,9,10,5,6,3,4)], 
                   scale = "row", show_rownames = F,
                   cluster_cols = F, cluster_rows = T)

breaksList = seq(0, 2, by = .1)
png("../figure/DAR_model_3group.heatmap.png", 2500, 3500, res = 500)
pheatmap::pheatmap((Down_normalized_counts)[rownames(DAR_down[order(DAR_down$delta),]), c(13,14,1,2,11,12,7,8,9,10,5,6,3,4)], 
                   scale = "none", show_rownames = F,
                   border_color = NA,
                   cluster_cols = F, cluster_rows = F,
                   gaps_row = c(sum(DAR_down$delta < -0.5), nrow(DAR_down[DAR_down$delta < 0, ])),
                   color = pals::gnuplot(length(breaksList)),
                   breaks = breaksList)
dev.off()


normalized_counts_group = read.csv("ATAC_normalized_counts_group.csv", header = T, row.names = 1)
normalized_counts_group = normalized_counts_group[rownames(DAR_down), ]
normalized_counts_group$group = DAR_down$group
normalized_counts_group = melt(normalized_counts_group)
normalized_counts_group$variable = factor(normalized_counts_group$variable,
                                          levels = c("WT", "Ctrl", "X0.3n", "X1n", "X3n", "X10n", "X100n"))
p1 = ggplot(normalized_counts_group, aes(variable, value, fill = group)) +
  geom_boxplot(outlier.size = -1) + ylim(0,20) +
  cowplot::theme_cowplot() + xlab("dTag KD") + ylab("ATAC signals") +
  scale_fill_manual(values = pals::gnuplot(10)[c(9,6,3)])

p2 = ggplot(normalized_counts_group, aes(variable, value, fill = group)) +
  geom_point(color = "white") + ylim(0,15) +
  cowplot::theme_cowplot() + xlab("dTag KD") + ylab("ATAC signals") + 
  geom_smooth(aes(as.numeric(variable), value, color = group), method = "loess") + #ylim(0,10) +
  scale_color_manual(values = pals::gnuplot(10)[c(9,6,3)])

p1/p2
ggsave("../figure/Model_3Group_ATAC_changes.barplot.pdf", width = 5, height = 8)

normalized_counts_group = read.csv("ATAC_normalized_counts_group.csv", header = T, row.names = 1)
normalized_counts_group = normalized_counts_group[rownames(DAR_down), ]
colnames(normalized_counts_group) = c(61,3,9,42,22,78,100)
normalized_counts_group$group = DAR_down$group
normalized_counts_group = melt(normalized_counts_group)
normalized_counts_group$variable = as.numeric(as.character(normalized_counts_group$variable))

ggplot(normalized_counts_group, aes(variable, value, fill = group)) +
  cowplot::theme_cowplot() + xlab("Brg1 Expr (%)") + ylab("ATAC signals") + 
  geom_smooth(aes(as.numeric(variable), value, color = group), method = "loess", fill = "gray80") + #ylim(0,10) +
  scale_color_manual(values = pals::gnuplot(10)[c(9,6,3)]) +
  scale_x_reverse()
ggsave("../figure/Model_3Group_ATAC_changes.linePlot.pdf", width = 5, height = 4)


## peak density
normalized_counts_group = read.csv("ATAC_normalized_counts_group.csv", header = T, row.names = 1)
normalized_counts_group = normalized_counts_group[rownames(DAR_down), "Ctrl", drop = F]
normalized_counts_group$group = DAR_down$group
normalized_counts_group = melt(normalized_counts_group)

ggplot(normalized_counts_group, aes(group, log(value), fill = group)) +
  cowplot::theme_cowplot() + xlab("") + ylab("log(ATAC signals)") + 
  geom_boxplot( outlier.shape = NA) +
  scale_fill_manual(values = pals::gnuplot(10)[c(9,6,3)])
ggsave("../figure/Model_3Group_PeakDensity.pdf", width = 5, height = 5)


## GC content
DAR_GC = read.table("DAR_down_delta.GCcontent.bed", header = T, comment.char = "?")
rownames(DAR_GC) = paste(DAR_GC$X.1_usercol, DAR_GC$X2_usercol, DAR_GC$X3_usercol)

DAR_GC$group = DAR_down$group
DAR_GC$delta = DAR_down$delta

ggplot(DAR_GC, aes(x = X5_pct_gc, fill = group)) +
  geom_density(alpha = .3) +
  scale_fill_manual(values = pals::gnuplot(10)[c(9,6,3)]) +
  cowplot::theme_cowplot() + xlab("GC content")

ggplot(DAR_GC, aes(delta, X5_pct_gc, color = group)) +
  geom_point(alpha = .3) +
  cowplot::theme_cowplot() + xlab("GC content")

t.test(DAR_GC[DAR_GC$group == "Group_3", "X5_pct_gc"], DAR_GC[DAR_GC$group == "Group_1", "X5_pct_gc"])

ggplot(DAR_GC, aes(group, X12_seq_len, fill = group)) +
  geom_boxplot(outlier.shape = NA ) +
  cowplot::theme_cowplot() +
  scale_fill_manual(values = pals::gnuplot(10)[c(9,6,3)]) + ylim(100, 1500) + ylab("Peaks length") + xlab("")
ggsave("../figure/Model_3Group_PeakLength.pdf", width = 5, height = 5)

t.test(DAR_GC[DAR_GC$group == "Group_3", "X12_seq_len"], DAR_GC[DAR_GC$group == "Group_2", "X12_seq_len"]) #p-value < 2.2e-16
t.test(DAR_GC[DAR_GC$group == "Group_1", "X12_seq_len"], DAR_GC[DAR_GC$group == "Group_2", "X12_seq_len"]) #p-value < 2.287e-12

ggplot(DAR_GC, aes(delta, X12_seq_len, color = group)) +
  geom_point() +
  cowplot::theme_cowplot() 


library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

files = list(G1 = "./DAR_down_delta.Group_1.bed",
             G2 = "./DAR_down_delta.Group_2.bed",
             G3 = "./DAR_down_delta.Group_3.bed")

peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb,
                       annoDb="org.Mm.eg.db",
                       tssRegion=c(-1000, 1000), verbose=FALSE)

plotAnnoBar(peakAnnoList)

genes= lapply(peakAnnoList, function(i) as.data.frame(i)$SYMBOL)
vennplot(genes)
library(clusterProfiler)
compGO <- compareCluster(geneCluster   = genes,
                           fun           = "enrichGO",
                           keyType = "SYMBOL",
                         OrgDb = org.Mm.eg.db,
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "BH")
dotplot(compGO, showCategory = 10, title = "GO Enrichment Analysis")




###### Only Brg1 peaks

normalized_counts = read.csv("ATAC_at_brg1_normalized_counts_unlog.csv", header = T, row.names = 1)
DAR = read.table("Brg1_DEGs_100n_Ctrl.txt", header = T, row.names = 1, sep = "\t")
DAR_down = DAR[DAR$sig == "Down", ]


Down_normalized_counts = normalized_counts[rownames(DAR_down), ]
#Brg1_level = c(.78, .78, .03, .03, .09, .09, .42, .42, .22, .22, .61, .61, 1, 1)

# remove WT
Down_normalized_counts = Down_normalized_counts[, !colnames(Down_normalized_counts) %in% c("WT_rep1", "WT_rep2")]
Brg1_level = c(.78, .78, .03, .03, .09, .09, .42, .42, .22, .22, .61, .61)
Brg1_level = Brg1_level/0.78
  
Down_normalized_counts = Down_normalized_counts/(rowMeans(Down_normalized_counts[,c("Ctrl_rep1", "Ctrl_rep2")]+0.01))
Zscores = apply(Down_normalized_counts, 1, function(x) {
  x = as.numeric(x)
  (x-mean(x))/sd(x)
  })

Zscores = t(Zscores)
summary(rowMax(Zscores))
Down_normalized_counts = Down_normalized_counts[rowMax(Zscores) <= 3, ]
DAR_down = DAR_down[rowMax(Zscores) <= 3, ]

res = apply(Down_normalized_counts, 1, function(i){
  dis1(1-Brg1_level, i) - dis2(1-Brg1_level, i)
})

DAR_down$delta = res
quan = quantile(DAR_down$delta, probs = seq(0,1,0.2))

for(i in 1:(length(quan)-1)){
  DAR_down$group[DAR_down$delta <= quan[i+1] & DAR_down$delta >= quan[i]] <- paste0("Group_", i)
}

#DAR_down$group = "Group_1"
#DAR_down$group[DAR_down$delta > -0.5 & DAR_down$delta < 0] = "Group_2"
#DAR_down$group[DAR_down$delta >= 0] = "Group_3"

write.csv(DAR_down, "Brg1_DAR_down_delta.csv", quote = F)


ggplot(DAR_down, aes(logCPM, delta))+
  geom_point() +
  geom_smooth(method = "loess", fill = "gray80")
  

ggplot(DAR_down, aes(x=delta)) + 
  geom_histogram(aes(y=..count..), colour="black",  fill = "gray80")+
  #geom_density(alpha=.2, fill="blue") +
  geom_vline(xintercept = 0, color = "red3") +
  cowplot::theme_cowplot() 

ggsave("../figure/Brg1_Model_bias.barplot.pdf", width = 5, height = 5)
dim(DAR_down[DAR_down$delta < 0,])

### only focuse on Enhancer
DAR_enhancer = read.table("Brg1_peak_ann.Enhancer.bed")
DAR_enhancer$name = paste(DAR_enhancer$V1, DAR_enhancer$V2, DAR_enhancer$V3)
DAR_down = DAR_down[rownames(DAR_down) %in% DAR_enhancer$name, ]
##


breaksList = seq(0, 1.5, by = .1)
pdf("../figure/Brg1_DAR_model_Enhancer_5group.heatmap.pdf", 5, 7)
pheatmap::pheatmap((Down_normalized_counts)[rownames(DAR_down[order(DAR_down$delta),]), c(1,2,11,12,7,8,9,10,5,6,3,4)], 
                   scale = "none", show_rownames = F,
                   border_color = NA,
                   cluster_cols = F, cluster_rows = F,
                   gaps_row = seq(1, nrow(DAR_down), length.out = 6)[-1],
                   color = pals::gnuplot(length(breaksList)),
                   breaks = breaksList)
dev.off()



normalized_counts1 = normalized_counts[rownames(DAR_down), ]
pheatmap::pheatmap(log(normalized_counts1+1)[rownames(DAR_down[order(DAR_down$delta),]), c(13,14,1,2,11,12,7,8,9,10,5,6,3,4)], 
                   scale = "none", show_rownames = F,
                   border_color = NA,
                   cluster_cols = F, cluster_rows = F,
                   gaps_row = c(sum(DAR_down$delta < -0.5), nrow(DAR_down[DAR_down$delta < 0, ])),
                   color = pals::jet(50))



delta_density = as.data.frame(DAR_down[order(DAR_down$delta),])
delta_density$rank = 1:nrow(delta_density)
ggplot(delta_density, aes(rank, delta)) +
  geom_bar(stat = "identity", fill = ifelse(delta_density$delta > 0, "red", "blue")) +
  coord_flip() +
  theme_classic()



normalized_counts1 = normalized_counts

normalized_counts_group = normalized_counts1[, seq(1,14,2)] + normalized_counts1[, seq(2,14,2)] 
normalized_counts_group = normalized_counts_group/2
colnames(normalized_counts_group) = stringr::str_remove(colnames(normalized_counts_group), "_rep[1|2]")
normalized_counts_group$group = DAR_down[rownames(normalized_counts_group), ]$group

#normalized_counts_group$group[is.na(normalized_counts_group$group)] <- "Constant"
normalized_counts_group = na.omit(normalized_counts_group)
normalized_counts_group1 = normalized_counts_group %>% group_by(group) %>% summarise_all("median")
normalized_counts_group1 = as.data.frame(normalized_counts_group1)
rownames(normalized_counts_group1) = normalized_counts_group1$group
normalized_counts_group1 = normalized_counts_group1[,-1]
normalized_counts_group1 = normalized_counts_group1/normalized_counts_group1$Ctrl

pdf("../figure/Brg1_DAR_model_Enhancer_5group_mean.heatmap.pdf")
pheatmap::pheatmap(log(normalized_counts_group1[,c(1,6,4,5,3,2)]+1), scale = "none", cluster_cols = F, cluster_rows = F,
                   border_color = NA,  color = pals::gnuplot(20))
dev.off()

normalized_counts_group = melt(normalized_counts_group)
normalized_counts_group$variable = factor(normalized_counts_group$variable,
                                          levels = c("WT", "Ctrl", "X0.3n", "X1n", "X3n", "X10n", "X100n"))


ggplot(normalized_counts_group, aes(variable, log(value+1), fill = group)) +
  geom_boxplot(outlier.size = -1) + 
  cowplot::theme_cowplot() + xlab("dTag KD") + ylab("log(ATAC signals)") +
  scale_fill_manual(values = pals::gnuplot(9)[4:9])
#ggsave("../figure/Brg1_Model_5Group_ATAC_changes.barplot.pdf", width = 5, height = 4)


ggplot(normalized_counts_group, aes(group, log(value+1), fill = variable)) +
  geom_boxplot(outlier.size = -1) + 
  cowplot::theme_cowplot() + xlab("") + ylab("log(ATAC signals)") + ylim(0,8)+
  scale_fill_manual(values = pals::gnuplot(9)[3:9])
#ggsave("../figure/Brg1_Model_5Group_ATAC_changes1.barplot.pdf", width = 5, height = 4)

ggplot(normalized_counts_group, aes(group, value, fill = variable)) +
  geom_boxplot(outlier.size = -1) + 
  cowplot::theme_cowplot() + xlab("") + ylab("ATAC signals") + ylim(0,8)+
  scale_fill_manual(values = pals::gnuplot(9)[3:9])+ ylim(0,100)
#ggsave("../figure/Brg1_Model_5Group_ATAC_changes1_unlog.barplot.pdf", width = 5, height = 4)



#
normalized_counts1 = normalized_counts[, 1:12]
normalized_counts_group = normalized_counts1[, seq(1,12,2)]/normalized_counts1$Ctrl_rep1 + normalized_counts1[, seq(2,12,2)]/normalized_counts1$Ctrl_rep2
normalized_counts_group = normalized_counts_group/2

colnames(normalized_counts_group) = Brg1_level[seq(1,12,2)]
normalized_counts_group$group = DAR_down[rownames(normalized_counts_group), ]$group
#normalized_counts_group$group[is.na(normalized_counts_group$group)] <- "Constant"
normalized_counts_group = na.omit(normalized_counts_group)
normalized_counts_group = melt(normalized_counts_group)
normalized_counts_group$variable = as.numeric(as.character(normalized_counts_group$variable))

library(locfit)
ggplot(normalized_counts_group, aes(variable, value, fill = group)) +
  geom_smooth(aes(variable, value, color = group), method = "locfit",method.args = list(deg=2, alpha=1), fill = "gray80") + #ylim(0,10) +
  cowplot::theme_cowplot() + xlab("Brg1 Expr (%)") + ylab("Normalized ATAC signals") + 
  scale_color_manual(values = c( pals::gnuplot(9)[4:9])) +
  scale_x_reverse()
ggsave("../figure/Brg1_Model_5Group_Enhancer_ATAC_changes.linePlot.pdf", width = 5, height = 4)

ggplot(normalized_counts_group[sample(1:nrow(normalized_counts_group), 20000),], aes(variable, value, fill = group)) +
  cowplot::theme_cowplot() + xlab("Brg1 Expr (%)") + ylab("ATAC signals") + 
  geom_smooth(aes(variable, value, color = group), method = "loess", fill = "gray80") + #ylim(0,10) +
  scale_color_manual(values = pals::gnuplot(9)[4:9]) +
  scale_x_reverse()



## peak density
normalized_counts_group = read.csv("ATAC_at_brg1_normalized_counts_group_unlog.csv", header = T, row.names = 1)
normalized_counts_group = normalized_counts_group[rownames(DAR_down), "Ctrl", drop = F]
normalized_counts_group$group = DAR_down$group
normalized_counts_group = melt(normalized_counts_group)

ggplot(normalized_counts_group, aes(group, log(value), fill = group)) +
  cowplot::theme_cowplot() + xlab("") + ylab("log(ATAC signals)") + 
  geom_boxplot( outlier.shape = NA) +
  scale_fill_manual(values = pals::gnuplot(9)[4:9])
ggsave("../figure/Brg1_Model_5Group_PeakDensity.pdf", width = 5, height = 5)


## GC content
# bedtools nuc -fi /nfs4/chaozhang/database/mm10/Sequence/WholeGenomeFasta/genome.fa -bed Brg1_DEGs_0.3n_Ctrl.bed > Brg1_peaks_GCcontent.bed
DAR_GC = read.table("Brg1_peaks_GCcontent.bed", header = T, comment.char = "?")
rownames(DAR_GC) = paste(DAR_GC$X.1_usercol, DAR_GC$X2_usercol, DAR_GC$X3_usercol)
DAR_GC = DAR_GC[rownames(DAR_down), ]

DAR_GC$group = DAR_down$group
DAR_GC$delta = DAR_down$delta

ggplot(DAR_GC, aes(x = X5_pct_gc, color = group)) +
  geom_density(size = 1) +
  scale_color_manual(values = pals::gnuplot(9)[4:9]) +
  cowplot::theme_cowplot() + xlab("GC content")
ggsave("../figure/Brg1_Model_5Group_Enhaancer_GC_content.pdf", width = 5, height = 5)

ggplot(DAR_GC, aes(group, X5_pct_gc, fill = group)) +
  geom_boxplot() +
  scale_fill_manual(values = pals::gnuplot(9)[4:9]) +
  cowplot::theme_cowplot() + ylab("GC content")
ggsave("../figure/Brg1_Model_5Group_Enhaancer_GC_content_boxplot.pdf", width = 5, height = 5)

wilcox.test(DAR_GC[DAR_GC$group == "Group_1", "X5_pct_gc"], DAR_GC[DAR_GC$group == "Group_2", "X5_pct_gc"]) # 0.0002542

ggplot(DAR_GC, aes(delta, X5_pct_gc, color = group)) +
  geom_point(alpha = .3) +
  cowplot::theme_cowplot() + xlab("delta")


ggplot(DAR_GC, aes(group, X12_seq_len, fill = group)) +
  geom_boxplot(outlier.shape = NA ) +
  cowplot::theme_cowplot() +
  scale_fill_manual(values = pals::gnuplot(9)[4:9]) + ylim(100, 1000) + ylab("Peaks length") + xlab("")
ggsave("../figure/Brg1_Model_5Group_Enhancer_PeakLength.pdf", width = 5, height = 5)

wilcox.test(DAR_GC[DAR_GC$group == "Group_3", "X12_seq_len"], DAR_GC[DAR_GC$group == "Group_2", "X12_seq_len"]) #p-value = 9.72e-05
wilcox.test(DAR_GC[DAR_GC$group == "Group_1", "X12_seq_len"], DAR_GC[DAR_GC$group == "Group_2", "X12_seq_len"]) #p-value = 0.5877

ggplot(DAR_GC, aes(delta, X12_seq_len, color = group)) +
  geom_point() +
  cowplot::theme_cowplot() 


#
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

lapply(paste0("Group_", 1:5), function(i){
  g = DAR_down[DAR_down$group == i,]
  p = rownames(g)
  p = as.data.frame(stringr::str_split_fixed(p, " ", 3))
  p$g = g$group
  write.table(p, paste0("./Brg1_DAR_down_Enhancer_delta.",i,".bed"), sep = "\t", col.names = F, row.names = F,quote = F)
})

files = list(G1 = "./Brg1_DAR_down_delta.Group_1.bed",
             G2 = "./Brg1_DAR_down_delta.Group_2.bed",
             G3 = "./Brg1_DAR_down_delta.Group_3.bed",
             G4 = "./Brg1_DAR_down_delta.Group_4.bed",
             G5 = "./Brg1_DAR_down_delta.Group_5.bed")

peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb,
                       annoDb="org.Mm.eg.db",
                       tssRegion=c(-1000, 1000), verbose=FALSE)

plotAnnoBar(peakAnnoList)
ggsave("../figure/Brg1_Model_5Group_Annotation.pdf", width = 7, height = 4)


genes= lapply(peakAnnoList, function(i) {
  i = as.data.frame(i)
  i = i[i$annotation == "Promoter", ]
  i$SYMBOL
})
vennplot(genes)
library(clusterProfiler)
compGO <- compareCluster(geneCluster   = genes,
                         fun           = "enrichGO",
                         keyType = "SYMBOL",
                         OrgDb = org.Mm.eg.db,
                         pvalueCutoff  = 0.5,
                         pAdjustMethod = "BH")
dotplot(compGO, showCategory = 10, title = "GO Enrichment Analysis")

# promoter
peak_ann = annotatePeak("/nfs4/chaozhang/proj/Embryo/Yota/BRG1_CutRun_Apr2024/peaks/Brg1_Ctrl_peaks.bed",
                        annoDb="org.Mm.eg.db",
                        TxDb=txdb,
                        tssRegion=c(-1000, 1000))
peak_ann = as.data.frame(peak_ann)
peak_ann$ann = stringr::str_replace_all(peak_ann$annotation, " .*", "")
peak_ann_promoter = peak_ann[peak_ann$annotation == "Promoter", ]
peak_ann_promoter$start = peak_ann_promoter$start - 1
peak_ann_promoter = paste(peak_ann_promoter$seqnames, peak_ann_promoter$start, peak_ann_promoter$end)

for (i in files){
  a = read.table(i)
  a.p = a[paste(a$V1, a$V2, a$V3) %in% peak_ann_promoter, ]
  write.table(a.p, stringr::str_replace(i, ".bed", ".Promoter.bed"), quote = F, sep = "\t", col.names = F, row.names = F)
  a.e = a[!paste(a$V1, a$V2, a$V3) %in% peak_ann_promoter, ]
  write.table(a.e, stringr::str_replace(i, ".bed", ".Enhancer.bed"), quote = F, sep = "\t", col.names = F, row.names = F)
}





# gene expr
normalized_rpkm = read.table("/nfs4/chaozhang/proj/Embryo/Yota/RNAseq_Apr2024/script/normalized_rpkm.tab")

genes = lapply(peakAnnoList, function(i){
  g = as.data.frame(i)
  g = g[abs(g$distanceToTSS) < 1000, ]
  g = g$SYMBOL
  unique(g)
})


genes = lapply(1:5, function(i){
  genes[[i]] = setdiff(genes[[i]], unlist(genes[!1:5 %in% i]))
})

genes_df = matrix(nrow = 0, ncol = 2)
genes_df = lapply(1:5, function(i){
  rbind(genes_df, cbind(genes[[i]], i))
})
genes_df = do.call(rbind, genes_df)
genes_df = as.data.frame(genes_df)
colnames(genes_df) = c("gene", "group")

genes_df = cbind(genes_df, normalized_rpkm[genes_df$gene, ])
genes_df = na.omit(genes_df)

#
genes_df = melt(genes_df)

genes_df_avg = genes_df %>%
  dplyr::group_by(group, variable) %>%
  dplyr::summarise(avg = mean(value))

genes_df_avg$samples = sub("_total_rep[1|2]", "", genes_df_avg$variable)
genes_df_avg$samples = factor(genes_df_avg$samples, 
                          levels = c("Brg1_Ctrl", "Brg1_0.3n", "Brg1_1n", "Brg1_3n", "Brg1_10n", "Brg1_100n"))


ggplot(genes_df_avg, aes(samples, avg, fill = samples))+
  geom_bar(stat = "identity", position = "dodge")+
  geom_point()+
  facet_wrap(.~group, ncol = 1) +
  scale_fill_brewer(palette="YlOrRd", direction=-1) +
  cowplot::theme_cowplot()

#

genes_df$variable = sub("_total_rep[1|2]", "", genes_df$variable)
genes_df$variable = factor(genes_df$variable, 
                                   levels = c("Brg1_Ctrl", "Brg1_0.3n", "Brg1_1n", "Brg1_3n", "Brg1_10n", "Brg1_100n"))

df.matrix = dcast(genes_df, group )
ggplot(genes_df, aes(variable, log(value+1), fill = variable))+
  geom_boxplot(outlier.size = -1) + ylim(0,5)+
  facet_wrap(.~group, ncol = 1) +
  scale_fill_brewer(palette="YlOrRd", direction=-1) +
  cowplot::theme_cowplot()
ggsave("../figure/Gene_expr_BRG1_dependent.pdf", width = 12, height = 8)

genes_df_avg_hp = dcast(genes_df, group ~ samples, median, margins = "value")
rownames(genes_df_avg_hp) = genes_df_avg_hp$group
genes_df_avg_hp = genes_df_avg_hp[,-1]
#genes_df_avg_hp = genes_df_avg_hp/genes_df_avg_hp[,1]
pheatmap((genes_df_avg_hp), 
         cluster_rows = F, cluster_cols = F,
         border_color = NA,
         color = pals::)


lapply(1:5, function(i){
  t.test(genes_df[genes_df$group == i & genes_df$variable == "Brg1_Ctrl", "value"],
         genes_df[genes_df$group == i & genes_df$variable == "Brg1_100n", "value"], paired = T)
})




# GO
library(clusterProfiler)
compGO <- compareCluster(geneCluster   = genes,
                         fun           = "enrichGO",
                         keyType = "SYMBOL",
                         ont = "BP",
                         OrgDb = org.Mm.eg.db,
                         pAdjustMethod = "BH")
dotplot(compGO, showCategory = 10, title = "GO Enrichment Analysis")
ggsave("../figure/Gene_expr_BRG1_dependent_GOenrich.pdf", width = 8, height = 6)




#motif
top10 = lapply(list.files("/nfs4/chaozhang/proj/Embryo/Yota/BRG1_ATAC_0901/script/motif_at_Brg1", "knownResults.txt", recursive =T, full.names = T), function(i){
  f = read.table(i,sep="\t", stringsAsFactors = F, comment.char = "*", header = T)
  f$Motif.Name[1:10]
})
top10 = unique(unlist(top10))

top10_d = lapply(list.files("/nfs4/chaozhang/proj/Embryo/Yota/BRG1_ATAC_0901/script/motif_at_Brg1", "knownResults.txt", recursive =T, full.names = T), function(i){
  f = read.table(i,sep="\t", stringsAsFactors = F, comment.char = "*", header = T)
  f = f[!duplicated(f$Motif.Name), ]
  rownames(f) = f$Motif.Name
  f = f[top10, c("P.value", "X..of.Target.Sequences.with.Motif")]
  f$X..of.Target.Sequences.with.Motif = as.numeric(sub("%", "", f$X..of.Target.Sequences.with.Motif))
  f
})

top10_p = do.call(rbind, top10_d)
colnames(top10_p) = c("p", "percentage")
top10_p$group = rep(1:5, each = length(top10))
top10_p$Motif = stringr::str_split_fixed(top10, "/", 3)[,1]
top10_p$p[top10_p$p==0] <- 1e-300
ggplot(top10_p, aes(group, Motif, size = percentage, color = -log10(p)))+
  geom_point()+
  cowplot::theme_cowplot()+
  scale_color_gradientn(colours = c("blue4", "blue", "gray90", "yellow", "red"))
ggsave("../figure//Brg1_5Group_Motif_enrich.pdf", width = 5, height = 5)

#
top10_p = do.call(cbind, top10_d)
top10_p$diff = abs(log(top10_p[,2] / top10_p[,4]))
top10_p = top10_p[order(top10_p$diff, decreasing = T), ]
top10_p = rownames(top10_p[1:20,])

top10_d = lapply(list.files("/nfs4/chaozhang/proj/Embryo/Yota/BRG1_ATAC_0901/script/motif_at_Brg1", "knownResults.txt", recursive =T, full.names = T), function(i){
  f = read.table(i,sep="\t", stringsAsFactors = F, comment.char = "*", header = T)
  f = f[!duplicated(f$Motif.Name), ]
  rownames(f) = f$Motif.Name
  f = f[top10_p, c("P.value", "X..of.Target.Sequences.with.Motif")]
  f$X..of.Target.Sequences.with.Motif = as.numeric(sub("%", "", f$X..of.Target.Sequences.with.Motif))
  f
})

top10_p = do.call(rbind, top10_d)
colnames(top10_p) = c("p", "percentage")
top10_p$group = rep(1:5, each = 20)
top10_p$Motif = stringr::str_split_fixed(rownames(top10_p), "/", 3)[,1]
top10_p$p[top10_p$p==0] <- 1e-300
ggplot(top10_p, aes(group, Motif, size = percentage, color = -log10(p)))+
  geom_point()+
  cowplot::theme_cowplot()+
  scale_color_gradientn(colours = c("blue4", "blue", "gray90", "yellow", "red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
ggsave("../figure//Brg1_independent_Motif_enrich.pdf", width = 5, height = 5)







x = seq(0,1,0.001)
y = x/sqrt(1+x**2) * 1.4
y = 1-x**3
#y = -x**2 + 2*x
pdf("../figure/Sen_Buf_Model.pdf", 5, 5,)
plot(x, y, ylim = c(0,1),pch = 16, cex = .5, xlab= "1-Brg1_level", ylab="Signals")
points(seq(0,1,.001), rev(seq(0,1,.001)),pch = 16, cex = .5, col = "red")
dev.off()

plot(x, y, ylim = c(0,1.3),pch = 16, cex = .5, xlab= "1-Brg1_level", ylab="ATAC_signals")
points(1-Brg1_level,as.numeric(Down_normalized_counts["chr17 79133751 79134071",]),  col = "red")

plot(seq(0,1,.001), rev(seq(0,1,.001)),pch = 16, cex = .5, xlab= "1-Brg1_level", ylab="ATAC_signals")
points(1-Brg1_level,as.numeric(Down_normalized_counts["chr7 82809189 82809323",]), col = "blue")

a = rev(c(0.1,0.3,0.5,0.7,1))
b = c(0.1, 0.5, 0.7, 0.8, 1.1)
points(a, b, col = "red")
points(seq(0,1,.001), rev(seq(0,1,.001)), col = "blue")

dat = as.data.frame(cbind(a, b))
fit = lm(b ~ a, dat)

fit1 = lm(b ~ poly(a, degree = 2), dat)

abline(fit1)

dis1 = function(a, b){
  #exp_b = a/sqrt(1+a**2) * 1.4
  exp_b = 1-a**3
  sqrt(sum((exp_b - b)**2))
}

dis2 = function(a, b){
  exp_b = 1 - a
  sqrt(sum((exp_b - b)**2))
}

dis1(a,b)
dis2(a,b)

dis1(1-Brg1_level,as.numeric(Down_normalized_counts["chr5 137933717 137934542",]))
