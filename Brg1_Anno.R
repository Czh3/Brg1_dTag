# Brg1 peak annotation

setwd("/nfs4/chaozhang/proj/Embryo/Yota/BRG1_ATAC_0901/script")

Brg1_peak_ann = read.table("Brg1_peak_ann.bed", quote = "?", sep = "\t")

Brg1_peak_ann$ann = stringr::str_replace_all(Brg1_peak_ann$V26, " .*", "")
Brg1_peak_ann$V2 = Brg1_peak_ann$V2 - 1
Brg1_peak_ann$name = paste(Brg1_peak_ann$V1, Brg1_peak_ann$V2, Brg1_peak_ann$V3)
Brg1_peak_ann$ann[Brg1_peak_ann$ann != "Promoter"] <- "Enhancer"

CTCF = read.table("/nfs4/chaozhang/proj/Embryo/Yota/BRG1_CutRun_Apr2024/script/Brg1_Ctcf.bed")
CTCF$name = paste(CTCF$V1, CTCF$V2, CTCF$V3)

Brg1_peak_ann[Brg1_peak_ann$name %in% CTCF$name & Brg1_peak_ann$ann != "Promoter", "ann"] <- "CTCF"

table(Brg1_peak_ann$ann)

write.table(Brg1_peak_ann[,c(1:3, 38)], "Brg1_peak_ann.all.bed",col.names = F, row.names = F,quote = F, sep = "\t")

write.table(Brg1_peak_ann[Brg1_peak_ann$ann == "Promoter",c(1:3, 38)], "Brg1_peak_ann.Promoter.bed",col.names = F, row.names = F,quote = F, sep = "\t")
write.table(Brg1_peak_ann[Brg1_peak_ann$ann == "Enhancer",c(1:3, 38)], "Brg1_peak_ann.Enhancer.bed",col.names = F, row.names = F,quote = F, sep = "\t")
write.table(Brg1_peak_ann[Brg1_peak_ann$ann == "CTCF",c(1:3, 38)], "Brg1_peak_ann.Insulator.bed",col.names = F, row.names = F,quote = F, sep = "\t")


data = as.data.frame(table(Brg1_peak_ann$ann))
ggplot(data, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() +
  scale_fill_manual(values = c("royalblue", "yellow2", "red"))
ggsave("../figure/Brg1_peak_ann_pie.pdf", width = 5, height = 4)


#
Brg1_CR = read.table("/nfs4/chaozhang/proj/Embryo/Yota/BRG1_CutRun_Apr2024/peaks/Brg1.peaks.read.bed")
Brg1_CR$name = paste(Brg1_CR$V1, Brg1_CR$V2, Brg1_CR$V3)
Brg1_CR = Brg1_CR[,10:22]
colnames(Brg1_CR) =  c("0.3n_rep1","0.3n_rep2", "100n_rep1","100n_rep2","10n_rep1","10n_rep2", "1n_rep1","1n_rep2","3n_rep1","3n_rep2","Ctrl_rep1", "Ctrl_rep2", "name")

Brg1_CR$ann = Brg1_peak_ann[match(Brg1_CR$name, Brg1_peak_ann$name), "ann"]

Brg1_CR[,seq(1, 12, by = 2)] = (Brg1_CR[,seq(1, 12, by = 2)]/Brg1_CR$Ctrl_rep1 + Brg1_CR[,seq(2, 12, by = 2)]/Brg1_CR$Ctrl_rep2)/2
Brg1_CR = Brg1_CR[, -seq(2, 12, by = 2)] 

Brg1_CR.m = melt(Brg1_CR)

Brg1_CR.m = aggregate(Brg1_CR.m, list(ann = Brg1_CR.m$ann, sample = Brg1_CR.m$variable), mean)

Brg1_CR.m$groups = stringr::str_replace(Brg1_CR.m$sample, "_rep[1|2]", "")
Brg1_CR.m$groups  = factor(Brg1_CR.m$groups , levels = c("WT", "Ctrl", "0.3n", "1n", "3n", "10n", "100n"))
Brg1_CR.m = Brg1_CR.m[,c(1,2,6:7)]
Brg1_CR.m$ann = factor(Brg1_CR.m$ann, levels = c("Promoter", "Enhancer", "CTCF"))

ggplot(Brg1_CR.m, aes(groups, value)) +
  geom_point(color = "red2") +
  geom_smooth(aes(as.numeric(groups)-1, value),  color = "red2", method = "lm", se = F)+
  facet_wrap( ann ~., ncol = 1) +
  cowplot::theme_cowplot() 
  
ggsave("../figure/Brg1_peak_ann_CR_line.pdf", width = 5, height = 10)



Brg1_ATAC = read.table("../peaks/Brg_new.ATAC_reads.bed")
Brg1_ATAC$name = paste(Brg1_ATAC$V1, Brg1_ATAC$V2, Brg1_ATAC$V3)
Brg1_ATAC = Brg1_ATAC[,c(10:21, 24)]
colnames(Brg1_ATAC) = c("Ctrl_rep1", "Ctrl_rep2", "100n_rep1","100n_rep2","10n_rep1","10n_rep2", "1n_rep1","1n_rep2","3n_rep1","3n_rep2","0.3n_rep1","0.3n_rep2", "name")

Brg1_ATAC$ann = Brg1_peak_ann[match(Brg1_ATAC$name, Brg1_peak_ann$name), "ann"]
#size_f = as.numeric(apply(Brg1_ATAC[,1:12], 2, mean))
size_f = colSums(Brg1_ATAC[,1:12])
Brg1_ATAC[,1:12] = t(apply(Brg1_ATAC[,1:12], 1, function(x) {
  x / size_f * 10^6
}))
Brg1_ATAC[,seq(1, 12, by = 2)] = (Brg1_ATAC[,seq(1, 12, by = 2)]/(Brg1_ATAC$Ctrl_rep1+ 1) + Brg1_ATAC[,seq(2, 12, by = 2)]/(Brg1_ATAC$Ctrl_rep2+ 1))/2
Brg1_ATAC = Brg1_ATAC[, -seq(2, 12, by = 2)] 

Brg1_ATAC.m = melt(Brg1_ATAC)

Brg1_ATAC.m = aggregate(Brg1_ATAC.m, list(ann = Brg1_ATAC.m$ann, sample = Brg1_ATAC.m$variable), mean)

Brg1_ATAC.m$groups = stringr::str_replace(Brg1_ATAC.m$sample, "_rep[1|2]", "")
Brg1_ATAC.m$groups  = factor(Brg1_ATAC.m$groups , levels = c("WT", "Ctrl", "0.3n", "1n", "3n", "10n", "100n"))
Brg1_ATAC.m = Brg1_ATAC.m[,c(1,2,6:7)]
Brg1_ATAC.m$ann = factor(Brg1_ATAC.m$ann, levels = c("Promoter", "Enhancer", "CTCF"))

ggplot(Brg1_ATAC.m, aes(groups, value)) +
  geom_point(color = "blue") +
  geom_smooth(aes(as.numeric(groups)-1, value),  color = "blue", method = "lm", se = F)+
  facet_wrap( ann ~., ncol = 1) +
  cowplot::theme_cowplot() 

ggsave("../figure/Brg1_peak_ann_ATAC_line.pdf", width = 5, height = 10)


Brg1_CR.m$seq = "CR"
Brg1_ATAC.m$seq = "ATAC"

df = as.data.frame(rbind(Brg1_CR.m, Brg1_ATAC.m))
ggplot(df, aes(groups, value, group = seq, color = seq)) +
  geom_point() +
  geom_smooth(aes(as.numeric(groups)-1, value, group = seq), method = "lm", se = F)+
  facet_wrap( ann ~., ncol = 1) +
  scale_color_manual(values = c(  "blue", "red")) +
  cowplot::theme_cowplot() 

ggsave("../figure/Brg1_peak_ann_ATAC_CR_line.pdf", width = 5, height = 10)


# 100n vs ctrl
DEG_100n = read.table("Brg1_DEGs_100n_Ctrl.txt", sep = "\t")
DEG_100n$ann = Brg1_peak_ann[match(rownames(DEG_100n),Brg1_peak_ann$name), "ann"]

table(DEG_100n$ann, DEG_100n$sig)
y1_normalized_counts = read.csv("ATAC_at_brg1_normalized_counts.csv", row.names = 1)

ggplot(as.data.frame(y1_normalized_counts), aes(Ctrl_rep1, X100n_rep1)) +
  geom_point(size = .1, alpha = .3) +
  geom_point(data = as.data.frame(y1_normalized_counts[rownames(DEG_100n)[DEG_100n$sig == "Down"],]), size = .5, color = "blue") +
  geom_point(data = as.data.frame(y1_normalized_counts[rownames(DEG_100n)[DEG_100n$sig == "Up"],]), size = .5, color = "red") +
  geom_abline(slope = 1) + xlim(1,10) + ylim(1,10) + xlab("Ctrl") + ylab("100n")+
  scale_fill_viridis_c() + theme_classic(base_size = 15) 

y1_normalized_counts$ann = Brg1_peak_ann[match(rownames(y1_normalized_counts),Brg1_peak_ann$name), "ann"]

ggplot(na.omit(as.data.frame(y1_normalized_counts)), aes(Ctrl_rep1, X100n_rep1)) +
  geom_point(size = .1) +
  geom_point(data = as.data.frame(y1_normalized_counts[rownames(DEG_100n)[DEG_100n$sig == "Down"],]), size = .2, color = "blue") +
  geom_point(data = as.data.frame(y1_normalized_counts[rownames(DEG_100n)[DEG_100n$sig == "Up"],]), size = .2, color = "red") +
  geom_abline(slope = 1) + xlim(1,10) + ylim(1,10) + xlab("Ctrl") + ylab("100n")+
  theme_classic(base_size = 15) +
  facet_wrap(ann ~ .)
ggsave("../figure/Brg1_peak_ann_DAR.pair-wise_scatterplot.pdf", width = 15, height = 5)

table(DEG_100n$ann, DEG_100n$sig)
#Down none   Up
#CTCF      400 1836   16
#Enhancer 5813 7282   69
#Promoter  776 3870  184

df = as.data.frame(table(DEG_100n$ann, DEG_100n$sig))

p1 = ggplot(df[df$Var1 == "Promoter",], aes(x="", y=Freq, fill=Var2)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  geom_text(aes(label = scales::percent(Freq/sum(Freq), .01)),
            size = 3,
            position = position_stack(vjust = 0.5)) +
  theme_void() + ggtitle("Promoter") +
scale_fill_manual(values = c("blue", "gray", "red"))

p2 = ggplot(df[df$Var1 == "Enhancer",], aes(x="", y=Freq, fill=Var2)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  geom_text(aes(label = scales::percent(Freq/sum(Freq), .01)),
            size = 3,
            position = position_stack(vjust = 0.5)) +
  theme_void() + ggtitle("Enhancer") +
scale_fill_manual(values = c("blue", "gray", "red"))
p3 = ggplot(df[df$Var1 == "CTCF",], aes(x="", y=Freq, fill=Var2)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  geom_text(aes(label = scales::percent(Freq/sum(Freq), .01)),
            size = 3,
            position = position_stack(vjust = 0.5)) +
  theme_void() + ggtitle("CTCF") +
  scale_fill_manual(values = c("blue", "gray", "red"))
p1|p2|p3
ggsave("../figure/Brg1_peak_ann_DAR.pie.pdf", width = 12, height = 4)



df = as.data.frame(table(DEG_100n[DEG_100n$sig == "Down", ]$ann))

ggplot(df, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  geom_text(aes(label = scales::percent(Freq/sum(Freq), .01)),
            size = 3,
            position = position_stack(vjust = 0.5)) +
  theme_void() + ggtitle("All Down DAR") +
  scale_fill_manual(values = c("royalblue", "yellow2", "red"))
ggsave("../figure/Brg1_peak_ann_Down_DAR.pie.pdf", width = 5, height = 4)


# motif
DEG_100n = na.omit(DEG_100n)
for (i in c("Promoter", "Enhancer", "CTCF")) {
  write.table(stringr::str_split_fixed(rownames(DEG_100n[DEG_100n$ann == i & DEG_100n$sig == "Down" , ]), pattern = " ", 3),
              paste0("motif_dependent_3elements/DAR_Down_", i, ".bed"), 
              sep = "\t", quote = F, row.names = F, col.names = F)
  write.table(stringr::str_split_fixed(rownames(DEG_100n[DEG_100n$ann == i & DEG_100n$sig == "none" , ]), pattern = " ", 3),
              paste0("motif_dependent_3elements/DAR_none_", i, ".bed"), 
              sep = "\t", quote = F, row.names = F, col.names = F)
}


motifs = lapply(list.files("/nfs4/chaozhang/proj/Embryo/Yota/BRG1_ATAC_0901/script/motif_dependent_3elements/", "knownResults.txt", recursive = T, full.names = T), function(f){
  motif_file = read.table(f, sep = "\t", comment.char = "", header = T)
  motif_file[,7] = as.numeric(stringr::str_replace(motif_file[,7], "%", ""))
  motif_file = motif_file[order(motif_file$Motif.Name), ]
  motif_file$P.value = motif_file$P.value + 1e-300
  motif_file[,c(1,3,7)]
  
})
motifs = do.call(cbind, motifs)
#motifs = motifs[apply(motifs[,seq(2,18,3)], 1, mean) < 0.05, ]

motifs_top10 = lapply(list.files("/nfs4/chaozhang/proj/Embryo/Yota/BRG1_ATAC_0901/script/motif_dependent_3elements", "knownResults.txt", recursive = T, full.names = T), function(f){
  motif_file = read.table(f, sep = "\t", comment.char = "", header = T)
  motif_file[1:10 ,1]
})
motifs_top10 = unique(do.call("c", motifs_top10))

motifs_plot = motifs[motifs$Motif.Name %in% motifs_top10,]

motifs_plot_hp = data.frame(motifs_plot[,seq(2,18,3)])
rownames(motifs_plot_hp) = stringr::str_replace_all(motifs_plot$Motif.Name, "/.*", "")
colnames(motifs_plot_hp) = c("Down_Insulator", "Down_Enhancer", "Down_promoter","none_Insulator", "none_Enhancer", "none_promoter")
pheatmap::pheatmap(-log10(motifs_plot_hp), cluster_cols = F)

motifs_plot = data.frame(unlist(motifs_plot[,seq(1,18,3)]), unlist(motifs_plot[,seq(2,18,3)]), unlist(motifs_plot[,seq(3,18,3)]))
colnames(motifs_plot) = c("Motif", "pval", "percentage")

motifs_plot$Group = rep(c("Down_Insulator", "Down_Enhancer", "Down_promoter","none_Insulator", "none_Enhancer", "none_promoter"), each = nrow(motifs_plot)/6)
motifs_plot$Group = factor(motifs_plot$Group, levels = rev(c("none_promoter", "Down_promoter", "none_Enhancer", "Down_Enhancer",  "none_Insulator", "Down_Insulator")))
motifs_plot$Motif = stringr::str_replace_all(motifs_plot$Motif, "\\(.*", "")
#motifs_plot$pval[motifs_plot$pval < 1e-100] <- 1e-100

ggplot(motifs_plot[motifs_plot$Motif %in% c("CTCF", "BORIS", "NFY", "Ronin", "Sox2", "Oct4", "Esrrb", "NANOG"), ],
       aes(Motif, Group)) +
  geom_point(aes(color = -log10(pval),  size = percentage )) +
  scale_color_gradientn(colours = c("gray95", "red", "red4")) +
  #scale_color_gradientn(colours = c("gray95", "yellow2", "orange", "red", "red4")) +
  theme_minimal(base_size = 15) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))

ggsave("../figure/Brg1_peak_ann_Motif_enrich.pdf", width = 8, height = 6)

ggplot(motifs_plot,
       aes(Motif, Group)) +
  geom_point(aes(color = -log10(pval),  size = percentage )) +
  scale_color_gradientn(colours = c("gray95", "yellow2", "orange", "red", "red4")) +
  theme_minimal(base_size = 20) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))

ggsave("../figure/Brg1_peak_ann_Motif_enrich_all.pdf", width = 30, height = 6)


# modeling

DAR_down = read.csv("Brg1_DAR_down_delta.csv", row.names = 1)

ggplot(DAR_down, aes(x=delta)) + 
  geom_histogram(aes(y=..count..), colour="black",  fill = "gray80")+
  #geom_density(alpha=.2, fill="blue") +
  geom_vline(xintercept = 0, color = "red3") +
  cowplot::theme_cowplot() 

ggsave("../figure/Brg1_Model_bias.barplot.pdf", width = 5, height = 5)

DAR_down$ann = Brg1_peak_ann[match(rownames(DAR_down), Brg1_peak_ann$name), "ann"]

ggplot(DAR_down, aes(x=delta)) + 
  geom_histogram(aes(y=..count..), colour="black",  fill = "gray80")+
  #geom_density(alpha=.2, fill="blue") +
  geom_vline(xintercept = 0, color = "red3") +
  cowplot::theme_cowplot() +
  facet_wrap(ann ~ ., scales = "free")

ggsave("../figure/Brg1_Model_bias_ann.barplot.pdf", width = 12, height = 4)

DAR_down$sen = ifelse(DAR_down$delta > 0, "Sensitive", "Buffered")
table(DAR_down$ann, DAR_down$sen)
#Buffered Sensitive
#CTCF          294       105
#Enhancer     2596      3210
#Promoter      552       224


# super enhancer
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

G1.ann = annotatePeak("Brg1_DAR_down_delta.Group_1.SE.bed", TxDb=txdb,
             annoDb="org.Mm.eg.db",
             tssRegion=c(-1000, 1000))

G1.ann = as.data.frame(G1.ann)

library(clusterProfiler)
G1.gene = unique(G1.ann$SYMBOL)
ego = enrichGO(gene = G1.gene,
               OrgDb = org.Mm.eg.db,
               keyType = "SYMBOL",
               ont = "BP")
barplot(ego, showCategory = 20)

ggsave("../figure/SE_at_G1.GO.barplot.pdf", width = 10, height = 6)

normalized_rpkm = read.table("/nfs4/chaozhang/proj/Embryo/Yota/RNAseq_Apr2024/script/normalized_rpkm.tab")
normalized_rpkm.p = na.omit(normalized_rpkm[G1.gene,])
normalized_rpkm.p = normalized_rpkm.p[rowMeans(normalized_rpkm.p)!=0,]
pheatmap::pheatmap(log(normalized_rpkm.p+1)[,c(11,12,1,2,7,8,9,10,5,6,3,4)],scale = "row", cluster_cols = F)

normalized_rpkm.m = melt(normalized_rpkm[G1.gene,])
normalized_rpkm.m$samples = sub("_total_rep[1|2]", "", normalized_rpkm.m$variable)
normalized_rpkm.m$samples = factor(normalized_rpkm.m$samples, 
                              levels = c("Brg1_Ctrl", "Brg1_0.3n", "Brg1_1n", "Brg1_3n", "Brg1_10n", "Brg1_100n"))

ggplot(na.omit(normalized_rpkm.m), aes(samples, log(value+1), fill = samples))+
  geom_boxplot()+
  scale_fill_brewer(palette="YlOrRd", direction=-1) +
  cowplot::theme_cowplot()
ggsave("../figure/SE_at_G1.expr.barplot.pdf", width = 7, height = 5)

#
CR_delta = read.table("/nfs4/chaozhang/proj/Embryo/Yota/BRG1_CutRun_Apr2024/script/Brg1_CR_delta.csv")
CR_delta$name = paste(CR_delta$V1, CR_delta$V2, CR_delta$V3)
CR_delta$ann = Brg1_peak_ann[match(CR_delta$name, Brg1_peak_ann$name), "ann"]

ggplot(as.data.frame(CR_delta), aes(x=V4)) + 
  geom_histogram(aes(y=..count..), colour="black",  fill = "gray80")+
  geom_vline(xintercept = 0, color = "red3") + xlab("")+ #ggtitle(sprintf("Oct4 peaks (%.2f %% Sensitive)", sum(Oct4_delta$V4>0) / nrow(Oct4_delta) * 100)) +
  cowplot::theme_cowplot() +
  facet_wrap(ann ~ ., scales = "free")
ggsave("../figure/Brg1_CR_Pro_Enh.delta.pdf", width = 12, height = 4)


CR_delta$sen = ifelse(CR_delta$V4 > 0, "Sen", "Dep")
table(CR_delta$ann, CR_delta$sen)
#           Dep   Sen
#CTCF       369  1883
#Enhancer   578 12586
#Promoter   640  4190

ggplot(as.data.frame(CR_delta), aes(x=V4)) + 
  geom_histogram(aes(y=..count..), colour="black",  fill = "gray80")+
  geom_vline(xintercept = 0, color = "red3") + xlab("")+ #ggtitle(sprintf("Oct4 peaks (%.2f %% Sensitive)", sum(Oct4_delta$V4>0) / nrow(Oct4_delta) * 100)) +
  cowplot::theme_cowplot() 
ggsave("../figure/Brg1_CR_All.delta.pdf", width = 5, height = 5)

table(CR_delta$sen)
#Dep   Sen 
#1587 18659 

## OCT4 
#intersectBed -a Brg1_CR_delta.bed -b /nfs4/chaozhang/proj/Embryo/Yota/public/histone/peaks/Oct4/Oct4_peaks.narrowPeak -u > Brg1_CR_OCT4.delta.bed
#

Oct4_delta = read.table("/nfs4/chaozhang/proj/Embryo/Yota/BRG1_CutRun_Apr2024/script/Brg1_CR_OCT4.delta.bed")

p1 = ggplot(as.data.frame(Oct4_delta), aes(x=V4)) + 
  geom_histogram(aes(y=..count..), colour="black",  fill = "gray80")+
  geom_vline(xintercept = 0, color = "red3") + xlab("")+ ggtitle(sprintf("Oct4 peaks (%.2f %% Sensitive)", sum(Oct4_delta$V4>0) / nrow(Oct4_delta) * 100)) +
  cowplot::theme_cowplot() 

Esrrb_delta = read.table("/nfs4/chaozhang/proj/Embryo/Yota/BRG1_CutRun_Apr2024/script/Brg1_CR_Esrrb.delta.bed")

p2 = ggplot(as.data.frame(Esrrb_delta), aes(x=V4)) + 
  geom_histogram(aes(y=..count..), colour="black",  fill = "gray80")+
  geom_vline(xintercept = 0, color = "red3") + xlab("")+ ggtitle(sprintf("Esrrb peaks (%.2f %% Sensitive)", sum(Esrrb_delta$V4>0) / nrow(Esrrb_delta) * 100)) +
  cowplot::theme_cowplot() 

Ctcf_delta = read.table("/nfs4/chaozhang/proj/Embryo/Yota/BRG1_CutRun_Apr2024/script/Brg1_CR_CTCF.delta.bed")

p3 = ggplot(as.data.frame(Ctcf_delta), aes(x=V4)) + 
  geom_histogram(aes(y=..count..), colour="black",  fill = "gray80")+
  geom_vline(xintercept = 0, color = "red3") + xlab("")+ ggtitle(sprintf("Ctcf peaks (%.2f %% Sensitive)", sum(Ctcf_delta$V4>0) / nrow(Ctcf_delta) * 100)) +
  cowplot::theme_cowplot() 

NFYA_delta = read.table("/nfs4/chaozhang/proj/Embryo/Yota/BRG1_CutRun_Apr2024/script/Brg1_CR_NFYA.delta.bed")

p4 = ggplot(as.data.frame(NFYA_delta), aes(x=V4)) + 
  geom_histogram(aes(y=..count..), colour="black",  fill = "gray80")+
  geom_vline(xintercept = 0, color = "red3") + xlab("")+ ggtitle(sprintf("Nfya peaks (%.2f %% Sensitive)", sum(NFYA_delta$V4>0) / nrow(NFYA_delta) * 100)) +
  cowplot::theme_cowplot() 

p1|p2|p3|p4

ggsave("../figure/Brg1_CR_TFs.delta.pdf", width = 16, height = 4)


K27ac_delta = read.table("/nfs4/chaozhang/proj/Embryo/Yota/BRG1_CutRun_Apr2024/script/Brg1_H3K27ac.delta.bed")

ggplot(as.data.frame(K27ac_delta), aes(x=V4)) + 
  geom_histogram(aes(y=..count..), colour="black",  fill = "gray80")+
  geom_vline(xintercept = 0, color = "red3") + xlab("")+ ggtitle(sprintf("H3K27ac peaks (%.2f %% Sensitive)", sum(K27ac_delta$V4>0) / nrow(K27ac_delta) * 100)) +
  cowplot::theme_cowplot() 
ggsave("../figure/Brg1_H3K27ac.delta.pdf", width = 5, height = 5)


# motif

motifs = lapply(list.files("/nfs4/chaozhang/proj/Embryo/Yota/BRG1_CutRun_Apr2024/script/motif/", "knownResults.txt", recursive = T, full.names = T), function(f){
  motif_file = read.table(f, sep = "\t", comment.char = "", header = T)
  motif_file[,7] = as.numeric(stringr::str_replace(motif_file[,7], "%", ""))
  motif_file = motif_file[order(motif_file$Motif.Name), ]
  motif_file$P.value = motif_file$P.value + 1e-300
  motif_file[,c(1,3,7)]
  
})
motifs = do.call(cbind, motifs)

motifs_plot = data.frame(unlist(motifs[,seq(1,6,3)]), unlist(motifs[,seq(2,6,3)]), unlist(motifs[,seq(3,6,3)]))
colnames(motifs_plot) = c("Motif", "pval", "percentage")

motifs_plot$Group = rep(c("Buffered", "Sensitive"), each = nrow(motifs_plot)/2)
motifs_plot$Motif = stringr::str_replace_all(motifs_plot$Motif, "\\(.*", "")
#motifs_plot$pval[motifs_plot$pval < 1e-100] <- 1e-100

ggplot(motifs_plot[motifs_plot$Motif %in% c("CTCF", "BORIS", "NFY", "Ronin", "Sox2", "Oct4", "Esrrb", "Klf4"), ],
       aes(Motif, Group)) +
  geom_point(aes(color = -log10(pval),  size = percentage )) +
  scale_color_gradientn(colours = c("gray95", "red", "red4")) +
  #scale_color_gradientn(colours = c("gray95", "yellow2", "orange", "red", "red4")) +
  theme_minimal(base_size = 15) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))

ggsave("../figure/Brg1_CR_delta_Motif_enrich.pdf", width = 10, height = 6)

