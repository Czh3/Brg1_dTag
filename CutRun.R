setwd("/nfs4/chaozhang/proj/Embryo/Yota/BRG1_CutRun_Apr2024/script")



expr = read.table("../peaks/Brg1.peaks.read.bed", header = F)
expr = expr[!duplicated(expr$V4),]
expr = expr[expr$V1 != "chrM",]
rownames(expr) = paste(expr$V1, expr$V2, expr$V3)

expr = expr[, -c(1:9)]
#colnames(expr) = gsub(".mm.sort.rmdup.bam", "", list.files("../mapping", ".mm.sort.rmdup.bam$"))
colnames(expr) = c("0.3n_rep1","0.3n_rep2", "100n_rep1","100n_rep2","10n_rep1","10n_rep2", "1n_rep1","1n_rep2","3n_rep1","3n_rep2","Ctrl_rep1", "Ctrl_rep2")

grp = stringr::str_replace(colnames(expr), "_rep[1|2]", "")

y <- DGEList(counts=expr, group=grp)
y$samples

libSize = sapply(list.files("../mapping/", ".rmdup.metrics", full.names = T), function(i){
  f = read.table(i, comment.char = "#", sep = "?")
  as.numeric(stringr::str_split_fixed(f[2,], "\t", n=10)[,3])
})

y$samples$lib.size = libSize

y <- calcNormFactors(y, method = "none")
y$samples

normalized_counts = cpm(y, normalized.lib.sizes=T, log=F, prior.count = 1)
normalized_counts_group = cpmByGroup(y, normalized.lib.sizes=T, log=T, prior.count = 1)

pdf("../plots/CunRun_scatter_cor.pdf", 5, 4)
lapply(seq(1, ncol(normalized_counts), 2), function(i){
  mat = log(normalized_counts+1)
  smoothScatter(mat[,i:(i+1)], main = sprintf("cor: %.3f", c = cor(mat[,i], mat[,i+1])))
  abline(0, 1, lty = 2, col = "red")
})
dev.off()


pca = prcomp(t(normalized_counts), scale = T)
pca = as.data.frame(pca$x)

pca = cbind(grp, pca)
pca$grp = factor(pca$grp, levels = c("WT", "Ctrl", "0.3n", "1n", "3n", "10n", "100n"))
 
pca$sample = rownames(pca)
ggplot(pca, aes(PC1, PC2, label = sample)) +
  geom_point(aes(fill = grp),shape = 21, colour = "white", size=6)+
  cowplot::theme_cowplot() +
  scale_fill_brewer(palette="Dark2") +
  geom_text_repel(min.segment.length = Inf, box.padding = 0.5)
ggsave("../plots/CunRun_PCA.pdf", width = 7, height = 5)


# DE Peaks in Brg1 loc
DAR_call = function(c1, c2){
  c1 = "0.3n"
  c2 = "Ctrl"
  
  cond = matrix(c(c1, "rep1",
                  c1, "rep2",
                  c2, "rep1",
                  c2, "rep2"), ncol = 2, byrow = T)
  
  rownames(cond) = paste(cond[,1], cond[,2], sep = "_")
  colnames(cond) = c("condition", "rep")
  cond = as.data.frame(cond)
  

  y2 <- y[,rownames(cond)]
  y2 <- calcNormFactors(y2, method = "none")
  y2$samples
  
  cond$condition = factor(cond$condition, levels = c(c2, c1))
  design <- model.matrix(~cond$condition)
  
  y2 <- estimateDisp(y2, design)
  
  fit <- glmFit(y2,design)
  lrt <- glmLRT(fit)
  Brg1_DEGs = as.data.frame(topTags(lrt, n = nrow(y2$counts)))
  Brg1_DEGs.sig = Brg1_DEGs[Brg1_DEGs$FDR < 0.05, ]
  
  Brg1_DEGs$sig = "none"
  Brg1_DEGs$sig[Brg1_DEGs$logFC > log2(1.5) & rownames(Brg1_DEGs) %in% rownames(Brg1_DEGs.sig)] = "Up"
  Brg1_DEGs$sig[Brg1_DEGs$logFC < -log2(1.5) & rownames(Brg1_DEGs) %in% rownames(Brg1_DEGs.sig)] = "Down"
  
  
  ggplot(as.data.frame(normalized_counts_group), aes(get(c2), get(c1))) +
    geom_point(size = .1, alpha = .3) +
    geom_point(data = as.data.frame(normalized_counts_group[rownames(Brg1_DEGs)[Brg1_DEGs$sig == "Down"],]), size = .5, color = "blue") +
    geom_point(data = as.data.frame(normalized_counts_group[rownames(Brg1_DEGs)[Brg1_DEGs$sig == "Up"],]), size = .5, color = "red") +
    geom_abline(slope = 1) + xlim(1,10) + ylim(1,10) + xlab(c2) + ylab(c1)+
    scale_fill_viridis_c() + theme_classic(base_size = 15) 
}

p1 = DAR_call("0.3n", "Ctrl")
p2 = DAR_call("1n", "Ctrl")
p3 = DAR_call("3n", "Ctrl")
p4 = DAR_call("10n", "Ctrl")
p5 = DAR_call("100n", "Ctrl")
p1|p2|p3|p4|p5


# peak number
num = sapply(c("Ctrl", "0_3n",  "1n", "3n", "10n", "100n"), function(i){
  rep1 = paste0("../peaks/Brg1_", i ,"_CR_rep1/Brg1_", i,"_CR_rep1_peaks.broadPeak")
  rep2 = paste0("../peaks/Brg1_", i ,"_CR_rep2/Brg1_", i,"_CR_rep2_peaks.broadPeak")
  system(paste("/netscr/chaozhang/miniconda3/bin/intersectBed  -a ", rep1, " -b ", rep2,
               "| /netscr/chaozhang/miniconda3/bin/intersectBed  -a - -b /nfs4/chaozhang/proj/Aging/HSC/ATACseq/script/mm10-blacklist.v2.bed.gz -v -wa | wc -l", sep = " "),
         intern = TRUE)
})

pn = data.frame(num)
pn$sample = c("Ctrl", "0.3n", "1n", "3n", "10n", "100n")
pn$sample = factor(pn$sample, level = c("Ctrl", "0.3n", "1n", "3n", "10n", "100n"))
pn$num = as.numeric(pn$num)
ggplot(pn, aes(sample, num, fill = sample)) +
  geom_bar(stat = "identity") +
  #geom_text(aes(label=num), vjust= -0.5) +
  theme_classic()+ ylab("#Brg1 CunRun peaks")+
  cowplot::theme_cowplot()+
  scale_fill_brewer(palette = "OrRd", direction=-1)
ggsave("../plots/Brg1_Peak_number.pdf", width = 8, height = 7)



#model

Brg1_level = c(.61, .61, .03, .03, .09, .09, .42, .42, .22, .22,.78, .78 )
Brg1_level = Brg1_level/0.78

counts_norm_Ctrl = normalized_counts/rowMeans(normalized_counts[,c("Ctrl_rep1", "Ctrl_rep2")])


Zscores = apply(counts_norm_Ctrl, 1, function(x) (x-mean(x))/sd(x))
summary(rowMax(t(Zscores)))
counts_norm_Ctrl = counts_norm_Ctrl[rowMax(Zscores) <= 3, ]


dis1 = function(a, b){
  #exp_b = a/sqrt(1+a**2) * 1.4
  exp_b = 1-a**3
  sqrt(sum((exp_b - b)**2))
}

dis2 = function(a, b){
  exp_b = 1 - a
  sqrt(sum((exp_b - b)**2))
}


res = apply(counts_norm_Ctrl, 1, function(i){
  dis1(1-Brg1_level, i) - dis2(1-Brg1_level, i)
})

#

pheatmap::pheatmap(counts_norm_Ctrl[order(res), c(11,12,1,2,7,8,9,10,5,6,3,4)],
         cluster_cols = F, cluster_rows = F, show_rownames = F,
         scale = "row",
         color = pals::gnuplot(20))



ggplot(as.data.frame(res), aes(x=res)) + 
  geom_histogram(aes(y=..count..), colour="black",  fill = "gray80")+
  #geom_density(alpha=.2, fill="blue") +
  geom_vline(xintercept = 0, color = "red3") +
  cowplot::theme_cowplot() 

write.table(as.data.frame(res), "Brg1_CR_delta.csv", quote = F, col.names = F)

pheatmap(counts_norm_Ctrl[order(rowMeans(counts_norm_Ctrl[,1:2])), c(11,12,1,2,7,8,9,10,5,6,3,4)],
         cluster_cols = F, cluster_rows = F, show_rownames = F,
         scale = "row",
         color = pals::gnuplot(20))

breaksList = seq(0, 2, by = .1)
pheatmap(log(normalized_counts+1)[order(rowMeans(normalized_counts[,1:2])), c(11,12,1,2,7,8,9,10,5,6,3,4)],
         cluster_cols = F, cluster_rows = F, show_rownames = F,
         #scale = "row",
         color = pals::gnuplot(length(breaksList)),
         breaks = breaksList)


breaksList = seq(0, 3, by = .1)
pheatmap(normalized_counts[order(rowMeans(normalized_counts[,9:10])), c(11,12,1,2,7,8,9,10,5,6,3,4)],
         cluster_cols = F, cluster_rows = F, show_rownames = F,
         #scale = "row",
         color = pals::gnuplot(length(breaksList)),
         breaks = breaksList)


pheatmap(normalized_counts_group[order(normalized_counts_group[,4]), c(6,1,4,5,3,2)],
         cluster_cols = F, cluster_rows = F, show_rownames = F,
         scale = "row",
         color = rainbow(20))


xx = seq(0,1,0.001)
yy = 1-xx**3
plot(xx, yy, ylim = c(0,1.3),pch = 16, cex = .5, xlab= "1-Brg1_level", ylab="expr")
points(seq(0,1,.001), rev(seq(0,1,.001)),pch = 16, cex = .5, col = "gray")
points(1-Brg1_level,as.numeric(counts_norm_Ctrl["chr1 3482760 3483419",]),  col = "red")

points(1-Brg1_level,as.numeric(Down_normalized_rpkm["Serpinb1a",]), col = "blue")




# 3 group
Grp_n = nrow(normalized_counts)/3
Brg1_grp = rownames(normalized_counts)[order(rowMeans(normalized_counts[,1:2]))]
Brg1_grp = lapply(1:3, function(i){
  Brg1_grp[((i-1)*Grp_n):(i*Grp_n-1)]
})

write.table(rownames(normalized_counts)[order(rowMeans(normalized_counts[,1:2]))],
            "order_normalized_counts.bed", quote = F, col.names = F, row.names = F)


breaksList = seq(0, 2, by = .1)
pheatmap(log(normalized_counts+1)[order(rowMeans(normalized_counts[,1:2])), c(11,12,1,2,7,8,9,10,5,6,3,4)],
         cluster_cols = F, cluster_rows = F, show_rownames = F,
         gaps_row = seq(1, nrow(normalized_counts), length.out = 4)[-1],
         color = pals::gnuplot(length(breaksList)),
         breaks = breaksList)


normalized_counts_group = cpmByGroup(y, normalized.lib.sizes=T, log=F, prior.count = 1)
normalized_counts_group = as.data.frame(normalized_counts_group)
normalized_counts_group$group = "G1"
normalized_counts_group[rownames(normalized_counts_group) %in% Brg1_grp[[2]], "group"] <- "G2"
normalized_counts_group[rownames(normalized_counts_group) %in% Brg1_grp[[3]], "group"] <- "G3"

normalized_counts_group.m = melt(normalized_counts_group)
normalized_counts_group.m$variable = factor(normalized_counts_group.m$variable, levels = c("Ctrl", "0.3n", "1n", "3n", "10n", "100n"))
ggplot(normalized_counts_group.m, aes(variable, value, fill = variable)) +
  geom_boxplot()+
  facet_wrap( group~., scales = "free", ncol = 1)

write.table(normalized_counts_group, "./group3/Grp3_weak_strong_peaks.bed", quote = F, col.names = F)
  
library(locfit)
normalized_counts_group.m$variable = as.numeric(normalized_counts_group.m$variable)

ggplot(normalized_counts_group.m, aes(variable, value, fill = group)) +
  geom_smooth(aes(variable, value, color = group), method = "locfit",method.args = list(deg=2, alpha=1), fill = "gray80") +
  cowplot::theme_cowplot() + xlab("Brg1 Expr (%)") + ylab("log(Brg1 signal)") + 
  scale_color_manual(values = c( pals::gnuplot(9)[4:9])) +
  facet_wrap( group~., scales = "free", ncol = 1)

ggsave("../plots//Brg1_3Group_changes.linePlot.pdf", width = 4, height = 12)


#
peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb,
                       annoDb="org.Mm.eg.db",
                       tssRegion=c(-1000, 1000), verbose=FALSE)


library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

peakAnno = annotatePeak("../peaks/Brg1_Ctrl_peaks.bed", TxDb=txdb,
             annoDb="org.Mm.eg.db",
             tssRegion=c(-1000, 1000),)

plotAnnoPie(peakAnno)
ggsave("../plots/Brg1_peaks_Features.pdf", width = 7, height = 5)

# chromHMM

Hmm = read.table("/nfs4/chaozhang/proj/Embryo/Yota/BRG1_ATAC_0901/script/chromHMM/Brg1_All_Ctrl.chrHMM.bed")
Hmm_shuff = read.table("/nfs4/chaozhang/proj/Embryo/Yota/BRG1_ATAC_0901/script/chromHMM/Brg1_All_shuff.chrHMM.bed")

df = cbind(as.data.frame(table(Hmm$V13)), as.data.frame(table(Hmm_shuff$V13)))
colnames(df) = c("Elements", "c1", "1", "c2")
df$ratio = log2(df$c1/df$c2)

df$Elements = stringr::str_split_fixed(df$Elements, "_", 2)[,2]
ggplot(df, aes(Elements, ratio, fill = Elements)) +
  geom_bar(stat = "identity") +
  cowplot::theme_cowplot()+ ylab("Enrich ratio: log2(obs/random)")+
  scale_fill_manual(values = as.character(pals::alphabet())) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position="none") 
  

# motif enrich in 3 Groups
top10 = lapply(list.files("group3", "knownResults.txt", recursive =T, full.names = T), function(i){
  f = read.table(i,sep="\t", stringsAsFactors = F, comment.char = "*", header = T)
  f$Motif.Name[1:10]
})
top10 = unique(unlist(top10))
 
top10_p = lapply(list.files("group3", "knownResults.txt", recursive =T, full.names = T), function(i){
  f = read.table(i,sep="\t", stringsAsFactors = F, comment.char = "*", header = T)
  f = f[!duplicated(f$Motif.Name), ]
  rownames(f) = f$Motif.Name
  f = f[top10, c("P.value", "X..of.Target.Sequences.with.Motif")]
  f$X..of.Target.Sequences.with.Motif = as.numeric(sub("%", "", f$X..of.Target.Sequences.with.Motif))
  f
})

top10_p = do.call(rbind, top10_p)
colnames(top10_p) = c("p", "percentage")
top10_p$group = rep(c("G1", "G2", "G3"), each = length(top10))
top10_p$Motif = stringr::str_split_fixed(top10, "/", 3)[,1]
ggplot(top10_p, aes(group, Motif, size = percentage, color = -log10(p)))+
  geom_point()+
  cowplot::theme_cowplot()+
  scale_color_gradientn(colours = c("blue4", "blue", "gray90", "yellow", "red"))
ggsave("../plots/Motif_enrich_3Groups.pdf", width = 7, height = 5)

# chromHMM on 3 group
HMM = lapply(1:3, function(i){
  f = read.table(paste0("group3/Grp3_weak_strong_peaks_G",i,".chrHMM.bed"))
  n = nrow(read.table(paste0("group3/Grp3_weak_strong_peaks_G",i,".bed")))
  table(f$V14)/n
})

df = as.data.frame(do.call(cbind, HMM))
colnames(df) = c( "G1", "G2", "G3")

pdf("../../figure/Brg1_down_group_chrHMM.pdf", 5, 5)
pheatmap::pheatmap(df[c(1,4:10,2,3),],cluster_rows = F, cluster_cols = F,
                   display_numbers = T,
                   border_color = "white",
                   number_format = "%.2f", 
                   fontsize = 15,
                   number_color = "black",
                   color = colorRampPalette(c("white", "yellow", "red"))(100))
dev.off()


#GO

files = list.files("group3", "Grp3_weak_strong_peaks_G.*.bed$", full.names = T)

peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb,
                       annoDb="org.Mm.eg.db",
                       tssRegion=c(-1000, 1000), verbose=FALSE)


plotAnnoBar(peakAnnoList)

genes= lapply(peakAnnoList, function(i) {
  i = as.data.frame(i)
  i = i[i$annotation == "Promoter", ]
  unique(i$SYMBOL)
})
names(genes) = c("G1", "G2", "G3")
vennplot(genes)

library(clusterProfiler)
compGO <- compareCluster(geneCluster   = genes,
                         fun           = "enrichGO",
                         keyType = "SYMBOL",
                         ont = "BP",
                         OrgDb = org.Mm.eg.db,
                         pAdjustMethod = "BH")
dotplot(compGO, showCategory = 10, title = "GO Enrichment Analysis")

# expr
expr = read.table("/nfs4/chaozhang/proj/Embryo/Yota/RNAseq_Apr2024/script/normalized_rpkm.tab")

expr_g3 = lapply(1:3, function(i){
  e = expr[rownames(expr) %in% genes[[i]], 11:12]
  e$group = paste0("G", i)
  e
})
expr_g3 = do.call(rbind, expr_g3)
expr_g3$expr = rowMeans(expr_g3[,1:2])

ggplot(expr_g3, aes(group, log(expr+1), fill = group))+
  geom_boxplot()+
  cowplot::theme_cowplot()+
  scale_fill_brewer(palette="Dark") 
ggsave("../plots/Expr_3Groups.pdf", width = 7, height = 5)

