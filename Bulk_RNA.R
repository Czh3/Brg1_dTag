# bulk RNAseq

setwd("/nfs4/chaozhang/proj/Embryo/Yota/RNAseq_Apr2024/script")
library(edgeR)
library(ggplot2)
library(pheatmap)
library(reshape2)
library(ggrepel)

read_count = read.table("../featureCount/counts.txt", row.names = 1, header = T)
gene_feature = read_count[, 1:5]
read_count = read_count[,-c(1:5)]

# set sample names
colnames(read_count) = gsub(".*mapping.", "", colnames(read_count))
colnames(read_count) = gsub("_S[0-9]*.sort.bam", "", colnames(read_count))

colnames(read_count)[1] = "Brg1_0.3n_total_rep1"
colnames(read_count)[2] = "Brg1_0.3n_total_rep2"

coldata = data.frame(names = colnames(read_count))
rownames(coldata) = colnames(read_count)
cond = strsplit(coldata$names, split = "_")
cond = do.call(rbind, cond)

coldata = as.data.frame(cbind(coldata, cond[,-3]))
colnames(coldata) = c("sample", "cell", "condition", "rep")


#
y <- DGEList(counts=read_count, group=coldata$condition)
y$samples


normalized_rpkm = rpkm(y, gene.length = gene_feature[rownames(y$counts), ]$Length,
                         normalized.lib.sizes=TRUE, prior.count = 1)
write.table(normalized_rpkm, "normalized_rpkm.tab", quote = F)
normalized_rpkm = read.table("normalized_rpkm.tab")
normalized_rpkm_group = rpkmByGroup(y, gene.length = gene_feature[rownames(y$counts), ]$Length,
                       normalized.lib.sizes=TRUE, prior.count = 1)
y_normalized_counts_group = cpmByGroup(y, normalized.lib.sizes=T, log=T, prior.count = 1,
                                        lib.size = y$samples$lib.size)

keep = rowMeans(normalized_rpkm[,11:12]) >= 1 # remove lowly expressed genes

y <- y[keep, , keep.lib.sizes=FALSE]

y <- calcNormFactors(y)
y$samples

normalized_counts = cpm(y, normalized.lib.sizes=TRUE, prior.count = 1)
#normalized_counts <- removeBatchEffect(normalized_counts, batch=coldata$rep)

pca = prcomp(t(normalized_counts), scale = T, center = T)
pca = as.data.frame(pca$x)

pca = cbind(coldata, pca)
pca$sample = paste(pca$condition, pca$rep, sep = "_")

ggplot(pca, aes(PC1, PC2, label = sample)) +
  geom_point(aes(fill = condition),shape = 21, colour = "white", size=6)+
  cowplot::theme_cowplot() +
  scale_fill_brewer(palette="Dark2") +
  geom_text_repel(min.segment.length = Inf, box.padding = 0.5)
ggsave("../plots/RNA_PCA.pdf", width = 7, height = 5)


pdf("../plots/RNA_scatter_cor.pdf", 5, 4)
lapply(seq(1, ncol(normalized_counts), 2), function(i){
  mat = log(normalized_counts+1)
  smoothScatter(mat[,i:(i+1)], main = sprintf("cor: %.3f", c = cor(mat[,i], mat[,i+1])))
  abline(0, 1, lty = 2, col = "red")
})
dev.off()


DAR_call = function(c1, c2){

  cond = matrix(c(c1, "rep1",
                  c1, "rep2",
                  c2, "rep1",
                  c2, "rep2"), ncol = 2, byrow = T)
  
  rownames(cond) = paste("Brg1", cond[,1], "total", cond[,2], sep = "_")
  colnames(cond) = c("condition", "rep")
  cond = as.data.frame(cond)

  y2 <- DGEList(counts=y$counts[, rownames(cond)], 
                group=cond$condition)
  y2$samples
  
  #y2 <- calcNormFactors(y2, method = "upperquartile", p = 0.95)
  y2 <- calcNormFactors(y2)
  
  cond$condition = factor(cond$condition, levels = c(c2, c1))
  design <- model.matrix(~cond$condition)
  
  y2 <- estimateDisp(y2, design)
  
  fit <- glmFit(y2,design)
  lrt <- glmLRT(fit)
  Brg1_DEGs = as.data.frame(topTags(lrt, n = nrow(y2$counts)))
  Brg1_DEGs.sig = Brg1_DEGs[Brg1_DEGs$PValue < 0.05, ]
  
  Brg1_DEGs$sig = "none"
  Brg1_DEGs$sig[Brg1_DEGs$logFC > log2(1.5) & rownames(Brg1_DEGs) %in% rownames(Brg1_DEGs.sig)] = "Up"
  Brg1_DEGs$sig[Brg1_DEGs$logFC < -log2(1.5) & rownames(Brg1_DEGs) %in% rownames(Brg1_DEGs.sig)] = "Down"
  
  
  #write.table(Brg1_DEGs, paste0("Brg1_DEGs_fc1.5_", c1, "_", c2, ".txt"), quote = F, sep = "\t")
  
  
  ggplot(as.data.frame(y_normalized_counts_group), aes(get(c2), get(c1))) +
    geom_point(size = .1, alpha = .3) +
    geom_point(data = as.data.frame(y_normalized_counts_group[rownames(Brg1_DEGs)[Brg1_DEGs$sig == "Down"],,drop=FALSE]), size = .5, color = "blue") +
    geom_point(data = as.data.frame(y_normalized_counts_group[rownames(Brg1_DEGs)[Brg1_DEGs$sig == "Up"],,drop=FALSE]), size = .5, color = "red") +
    geom_abline(slope = 1) + xlab(c2) + ylab(c1)+
    scale_fill_viridis_c() + theme_classic(base_size = 15)
}


p1 = DAR_call("0.3n", "Ctrl")
p2 = DAR_call("1n", "Ctrl")
p3 = DAR_call("3n", "Ctrl")
p4 = DAR_call("10n", "Ctrl")
p5 = DAR_call("100n", "Ctrl")

p1|p2|p3|p4|p5 
ggsave("../plots/RNA_DEGs_fc1.5_scatterPlot.pdf", width = 20, height = 4)


# DEGs number
dTag = c("0.3n", "1n", "3n", "10n", "100n")
Degs_num = lapply(dTag, function(i){
  f = paste0("Brg1_DEGs_fc1.5_", i ,"_Ctrl.txt")
  deg = read.table(f)
  x = as.data.frame(table(deg$sig))
  x$sample = i
  x
})
Degs_num = do.call(rbind, Degs_num)
Degs_num = Degs_num[Degs_num$Var1 != "none",]
Degs_num$sample = factor(Degs_num$sample, levels = dTag)
Degs_num[Degs_num$Var1 == "Down", "Freq"] = -Degs_num[Degs_num$Var1 == "Down", "Freq"]
ggplot(Degs_num, aes(sample, Freq, fill = Var1)) +
  geom_bar(stat = "identity") +
  cowplot::theme_cowplot() + ylab("# DEGs")+
  scale_fill_manual(values = c( "blue3", "red3"))
ggsave("../plots/RNA_DEGs_fc1.5_number.pdf", width = 7, height = 5)





Degs = lapply(dTag, function(i){
  f = paste0("Brg1_DEGs_fc1.5_", i ,"_Ctrl.txt")
  deg = read.table(f)
  deg = deg[deg$sig == "Down",]
  rownames(deg)
})
names(Degs) = dTag

library(UpSetR)
upset(fromList(Degs), order.by = "freq")


specific = function(i){
  a = Degs[[i]]
  b = unlist(Degs[!names(Degs) %in% i])
  setdiff(a, b)
}

sep_10n = specific("10n")
pheatmap(log(normalized_rpkm[sep_10n,]+1))


# GO

library(clusterProfiler)
compGO <- compareCluster(geneCluster   = Degs,
                         fun           = "enrichGO",
                         #universe = rownames(y$counts),
                         keyType = "SYMBOL",
                         OrgDb = org.Mm.eg.db,
                         ont = "BP",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH")
dotplot(compGO, showCategory = 10, title = "GO Enrichment Analysis")


# clustering
all_degs = unique(unlist(Degs))

pheatmap(log(normalized_rpkm+1)[all_degs, c(11,12,1,2,7,8,9,10,5,6,3,4)],
         scale = "row", 
         show_rownames = F, cluster_cols = F,
         color = rainbow(20),
         clustering_method = "ward.D2")



# Modeling


Brg1_level = c(.61, .61, .03, .03, .09, .09, .42, .42, .22, .22,.78, .78 )
Brg1_level = Brg1_level/0.78

Down_normalized_rpkm = normalized_rpkm/rowMeans(normalized_rpkm[,c("Brg1_Ctrl_total_rep1", "Brg1_Ctrl_total_rep2")])
Down_normalized_rpkm = Down_normalized_rpkm[Degs[["100n"]], ] # only focused on DEGs-down

Zscores = apply(Down_normalized_rpkm, 1, function(x) {
  x = as.numeric(x)
  (x-mean(x))/sd(x)
})

Zscores = t(Zscores)
summary(rowMax(Zscores))

Zscores = apply(Down_normalized_rpkm, 1, function(x) (x-mean(x))/sd(x))
summary(rowMax(t(Zscores)))
Down_normalized_rpkm = Down_normalized_rpkm[rowMax(Zscores) <= 3, ]


dis1 = function(a, b){
  #exp_b = a/sqrt(1+a**2) * 1.4
  exp_b = 1-a**3
  sqrt(sum((exp_b - b)**2))
}

dis2 = function(a, b){
  exp_b = 1 - a
  sqrt(sum((exp_b - b)**2))
}


res = apply(Down_normalized_rpkm, 1, function(i){
  dis1(1-Brg1_level, i) - dis2(1-Brg1_level, i)
})

res_df = as.data.frame(res)
colnames(res_df) = "delta"
write.csv(res_df, "Down_DEGs_delta.csv", quote = F)
ggplot(res_df, aes(x=delta)) + 
  geom_histogram(aes(y=..count..), colour="black", fill="gray80")+
  #geom_density(alpha=.2, fill="blue") +
  geom_vline(xintercept = 0, color = "red3") +
  cowplot::theme_cowplot() 
ggsave("../plots/Expr_Model_bias.barplot.pdf", width = 5, height = 5)

sum(res < 0)/length(res)

#
Down_rpkm = normalized_rpkm[names(res), ]
Down_rpkm.m = melt(Down_rpkm)

sum(res < 0)

pdf("../plots/DEGs_Buff_Sen.heatmap.pdf", 5,7)
pheatmap(Down_rpkm[order(res), c(11,12,1,2,7,8,9,10,5,6,3,4)],
         cluster_cols = F, cluster_rows = F, show_rownames = F,
        scale = "row",
        gaps_row = sum(res < 0),
        color = pals::gnuplot(20))
dev.off()

normalized_rpkm_group1 = normalized_rpkm_group/normalized_rpkm_group[,6]

pheatmap(normalized_rpkm_group1[names(res[order(res)]), c(6,1,4,5,3,2)],
         cluster_cols = F, cluster_rows = F, show_rownames = F,
         scale = "row",
         gaps_row = sum(res < 0),
         color = pals::gnuplot(20))


res = res[order(res)]

res.group = list(G1 = names(res[1:sum(res < 0)]),
                 G2 = names(res[(sum(res < 0)+1):length(res)]))


compGO <- compareCluster(geneCluster   = res.group,
                         fun           = "enrichGO",
                         #universe = rownames(y$counts),
                         keyType = "SYMBOL",
                         OrgDb = org.Mm.eg.db,
                         ont = "BP",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH")
dotplot(compGO, showCategory = 10, title = "GO Enrichment Analysis")


#
normalized_rpkm_group = rpkmByGroup(y, gene.length = gene_feature[rownames(y$counts), ]$Length,
                                      normalized.lib.sizes=TRUE, prior.count = 1)

normalized_rpkm_group = as.data.frame(normalized_rpkm_group)
normalized_rpkm_group = normalized_rpkm_group/(normalized_rpkm_group$Ctrl)
normalized_rpkm_group = normalized_rpkm_group[Degs[["100n"]], ] # only focused on DEGs-down

colnames(normalized_rpkm_group) = Brg1_level[seq(1,12,2)]

normalized_rpkm_group$group = "G1"
normalized_rpkm_group[rownames(normalized_rpkm_group) %in% res.group[[2]], "group"] <- "G2"


normalized_rpkm_group.m = melt(normalized_rpkm_group)

library(locfit)
normalized_rpkm_group.m$variable = as.numeric(as.character(normalized_rpkm_group.m$variable))

ggplot(normalized_rpkm_group.m, aes(variable, value, fill = group)) +
  geom_smooth(aes(variable, value, color = group), method = "locfit",method.args = list(deg=2, alpha=1), fill = "gray80") +
  cowplot::theme_cowplot() + xlab("Brg1 Expr (%)") + ylab("Gene expr") + 
  scale_color_manual(values = c( pals::gnuplot(9)[4:9])) +
  facet_wrap( group~., scales = "free", ncol = 1)+ ylim(0,1)+
  scale_x_reverse()
ggsave("../plots//Brg1_3Group_changes.linePlot.pdf", width = 4, height = 12)

ggplot(normalized_rpkm_group.m, aes(variable, value, fill = group)) +
  geom_smooth(aes(variable, value, color = group), method = "locfit",method.args = list(deg=2, alpha=1), fill = "gray80") + #ylim(0,10) +
  cowplot::theme_cowplot() + xlab("Brg1 Expr (%)") + ylab("Gene expr") + 
  scale_color_manual(values = c( pals::gnuplot(5)[2:5])) + ylim(0,1)+
  scale_x_reverse()

ggsave("../plots//Brg1_Buff_Sen_changes.linePlot.pdf", width = 5, height = 5)



xx = seq(0,1,0.001)
yy = 1-xx**3
plot(seq(0,1,.001), rev(seq(0,1,.001)),pch = 16, cex = .5, xlab= "1-Brg1_level", ylab="expr")
points(1-Brg1_level,as.numeric(Down_normalized_rpkm["Serpinh1",]),  col = "blue")

plot(xx, yy, ylim = c(0,1.3),pch = 16, cex = .5, xlab= "1-Brg1_level", ylab="expr")
points(1-Brg1_level,as.numeric(Down_normalized_rpkm["Serpinb1a",]), col = "red")



# WGCNA
library(WGCNA)
options(stringsAsFactors = FALSE);
enableWGCNAThreads()


Down_rpkm_cor = cor(t(Down_rpkm))
pheatmap(Down_rpkm_cor**3)



pca = prcomp((Down_rpkm), scale = T, center = T)
pca = as.data.frame(pca$x)

g=ggplot(pca, aes(PC1, PC2)) +
  geom_point()+
  cowplot::theme_cowplot() 
g



## compare two samples with closest dosages
# Fold change 1.5
p1 = DAR_call("0.3n", "Ctrl")
p2 = DAR_call("1n", "0.3n")
p3 = DAR_call("3n", "1n")
p4 = DAR_call("10n", "3n")
p5 = DAR_call("100n", "10n")

p1|p2|p3|p4|p5 


# DEGs number
dTag = c("Ctrl", "0.3n", "1n", "3n", "10n", "100n")
Degs_num = lapply(2:6, function(i){
  f = paste0("Brg1_DEGs_", dTag[i] ,"_", dTag[i-1],".txt")
  deg = read.table(f)
  x = as.data.frame(table(deg$sig))
  x$sample = paste0( dTag[i] ,"_vs_", dTag[i-1])
  x
})
Degs_num = do.call(rbind, Degs_num)
Degs_num = Degs_num[Degs_num$Var1 != "none",]
Degs_num$sample = factor(Degs_num$sample, levels = sapply(2:6, function(i) paste0( dTag[i] ,"_vs_", dTag[i-1])))
Degs_num[Degs_num$Var1 == "Down", "Freq"] = -Degs_num[Degs_num$Var1 == "Down", "Freq"]
ggplot(Degs_num, aes(sample, Freq, fill = Var1)) +
  geom_bar(stat = "identity") +
  cowplot::theme_cowplot() + ylab("# DEGs")+ xlab("Comparsion")+
  scale_fill_manual(values = c( "blue3", "red3"))
ggsave("../plots/RNA_DEGs_2closeDosage_number.pdf", width = 8, height = 5)


Degs = lapply(2:6, function(i){
  f = paste0("Brg1_DEGs_", dTag[i] ,"_", dTag[i-1],".txt")
  deg = read.table(f)
  deg = deg[deg$sig == "Down",]
  rownames(deg)
})
names(Degs) = sapply(2:6, function(i) paste0( dTag[i] ,"_vs_", dTag[i-1]))

upset(fromList(Degs), order.by = "freq")


compGO <- compareCluster(geneCluster   = Degs,
                         fun           = "enrichGO",
                         #universe = rownames(y$counts),
                         keyType = "SYMBOL",
                         OrgDb = org.Mm.eg.db,
                         ont = "BP",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH")
dotplot(compGO, showCategory = 10, title = "GO Enrichment Analysis")

all_degs = unique(unlist(Degs))
pheatmap(normalized_rpkm[all_degs, ], scale = "row", show_rownames = F,
         cluster_rows = F,
         clustering_method = "ward.D2")





### clustering

normalized_rpkm = read.table("normalized_rpkm.tab")


Degs = lapply(list.files(".", "Brg1_DEGs_fc1.5"), function(i){
  deg = read.table(i)
  deg = deg[deg$sig != "none",]
  rownames(deg)
})
Degs = unique(unlist(Degs))


pheatmap(normalized_rpkm[Degs, c(11,12,1,2,7,8,9,10,5,6,3,4)], show_rownames = F, scale = "row",
         cluster_cols = F,
         clustering_method = "ward.D")


normalized_rpkm = read.table("normalized_rpkm.tab")


Degs = lapply(list.files(".", "Brg1_DEGs_fc1.5"), function(i){
  deg = read.table(i)
  deg = deg[deg$sig != "none",]
  rownames(deg)
})
Degs = unique(unlist(Degs))


pheatmap(normalized_rpkm[Degs, c(11,12,1,2,7,8,9,10,5,6,3,4)],
         show_rownames = F, scale = "row",
         cluster_cols = F,
         clustering_method = "ward.D",
         color = colorRampPalette(c("blue", "white", "red"))(100)
         )

normalized_rpkm_merge = normalized_rpkm[,seq(1,12,2)] + normalized_rpkm[,seq(2,12,2)]

pheatmap(normalized_rpkm_merge[Degs, c(6,1,4,5,3,2)],
         show_rownames = F, scale = "row",
         cluster_cols = F,
         clustering_method = "ward.D",
         color = colorRampPalette(c("blue", "white", "red"))(100)
)


ph = pheatmap(normalized_rpkm[Degs, c(11,12,1,2,7,8,9,10,5,6,3,4)],
         show_rownames = F, scale = "row",
         cluster_cols = F,
         border_color = NA,
         cutree_rows = 5,
         clustering_method = "ward.D2",
         color = colorRampPalette(c("blue", "white", "red"))(100))
ph = cutree(ph$tree_row, k = 5)
ph[ph==1] <- "C4"
ph[ph==2] <- "C3"
ph[ph==3] <- "C1"
ph[ph==4] <- "C2"
ph[ph==5] <- "C5"

pheatmap(normalized_rpkm[names(ph[ph==5]), c(11,12,1,2,7,8,9,10,5,6,3,4)],
         show_rownames = F, scale = "row",
         cluster_cols = F,
         border_color = NA,
         clustering_method = "ward.D",
         color = colorRampPalette(c("blue", "white", "red"))(100))


colnames(normalized_rpkm_merge) = stringr::str_remove(colnames(normalized_rpkm_merge), "_total_rep1")
normalized_rpkm_merge1 = as.data.frame(normalized_rpkm_merge[Degs,])

normalized_rpkm_merge1 = as.data.frame(t(scale(t(normalized_rpkm_merge1))))

normalized_rpkm_merge1$cluster = ph[rownames(normalized_rpkm_merge1)]
normalized_rpkm_merge1$gene = rownames(normalized_rpkm_merge1)
normalized_rpkm_merge1.m = melt(normalized_rpkm_merge1, id.vars = c("cluster", "gene"))
normalized_rpkm_merge1.m$variable = factor(normalized_rpkm_merge1.m$variable,
                                           levels = c("Brg1_Ctrl", "Brg1_0.3n", "Brg1_1n", "Brg1_3n", "Brg1_10n", "Brg1_100n"))
ggplot(normalized_rpkm_merge1.m, aes(variable, value, group = gene))+
  geom_line(alpha = 0.1) +
  geom_smooth(data = normalized_rpkm_merge1.m, aes(variable, value, group = cluster), method='loess', color = "red", fill = "pink", se = T)+
  ylab("Normalized Expr.")+ xlab("")+
  cowplot::theme_cowplot()+
  facet_wrap(cluster ~ ., ncol = 1)


ggplot(normalized_rpkm_merge1.m, aes(variable, value))+
  geom_smooth(data = normalized_rpkm_merge1.m, aes(variable, value, group = cluster, color = cluster), method='loess', color = "red", fill = "pink", se = T)+
  ylab("Normalized Expr.")+ xlab("")+
  cowplot::theme_cowplot()


library(clusterProfiler)
ph1 = as.data.frame(ph)
ph1$gene = rownames(ph1)
compGO <- compareCluster(gene ~ ph,
                         ph1,
                         fun           = "enrichGO",
                         universe = rownames(normalized_rpkm_merge),
                         keyType = "SYMBOL",
                         OrgDb = org.Mm.eg.db,
                         ont = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.1)
dotplot(compGO, showCategory = 5, title = "GO Enrichment Analysis")

library(enrichR)
listEnrichrDbs()
dbs <- c("ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X", "GO_Biological_Process_2018", "WikiPathways_2016", "KEGG_2016")

enriched <- enrichr(ph1[ph1$ph == "C5", "gene"], dbs)
plotEnrich(enriched[[4]], showTerms = 20, numChar = 40)


base_exp =  as.data.frame(normalized_rpkm_merge[Degs, 6, drop = F])
base_exp$cluster = ph[rownames(base_exp)]


ggplot(base_exp, aes(cluster, log(Brg1_Ctrl+1))) +
  geom_boxplot() +
  cowplot::theme_cowplot() 


# Dosage response genes

normalized_rpkm = read.table("normalized_rpkm.tab")
normalized_rpkm = as.data.frame(t(normalized_rpkm))

normalized_rpkm$Brg1_level = 1 - Brg1_level

gene.cor = lapply(1:(ncol(normalized_rpkm)-1), function(g){
  a = cor.test(normalized_rpkm[,g], normalized_rpkm$Brg1_level)
  c(a$p.value, a$estimate)
})
gene.sig = do.call(rbind, gene.cor)
gene.sig = as.data.frame(gene.sig)
rownames(gene.sig) = colnames(normalized_rpkm[,-ncol(normalized_rpkm)])
colnames(gene.sig) = c("p_value", "cor")
gene.sig$p.adj = p.adjust(gene.sig$p_value, method = "BH")
gene.sig = gene.sig[gene.sig$p.adj < 0.05, ]
gene.sig = na.omit(gene.sig)
dim(gene.sig)

gene.sig$DEG = ifelse(gene.sig$cor < 0, "Down", "Up")
gene.sig$gene = rownames(gene.sig)

compGO <- compareCluster(gene ~ DEG,
                         gene.sig,
                         fun           = "enrichGO",
                         universe = rownames(y),
                         keyType = "SYMBOL",
                         OrgDb = org.Mm.eg.db,
                         ont = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.1)
dotplot(compGO, showCategory = 15, title = "GO Enrichment Analysis")

p = pheatmap::pheatmap(normalized_rpkm_group[rownames(gene.sig), c("Ctrl", "0.3n", "1n", "3n", "10n", "100n")], 
                   cluster_cols = F,
                   clustering_method = "ward.D2",
                   cutree_rows = 3,
                   scale = "row", show_rownames = F)
p$tree_row$order

gene_cluster = cutree(p$tree_row, k = 3)
pheatmap::pheatmap(normalized_rpkm_group[names(gene_cluster[gene_cluster==1]), c("Ctrl", "0.3n", "1n", "3n", "10n", "100n")], 
                   cluster_cols = F,
                   clustering_method = "ward.D2",
                   scale = "row", show_rownames = F)

gene_cluster_df = as.data.frame(gene_cluster)
gene_cluster_df$gene = rownames(gene_cluster_df)

compGO <- compareCluster(gene ~ gene_cluster,
                         gene_cluster_df,
                         fun           = "enrichGO",
                         universe = rownames(y),
                         keyType = "SYMBOL",
                         OrgDb = org.Mm.eg.db,
                         ont = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.1)
dotplot(compGO, showCategory = 25, title = "GO Enrichment Analysis")


plot(normalized_rpkm[, c("Spdl1", "Brg1_level")])

