library(splines)
library(reshape2)
library(topGO)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(reshape2)
library(scattermore)
library(viridis)

source('/home/whou10/scratch16/whou10/trajectory_variability/function/01_function.R')
Res <- readRDS('/home/whou10/scratch4/whou10/pardoll/nsclc/pseudotime_treg/mypt/res_tde/R3_NR3/pm/res_tde.rds')
stat <- Res$statistics
stat <- stat[stat[, grep('^fdr.*overall$', colnames(stat))] < 0.05, ]
stat <- stat[order(stat[, 1],-stat[, 3]),]
diffgene <- rownames(stat)
str(diffgene)

## get populationFit, and cluster genes
Res$populationFit <-
  getPopulationFit(Res, gene = diffgene, type = 'time')
fit <- Res$populationFit
mat.scale <- scalematrix(fit[diffgene, ,drop=F])
set.seed(12345)
Res$cluster <- mykmeans(mat.scale, maxclunum = 20)$cluster ## auto select k
# set.seed(12345)
# Res$cluster <- kmeans(mat.scale, 5, iter.max = 1000)$cluster
table(Res$cluster)

pdir <- '/home/whou10/scratch4/whou10/pardoll/nsclc/pseudotime_treg/mypt/plot_tde/R3_NR3/pm/'
stat$cluster = Res$cluster[rownames(stat)]
write.csv(stat, paste0(pdir, 'tde_diffgene.csv'))

library(RColorBrewer)
library(ggplot2)
library(gridExtra)
pdf(paste0(pdir, 'tde_hm.pdf'), height = 4, width = 9)
plotFitHm(Res, type = 'time', cellHeightTotal = 200, cellWidthTotal=150, subsampleCell=FALSE )
dev.off()

pdf(paste0(pdir, 'tde_cluster_pattern.pdf'), height = 3, width = 4)
plotClusterMean(testobj = Res,
                cluster = Res$cluster,
                type = 'time')
dev.off()

colnames(Res[[1]])[1] <- 'fdr.overall'
goRes <- GOEnrich(testobj = Res, type = 'time')
pdf(paste0(pdir, 'tde_GO.pdf'), height = 4.5, width = 7)
plotGOEnrich(goRes = goRes)
dev.off()

## plot individual genes
pdf(paste0(pdir, 'top_diffgene.pdf'), width = 8, height = 6)
plotGene(testobj = Res,
           gene = head(diffgene, 20)) 
dev.off()

for (i in unique(Res$cluster)){
  print(i)
  pdf(paste0(pdir, 'top_diffgene_cluster', i, '.pdf'), width = 8, height = 6)
  plotGene(testobj = Res,
           gene = head(diffgene[Res$cluster == i], min(20, sum(Res$cluster==i)))) 
  dev.off()
}

