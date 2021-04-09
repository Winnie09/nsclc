rm(list = ls())
library(here)
# setwd(here())
# source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
setwd('/Users/wenpinhou/Dropbox/pardoll/nsclc/')
source('/Users/wenpinhou/Dropbox/trajectory_variability/function/01_function.R')
for (path in c('path1', 'path2', 'path3')){
 #path = 'path2'
  rdir <- paste0('pseudotime_lung/result/', path, '/')
  pdir <- paste0('pseudotime_lung/plot/', path,'/')
  dir.create(pdir, showWarnings = F, recursive = T)
  
  # ---------------------------- #
  # downstream analysis pipeline #  
  # ---------------------------- #
  Res <- readRDS(paste0(rdir, '/testtime_res.rds'))
  
  diffgene = rownames(Res$statistics)[Res$statistics[,grep('^fdr.*overall$', colnames(Res$statistics))] < 0.05]
  str(diffgene)
  ## --------------
  ## population fit
  ## --------------
  Res$populationFit <- getPopulationFit(Res, gene = rownames(Res$statistics), type = 'time')
  
  ## -----------
  ## clustering
  ## -----------
  Res$cluster <- clusterGene(Res, gene = diffgene, type = 'time', k=3)
  
  ## --------------
  ## save diff gene
  ## --------------
  allg <- rownames(Res$statistics[Res$statistics[,1]<0.05,,drop=FALSE])
  res <- Res$statistics[allg, ]
  res <- res[order(res[,1], -res[,3]), ]
  write.csv(cbind(res, cluster = Res$cluster[rownames(res)]), paste0(pdir, 'testtime_differential_genes.csv'))
  
  ## ---------------
  ## plotClusterMean
  ## ----------------
  pdf(paste0(pdir, 'cluster_mean.pdf'), width = 5, height = 3.5)
  plotClusterMean(testobj=Res, cluster = Res$cluster, type = 'time')
  dev.off()
  
  ## -----------
  ## GO analysis
  ## -----------
  goRes <- GOEnrich(testobj = Res, type = 'time', sep = ':.*')
  saveRDS(goRes, paste0(pdir, 'goRes.rds'))
  
  nn <- sapply(1:length(goRes), function(i){
    tmp <- goRes[[i]]
    write.csv(tmp, paste0(pdir, 'cluster', i, '_GO.csv'))
    saveRDS(tmp, paste0(pdir, 'cluster', i, '_GO.rds'))
    tmp <- tmp[tmp[, 'FDR'] < 0.05, ]
    return(0)
    print(nrow(tmp))
  })
  
  if (sum(nn) > 0){
    pdf(paste0(pdir, 'hm_GO_term.pdf'), width = 8, height = 7)
    print(plotGOEnrich(goRes, n=5))
    dev.off()
  }
    
  # ------------------------------------------------------
  # compare original and fitted expression: not tested yet
  # ------------------------------------------------------
  colnames(Res$populationFit) <- colnames(Res$expr)
  png(paste0(pdir, 'fitHm.png'),width = 4000,height = 2500,res = 300)
  print(plotFitHm(Res, subsampleCell = F, cellHeightTotal = 300))
  dev.off()
  
  png(paste0(pdir, 'fitHm_with_genenames.png'),width = 12000,height = 10000,res = 300)
  print(plotFitHm(Res, showRowName = T, subsampleCell = F, cellWidthTotal = 300, cellHeightTotal = length(Res$cluster) * 10))
  dev.off()
  
  gene <- c('GZMA', 'CCL5', 'NKG7', 'GZMK', 'IL6R', 'SELL', 'CD74', 'CCR7', 'IFNGR2', 'TCF7', 'HLA-DPA1', 'IL2RA', 'EOMES', 'IFNG', 'GZMB')
  gene <- gene[gene %in% rownames(Res$populationFit)]
  png(paste0(pdir, 'example_genes.png'),width = 2000,height = 1800, res = 300)
  plotGene(Res, gene, plot.point = T, point.size =0.2)
  dev.off()
  
  png(paste0(pdir, 'example_genes_pullplot.png'),width = 2000,height = 1800, res = 300)
  plotGeneCellAndPopulation(Res, gene = gene)
  dev.off()
}
  


