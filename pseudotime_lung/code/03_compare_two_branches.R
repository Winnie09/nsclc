## two branch --------------------
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(slingshot))
suppressMessages(library(tradeSeq))
library(ggplot2)
library(RColorBrewer)
library(here)
here()
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
pdir <- here('pseudotime_lung', 'plot')
rdir <- here('pseudotime_lung', 'result')

## read in data
imp <- readRDS(here('data','ManaT','count.rds')) ## use this to read in data
cnt <- readRDS(here('pseudotime_lung', 'data','diffusion.antigen.rds'))
cnt <- as.matrix(cnt@assays$RNA@counts)

ord1 <- readRDS(paste0(rdir, '/path_znf683_ord.rds'))
ord2 <- readRDS(paste0(rdir, '/path_il7r_ord.rds'))
expr <- cnt[, c(ord1, ord2)]


pdt <- data.frame(curve1 = c(seq(1,length(ord1)), rep(0, length(ord2))), 
                  curve2 = c(rep(0, length(ord1)), seq(1, length(ord2))))

cellWeights <- data.frame(curve1 = c(rep(0.99, length(ord1)), rep(0.01, length(ord2))), 
                          curve2 = c(rep(0.01, length(ord1)), rep(0.99, length(ord2))))
rownames(cellWeights) <- rownames(pdt) <- colnames(expr)
# design = matrix(rep(1,8), nrow=8)
# dimnames(design) = list(paste0('BM',seq(1,8)), 'intercept')
# cellanno = data.frame(cell=colnames(expr), sample = sub(':.*','', colnames(expr)), stringsAsFactors = FALSE)

##
# pdt <- data.frame(curve1 = pseudotime, curve2 = pseudotime)
# rownames(pdt) <- names(pseudotime)
# pdt = pdt[colnames(expr), ]
# 
# v <- (cellanno$sample %in% paste0('BM',seq(1,8)) + 0)
# v <- ifelse(v==1, 0.99, 0.01)
# cellWeights <- data.frame(curve1 = v, curve2 = 1-v)
# rownames(cellWeights) <- colnames(expr)

set.seed(12345)
sce <- fitGAM(counts = expr, pseudotime = pdt, cellWeights = cellWeights,
            nknots = 6, verbose = FALSE,parallel=T)
# sce <- fitGAM(counts = round(expr), pseudotime = pdt, cellWeights = cellWeights,
#               nknots = 6, verbose = FALSE,parallel=TRUE, BPPARAM = MulticoreParam(2))
saveRDS(sce, paste0(rdir, method,'/clusterType', clusterType, '_', pctGene,'_sce.rds'))
Final <- list()
for (TestType in (c('startVsEndTest', 'associationTest'))){
  print(TestType)
  if (grepl('startVsEndTest', TestType)){
    Res <- startVsEndTest(sce)
  } else if (grepl('associationTest', TestType)){
    Res <- associationTest(sce)
  } 
  res <- data.frame(waldStat = Res[,'waldStat'], P.Value = Res[,'pvalue'] ,adj.P.Val = p.adjust(Res$pvalue, method='fdr'))
  row.names(res) <- row.names(Res)
  res <- res[order(res[,3], -res[,1]), ]
  Final[[TestType]] <- res
}
saveRDS(Final, paste0(rdir, method,'/clusterType', clusterType, '_', pctGene,'.rds'))  

