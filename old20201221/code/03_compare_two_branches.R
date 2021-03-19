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
cnt <- readRDS(here('pseudotime_lung', 'data','diffusion.antigen.rds'))
cnt <- as.matrix(cnt@assays$RNA@counts)

ord1 <- readRDS(paste0(rdir, '/path_znf683_ord.rds'))
ord2 <- readRDS(paste0(rdir, '/path_il7r_ord.rds'))
expr <- cnt[, c(ord1, ord2)]
expr <- expr[rowMeans(expr>0.01) > 0.01, ]

pdt <- data.frame(curve1 = c(seq(1,length(ord1)), rep(0, length(ord2))), 
                  curve2 = c(rep(0, length(ord1)), seq(1, length(ord2))))

cellWeights <- data.frame(curve1 = c(rep(0.99, length(ord1)), rep(0.01, length(ord2))), 
                          curve2 = c(rep(0.01, length(ord1)), rep(0.99, length(ord2))))
rownames(cellWeights) <- rownames(pdt) <- colnames(expr)

set.seed(12345)
sce <- fitGAM(counts = expr, pseudotime = pdt, cellWeights = cellWeights,
            nknots = 6, verbose = FALSE,parallel=T)
saveRDS(sce, paste0(rdir, '/compare_two_branches_sce.rds'))

Final <- list()
endRes <- diffEndTest(sce)  
endRes$fdr <- p.adjust(endRes[,3],'fdr')
endRes <- endRes[endRes$fdr<0.05,]
endRds <- endRes[order(endRes$fdr, -endRes$waldStat), ]
print(nrow(endRes))
write.csv(endRes, paste0(rdir, '/compare_two_branches_endRes.csv'))

pdf(paste0(pdir, '/endRes_fdr_waldStat.pdf'), width = 4, height = 4)
plot(endRes$fdr~endRes$waldStat, pch = 20, main = paste0(nrow(endRes),'genes; PCC=', round(cor(endRes$fdr, endRes$waldStat),2)))
dev.off()

patternRes <- patternTest(sce)  
patternRes$fdr <- p.adjust(patternRes[,3],'fdr')
patternRes <- patternRes[patternRes$fdr<0.05,]
patternRes <- patternRes[order(patternRes$fdr, -patternRes$waldStat), ]
print(nrow(patternRes))
write.csv(patternRes, paste0(rdir, '/compare_two_branches_patternRes.csv'))

pdf(paste0(pdir, '/patternRes_fdr_waldStat.pdf'), width = 4, height = 4)
plot(patternRes$fdr~patternRes$waldStat, pch = 20, main = paste0(nrow(patternRes),'genes; PCC=', round(cor(patternRes$fdr, patternRes$waldStat),2)))
dev.off()

earlyRes <- earlyDETest(sce, knots = c(1,2), global = TRUE)
earlyRes <- earlyRes[!is.na(earlyRes$pvalue), ]
earlyRes$fdr <- p.adjust(earlyRes[,3],'fdr')
earlyRes <- earlyRes[earlyRes$fdr<0.05,]
earlyRes <- earlyRes[order(earlyRes$fdr, -earlyRes$waldStat), ]
print(nrow(earlyRes))
write.csv(earlyRes, paste0(rdir, '/compare_two_branches_earlyRes.csv'))

pdf(paste0(pdir, '/earlyRes_fdr_waldStat.pdf'), width = 4, height = 4)
plot(earlyRes$fdr~earlyRes$waldStat, pch = 20, main = paste0(nrow(earlyRes),'genes; PCC=', round(cor(earlyRes$fdr, earlyRes$waldStat),2)))
dev.off()


## 
reslist <- list(endRes, patternRes, earlyRes)
for (id in 1:length(reslist)){
  plist <- list()
  res <- reslist[[id]]
  print(str(res))
  for (i in rownames(res)[1:20]){
    plist[[i]] <- plotSmoothers(sce, expr, i) + ggtitle(paste0(i, ';fdr=',round(res[i,'fdr'],2),'; waldStat=', round(res[i, 1],2)))
  }
  pdf(paste0(pdir, '/compare_two_branches_', c('endRes', 'patternRes', 'earlyRes')[id], '.pdf'), width = 16, height = 9)
  print(gridExtra::grid.arrange(grobs = plist, nrow=4))
  dev.off()
}
  
## pseudotime DM plot
dfs_meta <- readRDS(here('pseudotime_lung', 'data', 'diffusion.coordinate.meta.rds'))
dm <- dfs_meta[,c('DC1','DC2','CellType')]
dm <- dm[dm$CellType!='TRM-NK like',]

pdf(paste0(pdir, '/pseudotime.pdf'), width = 4, height = 3)
ggplot(data.frame(x = c(dm[ord1,1], dm[ord2,1]), y = c(dm[ord1,2], dm[ord2, 2]), time = c(seq(1, length(ord1)), seq(1, length(ord2))))) + 
        geom_point(aes(x = x, y = y, col = time), size = 0.5) + 
        scale_color_gradientn(colors = colorRampPalette(rev(brewer.pal(11, 'RdYlBu')))(length(c(ord1, ord2)))) + 
        theme_classic() + xlab('DM1') + ylab('DM2')
dev.off()

