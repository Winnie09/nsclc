library(reshape2)
library(TSCAN)
library(scattermore)
library(ggplot2)
library(RColorBrewer)
library(parallel)
library(here)
here()
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
pdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/pardoll/nsclc/pseudotime_lung/plot'
rdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/pardoll/nsclc/pseudotime_lung/result'
## read in data
tf <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/cistrome_allTFs.rds')
dfs_meta <- readRDS(here('pseudotime_lung', 'data', 'diffusion.coordinate.meta.rds'))
# imp <- readRDS('/home-4/zji4@jhu.edu/scratch/TCR/saver/combine/tumornormal.rds') ## very large (182 GB)
# int <- intersect(colnames(imp), dfs_meta$barcode)
# imp <- imp[, int]
# saveRDS(imp, '/scratch/users/whou10@jhu.edu/Wenpin/pardoll/nsclc/data/ManaT/saver.rds')
imp <- readRDS('/scratch/users/whou10@jhu.edu/Wenpin/pardoll/nsclc/data/ManaT/saver.rds') ## use this to read in data

## construct pseudotime
dm <- dfs_meta[,c('DC1','DC2','CellType')]
dm <- dm[dm$CellType!='TRM-NK like',]
v <- as.numeric(as.factor(dm[,3]))
names(v) <- rownames(dm)
mc <- exprmclust(t(dm[,-3]),cluster=v,reduce = T)
ord1 <- TSCANorder(mc, orderonly = T,MSTorder = c(4,1))
ord2 <- TSCANorder(mc, orderonly = T,MSTorder = c(3,2,1))
saveRDS(ord1, '/home-4/whou10@jhu.edu/scratch/Wenpin/pardoll/nsclc/pseudotime_lung/result/path_znf683_ord.rds')
saveRDS(ord2, '/home-4/whou10@jhu.edu/scratch/Wenpin/pardoll/nsclc/pseudotime_lung/result/path_il7r_ord.rds')

## difftest
expr1 = imp[,ord1]
expr1 <- expr1[rowMeans(expr1 > 0.01) > 0.01,]
expr2 = imp[,ord2]
expr2 <- expr2[rowMeans(expr2 > 0.01) > 0.01,]

res1 <- difftest(expr1, ord1)
res2 <- difftest(expr2, ord2)
saveRDS(res1, paste0(rdir, '/path_znf683_res.rds'))
saveRDS(res2, paste0(ridr, '/path_il7r_res.rds'))

res1 <- read.csv('/home-4/whou10@jhu.edu/scratch/Wenpin/pardoll/nsclc/pseudotime_lung/result/path_znf683_res.csv', as.is = T, row.names = 1)
res1$is_TF <- ifelse(rownames(res1) %in% tf[,1], 'Yes', 'No')
res1$has_target_gene <- sapply(rownames(res1), function(i) {
  if (res1[i,3] == 'No'){
    'NA'
  } else if (res1[i,3] == 'Yes' & !is.na(tf[tf[,1]==i, 2])){
    'TRUE'
  } else {
    'NA'
  }
})
write.csv(res1, '/home-4/whou10@jhu.edu/scratch/Wenpin/pardoll/nsclc/pseudotime_lung/result/path_znf683_res.csv')

res2 <- read.csv('/home-4/whou10@jhu.edu/scratch/Wenpin/pardoll/nsclc/pseudotime_lung/result/path_il7r_res.csv', as.is = T, row.names = 1)
res2$is_TF <- ifelse(rownames(res2) %in% tf[,1], 'Yes', 'No')
res2$has_target_gene <- sapply(rownames(res2), function(i) {
  if (res2[i,3] == 'No'){
    'NA'
  } else if (res2[i,3] == 'Yes' & !is.na(tf[tf[,1]==i, 2])){
    'TRUE'
  } else {
    'NA'
  }
})
write.csv(res2, '/home-4/whou10@jhu.edu/scratch/Wenpin/pardoll/nsclc/pseudotime_lung/result/path_il7r_res.csv')

## plot ----------
png('/home-4/whou10@jhu.edu/scratch/Wenpin/pardoll/nsclc/pseudotime_lung/plot/path_znf683_differential_genes_hm.png', width = 600, height = 800)
gl1 = rownames(res1[res1[,2]<0.05,])
# fit <-  get_spline_fit(trainData = expr1[gl, ], trainX = seq(1, ncol(expr1)), fit.min = 1, fit.max = ncol(expr1), num.base = 3)
fit1 <- get_gam_fit(expr1, ord1)  
colnames(fit1) <- seq(1, ncol(fit1))
print(mySTIP2(fit1[gl1,], gl1))
dev.off()

pdf('/home-4/whou10@jhu.edu/scratch/Wenpin/pardoll/nsclc/pseudotime_lung/plot/path_znf683_differential_genes_hm.pdf', width = 6, height = 8)
print(mySTIP2(fit1[gl1,], gl1))
dev.off()


png(paste0(pdir, '/path_il7r_differential_genes_hm.png'), width = 600, height = 1500)
gl2 = rownames(res2[res2[,2]<0.05,])
fit2 <- get_gam_fit(expr2, ord2)  
colnames(fit2) <- seq(1, ncol(fit2))
mySTIP2(fit2[gl2,], gl2)
dev.off()

pdf(paste0(pdir, '/path_il7r_differential_genes_hm.pdf'), width = 6, height = 15)
mySTIP2(fit2[gl2,], gl2)
dev.off()
  
pdf(paste0(pdir, '/path_znf683_pseudotime.pdf'), width = 4, height = 3)
ord = ord1
print(ggplot(data.frame(x = c(dm[ord,1], dm[!rownames(dm) %in% ord,1]), y = c(dm[ord,2], dm[!rownames(dm) %in% ord, 2]), time = c(seq(1, length(ord)), rep(NA, sum(!rownames(dm) %in% ord))))) + 
        geom_point(aes(x = x, y = y, col = time), size = 0.5) + 
        scale_color_gradientn(colors = c(colorRampPalette(rev(brewer.pal(11, 'RdYlBu')))(length(ord)), 'grey')) + 
        theme_classic() + xlab('DM1') + ylab('DM2'))
dev.off()

pdf(paste0(pdir, '/path_il7r_pseudotime.pdf'), width = 4, height = 3)
ord = ord2
print(ggplot(data.frame(x = c(dm[ord,1], dm[!rownames(dm) %in% ord,1]), y = c(dm[ord,2], dm[!rownames(dm) %in% ord, 2]), time = c(seq(1, length(ord)), rep(NA, sum(!rownames(dm) %in% ord))))) + 
        geom_point(aes(x = x, y = y, col = time), size = 0.5) + 
        scale_color_gradientn(colors = c(colorRampPalette(rev(brewer.pal(11, 'RdYlBu')))(length(ord)), 'grey')) + 
        theme_classic() + xlab('DM1') + ylab('DM2'))
dev.off()

## plot marker genes on both paths
ord = c(ord1, ord2)
gs <- c('PRDM1', 'ZNF683', 'ZNF684') ## 'BLIMP1', 'HOBIT' not in rownames
expr.tmp <- cbind(expr1[gs,], expr2[gs,])
expr.tmp <- expr.tmp[,ord]
for (g in gs){
  pdf(paste0(pdir, '/', g, '_expression_along_pseudotime.pdf'), width = 3.3, height = 2.2)
  v = expr.tmp[g, ]
  v[v > quantile(v, 0.95)] <- quantile(v, 0.95)
  v[v < quantile(v, 0.05)] <- quantile(v, 0.05)
  print(ggplot(data.frame(x = c(dm[ord,1], dm[!rownames(dm) %in% ord,1]), y = c(dm[ord,2], dm[!rownames(dm) %in% ord, 2]), expression = v)) + 
          geom_point(aes(x = x, y = y, col = expression), size = 0.5) + 
          scale_color_gradientn(colors = c(colorRampPalette(rev(brewer.pal(11, 'RdYlBu')))(ncol(expr.tmp)), 'grey')) + 
          theme_classic() + xlab('DM1') + ylab('DM2') +
          ggtitle(g))
  dev.off()
}

## GO analysis
res1 <- readRDS(paste0(rdir, '/path_znf683_res.rds'))
gl1 = rownames(res1[res1[,2]<0.05,])
fit1 <- get_gam_fit(expr1, ord1)  
colnames(fit1) <- seq(1, ncol(fit1))
geneorder1 <- mySTIP2(fit1[gl1,], gl1, ReturnGeneOrder = T)
clu1 <- geneorder1[seq(1,which(geneorder1 == 'IFNG'))]
clu2 <- setdiff(gl1, clu1)

res.go <- myGO(clu1, rownames(expr1))
res.go <- res.go[res.go$FDR<0.05, ]
res.go <- res.go[order(res.go$FDR, -res.go$FC), ]
print(res.go[1:min(20, nrow(res.go)),c(2,7,8)])
write.csv(res.go, paste0(rdir, '/path_znf683_GO_cluster1_increasing.csv'))


res.go <- myGO(clu2, rownames(expr1))
res.go <- res.go[res.go$FDR<0.05, ]
res.go <- res.go[order(res.go$FDR, -res.go$FC), ]
print(res.go[1:min(20, nrow(res.go)),c(2,7,8)])
write.csv(res.go, paste0(rdir, '/path_znf683_GO_cluster2_decreasing.csv'))



res2 <- readRDS(paste0(rdir, '/path_il7r_res.rds'))
gl2 = rownames(res2[res2[,2]<0.05,])
fit2 <- get_gam_fit(expr2, ord2)  
colnames(fit2) <- seq(1, ncol(fit2))
geneorder2 <- mySTIP2(fit2[gl2,], gl2, ReturnGeneOrder = T)
clu1 <- geneorder2[seq(1,which(geneorder2 == 'TMSB4X'))]
clu2 <- setdiff(gl2, clu1)

res.go <- myGO(clu1, rownames(expr2))
res.go <- res.go[res.go$FDR<0.05, ]
res.go <- res.go[order(res.go$FDR, -res.go$FC), ]
print(res.go[1:min(20, nrow(res.go)),c(2,7,8)])
write.csv(res.go, paste0(rdir, '/path_il7r_GO_cluster1_increasing.csv'))


res.go <- myGO(clu2, rownames(expr2))
res.go <- res.go[res.go$FDR<0.05, ]
res.go <- res.go[order(res.go$FDR, -res.go$FC), ]
print(res.go[1:min(20, nrow(res.go)),c(2,7,8)])
write.csv(res.go, paste0(rdir, '/path_il7r_GO_cluster2_decreasing.csv'))


