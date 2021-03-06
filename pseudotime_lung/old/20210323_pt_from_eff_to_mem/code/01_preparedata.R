library(reshape2)
library(TSCAN)
library(ggplot2)
library(RColorBrewer)
library(parallel)
library(here)
here()
rdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/pardoll/nsclc/pseudotime_lung/data/'
pdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/pardoll/nsclc/pseudotime_lung/plot/'
## ---------------------
## load and check data
## ---------------------
meta <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/pardoll/nsclc/data/ManaT2/blood.meta.rds')
saver <- readRDS('/home-4/whou10@jhu.edu/scratch/Wenpin/pardoll/nsclc/data/ManaT2/blood.saver.rds')
dm <- meta[, grepl('DC', colnames(meta))]
pt <- meta[,117]
names(pt) <- rownames(meta)

cellanno <- data.frame(cell = rownames(meta), sample = as.character(meta$tissue), stringsAsFactors = F)
rownames(cellanno) <- cellanno[,1]
saver <- saver[rowMeans(saver>0.01) > 0.01, ]
design <- matrix(rep(1, length(unique(cellanno[,2]))), nrow = length(unique(cellanno[,2])))
dimnames(design) <- list(unique(cellanno[,2]), 'intercept')
library(here)
here()
library(Seurat)
library(ggplot2)
library(RColorBrewer)
pdf(paste0(pdir, 'dm_tissuetime.pdf'), width = 4, height = 3)
ggplot(data = meta) + 
  geom_point(aes(x = DC1, y = DC2, color = tissue), size = 0.5) +
  theme_classic() 
dev.off()

## ---------------------
## construct MST 
## ---------------------
mc <- exprmclust(t(dm[,1:2]),reduce = T)

pdf(paste0(pdir, 'MST_cluster.pdf'), width = 4, height = 3)
plot(mc$MSTtree)
dev.off()
pdf(paste0(pdir, 'dm_cluster.pdf'), width = 4, height = 3)
ggplot(data=data.frame(d1=dm[names(mc$clusterid),1],d2=dm[names(mc$clusterid),2],cluster=as.character(mc$clusterid)),aes(x=d1,y=d2,col=cluster)) + geom_point(size=0.5) + theme_classic()+scale_color_brewer(palette = 'Set1') + xlab('DM1') + ylab('DM2')
dev.off()

## ---------------------
## construct pseudotime 
## ---------------------
ord <- TSCANorder(mc,orderonly = T,listbranch = T)
names(ord)
p1 <- ord[[1]]
p2 <- ord[[2]]

pdf(paste0(pdir, 'dm_pseudotime_path1.pdf'), width = 4, height = 3)
mycolor = c(colorRampPalette(rev(brewer.pal(11, 'RdYlBu')))(length(p1)), 'grey')
mycolor = c(rep(mycolor[1:30], each = 100), mycolor[31:length(mycolor)])
print(ggplot(data.frame(x = c(dm[p1,1], dm[!rownames(dm) %in% p1,1]), y = c(dm[p1,2], dm[!rownames(dm) %in% p1, 2]), time = c(seq(1, length(p1)), rep(NA, sum(!rownames(dm) %in% p1))))) + 
        geom_point(aes(x = x, y = y, col = time), size = 0.5) + 
        scale_color_gradientn(colors = mycolor) + 
        theme_classic() + xlab('DM1') + ylab('DM2'))
dev.off()

pdf(paste0(pdir, 'dm_pseudotime_path2.pdf'), width = 4, height = 3)
mycolor = c(colorRampPalette(rev(brewer.pal(11, 'RdYlBu')))(length(p2)), 'grey')
mycolor = c(rep(mycolor[1:30], each = 100), mycolor[31:length(mycolor)])
print(ggplot(data.frame(x = c(dm[p2,1], dm[!rownames(dm) %in% p2,1]), y = c(dm[p2,2], dm[!rownames(dm) %in% p2, 2]), time = c(seq(1, length(p2)), rep(NA, sum(!rownames(dm) %in% p2))))) + 
        geom_point(aes(x = x, y = y, col = time), size = 0.5) + 
        scale_color_gradientn(colors = mycolor) + 
        theme_classic() + xlab('DM1') + ylab('DM2'))
dev.off()


## ------------------------
## preparedata for testtime
## ------------------------
pt <- seq(1, length(p1))
names(pt) <- p1
saveRDS(pt, paste0(rdir, 'path1_pseudotime.rds'))
saveRDS(saver[, p1], paste0(rdir, 'path1_logsaver.rds'))
saveRDS(cellanno[p1,], paste0(rdir, 'path1_cellanno.rds'))
saveRDS(design, paste0(rdir, 'path1_design.rds'))


pt <- seq(1, length(p2))
names(pt) <- p2
saveRDS(pt, paste0(rdir, 'path2_pseudotime.rds'))
saveRDS(saver[, p2], paste0(rdir, 'path2_logsaver.rds'))
saveRDS(cellanno[p2,], paste0(rdir, 'path2_cellanno.rds'))
saveRDS(design, paste0(rdir, 'path2_design.rds'))


