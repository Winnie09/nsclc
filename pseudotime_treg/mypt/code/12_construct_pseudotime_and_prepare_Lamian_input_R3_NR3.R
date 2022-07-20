library(reshape2)
library(TSCAN)
library(ggplot2)
library(RColorBrewer)
library(parallel)
library(here)
here()
# ddir <- '/Users/wenpinhou/Dropbox/pardoll/nsclc/pseudotime_treg/data/'
# rdir <- '/Users/wenpinhou/Dropbox/pardoll/nsclc/pseudotime_treg/mypt/data/'
# pdir1 <- '/Users/wenpinhou/Dropbox/pardoll/nsclc/pseudotime_treg/mypt/plot/'
# pdir2 <- '/Users/wenpinhou/Dropbox/pardoll/nsclc/pseudotime_treg/mypt/plot/R3_NR3/'

ddir <- '/home/whou10/scratch4/whou10/pardoll/nsclc/pseudotime_treg/mypt/data/'
rdir <- '/home/whou10/scratch4/whou10/pardoll/nsclc/pseudotime_treg/mypt/data/'
pdir1 <- '/home/whou10/scratch4/whou10/pardoll/nsclc/pseudotime_treg/mypt/plot/'
pdir2 <- '/home/whou10/scratch4/whou10/pardoll/nsclc/pseudotime_treg/mypt/plot/R3_NR3/'


ser <- readRDS(paste0(ddir, 'treg.pseudo.rds'))
meta <- ser@meta.data
str(meta)
dim(meta)
## ---------------------
## load and check data
## ---------------------
# saver <- readRDS('/Users/wenpinhou/Dropbox/pardoll/nsclc/pseudotime_treg/data/saver_treg.rds')
saver = readRDS(paste0(ddir, 'saver_treg.rds'))
int <- intersect(colnames(saver), rownames(meta))
str(int)
meta = meta[colnames(saver), ]
dm <- meta[, grepl('DC', colnames(meta))]
dim(dm)
#####################
## plot current data
#####################
library(ggplot2)
library(RColorBrewer)
pdf(paste0(pdir1, 'dm_patient_alltreg.pdf'), width = 4, height = 3)
ggplot(data = meta) +
  geom_point(aes(x = DC1, y = DC2, color = patient_id), size = 0.5) +
  theme_classic() + xlab('DM1') + ylab('DM2')
dev.off()

pdf(paste0(pdir1, 'dm_patient_facet_alltreg.pdf'), width = 6, height = 4.5)
ggplot(data = meta) +
  geom_point(aes(x = DC1, y = DC2, color = patient_id), size = 0.1) +
  theme_classic() + xlab('DM1') + ylab('DM2') + facet_wrap(~patient_id)
dev.off()

pdf(paste0(pdir1, 'dm_celltype_alltreg.pdf'), width = 4.8, height = 3)
ggplot(data = meta) +
  geom_point(aes(x = DC1, y = DC2, color = celltype), size = 0.5) +
  theme_classic() + xlab('DM1') + ylab('DM2')
dev.off()

## ---------------------------------
## filter samples with too few cells
## ---------------------------------
## select two cell types (C0, C7) for pseudotime construction
ct = meta[, 'celltype']
names(ct) = rownames(meta)
ct = ct[ct %in% c('C0-Activated (1)', 'C7-Th1-like/cytoxic')]
cell.select = names(ct)
str(cell.select)

pdf(paste0(pdir1, 'dm_patient_c0_c7_facet.pdf'), width = 6, height = 4.5)
ggplot(data = meta[cell.select, ]) +
  geom_point(aes(x = DC1, y = DC2, color = patient_id), size = 0.5) +
  theme_classic() + facet_wrap(~patient_id) + xlab('DM1') + ylab('DM2')
dev.off()



## filtering out samples with <10 cells, finalizing cell.select
cellanno <- data.frame(cell = rownames(meta), patient = as.character(meta$patient_id), stringsAsFactors = F)
dim(cellanno)
rownames(cellanno) = rownames(meta)

cellanno = cellanno[cell.select,]
tab = table(cellanno[,2])
s.select = names(tab[tab >= 20])
str(s.select)
cell.select = rownames(cellanno[cellanno[,2]%in% s.select, ])
str(cell.select)  

## refresh data's cells
cellanno = cellanno[cell.select,]
str(cellanno)
ct = ct[cell.select]
expr = saver[, cell.select]
dim(expr)

expr <- expr[rowMeans(expr>0.01) > 0.01, ]
design <- matrix(rep(1, length(unique(cellanno[,2]))), nrow = length(unique(cellanno[,2])))
dimnames(design) <- list(unique(cellanno[,2]), 'intercept')
design

tmp = meta[, c('response', 'patient_id')]  
tmp = tmp[!duplicated(tmp),] 
str(tmp)

pdf(paste0(pdir2, 'dm_c0_c7_selected_patient_facet.pdf'), width = 6, height = 4.5)
ggplot(data = meta[cell.select, ]) +
  geom_point(aes(x = DC1, y = DC2, color = patient_id), size = 0.5) +
  theme_classic() + facet_wrap(~patient_id) + xlab('DM1') + ylab('DM2')
dev.off()

## ---------------------
## construct MST 
## ---------------------
ct.clu = rep(1, length(ct))
names(ct.clu) = names(ct)
ct.clu[ct == 'C7-Th1-like/cytoxic'] <- 2
table(ct.clu)

mc <- exprmclust(t(dm[cell.select, 1:2]),reduce = F, cluster = ct.clu)
pdf(paste0(pdir2, 'MST_cluster.pdf'), width = 4, height = 3)
plot(mc$MSTtree)
dev.off()

pdf(paste0(pdir2, 'dm_cluster.pdf'), width = 3, height = 2)
ggplot(data=data.frame(d1=dm[names(mc$clusterid),1],d2=dm[names(mc$clusterid),2],cluster=as.character(mc$clusterid)),aes(x=d1,y=d2,col=cluster)) + geom_point(size=0.2) + theme_classic()+scale_color_brewer(palette = 'Set1') + xlab('DM1') + ylab('DM2') +scale_color_brewer(palette = 'Set1')
dev.off()

## ---------------------
## construct pseudotime 
## ---------------------
ord <- TSCANorder(mc,orderonly = T,listbranch = F, startcluster = 1)
# # names(ord)
# #  $ backbone 6,8,2,1: chr [1:3503] 
# #  $ branch: 6,4,9,5 : chr [1:433] 
# #  $ branch: 6,4,3,7 :
# p1 <- ord[[1]] ## right
# p2 <- ord[[3]]  ## lower-left
# p3 <- ord[[2]] ## upper-left
p1 = ord
pt <- seq(1, length(p1))
names(pt) <- p1

pdf(paste0(pdir2, 'dm_pseudotime_path1.pdf'), width = 3, height = 2)
mycolor = c(colorRampPalette(rev(brewer.pal(11, 'RdYlBu')))(length(p1)), 'grey')
print(ggplot(data.frame(x = c(dm[p1,1], dm[!rownames(dm) %in% p1,1]), y = c(dm[p1,2], dm[!rownames(dm) %in% p1, 2]), time = c(seq(1, length(p1)), rep(NA, sum(!rownames(dm) %in% p1))))) + 
        geom_point(aes(x = x, y = y, col = time), size = 0.1) + 
        scale_color_gradientn(colors = mycolor) + 
        theme_classic() + xlab('DM1') + ylab('DM2'))
dev.off()

## ----------------------------------------
## preparedata for Lamian TDE and XDE input
## ----------------------------------------
r = rep(0, nrow(design))
names(r) = rownames(design)
r[tmp[tmp[,1]=='R',2]] <- 1
design = cbind(design, response = r[rownames(design)])
input = list(expr = expr[,p1], cellanno = cellanno[p1,], design = design, pseudotime = pt)
saveRDS(input, paste0(rdir, 'C0_to_C7_lamian_input_R3_NR3.rds'))


