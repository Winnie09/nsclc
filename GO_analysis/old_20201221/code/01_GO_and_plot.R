library(ggplot2)
library(RColorBrewer)
## GO analysis
ddir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/pardoll/nsclc/data/TF/'
rdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/pardoll/nsclc/GO_analysis/result/'
pdir <- '/home-4/whou10@jhu.edu/scratch/Wenpin/pardoll/nsclc/GO_analysis/plot/'

af <- list.files(ddir, pattern = 'csv')
af <- af[!grepl('all', af)]
ser = readRDS(paste0(ddir, 'ser.antigen.rds'))
gb <- rownames(ser@assays$RNA@counts)
source('/home-4/whou10@jhu.edu/scratch/Wenpin/resource/myfunc/01_function.R')
  
for (f in af){
  print(f)
  res <- read.csv(paste0(ddir, f), as.is = TRUE)
  res.bak = res
  res <- res.bak[res.bak[,2]>=0, ]
  gs <- res[res[,5] < 0.05, 6]
  res.go <- myGO(gs, gb)
  res.go <- res.go[res.go$FDR<0.05, ]
  res.go <- res.go[order(res.go$FDR, -res.go$FC), ]
  print(res.go[1:min(20, nrow(res.go)),c(2,7,8)])
  write.csv(res.go, paste0(rdir, 'GO_nonFilterGene_', strsplit(f, '\\.')[[1]][1], '_side_in_', f))
  
  pdf(paste0(pdir, 'GO_nonFilterGene_', strsplit(f, '\\.')[[1]][1], '_side_in_', sub('.csv','.pdf',f)), width = 10, height = 7)
  print(ggplot(data = res.go[1:20,]) +
    geom_bar(aes(x = Term, y = FC, fill = Term), stat="identity", width=0.5) +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, color= 'black', hjust = 1),
          legend.position = 'none')+
    scale_fill_manual(values = colorRampPalette(brewer.pal(8,'Set1'))(20)))
  dev.off()
    
  
  res <- res.bak[res.bak[,2]<0, ]
  gs <- res[res[,5] < 0.05, 6]
  res.go <- myGO(gs, gb)
  res.go <- res.go[res.go$FDR<0.05, ]
  res.go <- res.go[order(res.go$FDR, -res.go$FC), ]
  print(res.go[1:min(20, nrow(res.go)),c(2,7,8)])
  write.csv(res.go, paste0(rdir, 'GO_nonFilterGene_', strsplit(f, '\\.')[[1]][2], '_side_in_', f))
  
  pdf(paste0(pdir, 'GO_nonFilterGene_', strsplit(f, '\\.')[[1]][2], '_side_in_', sub('.csv','.pdf',f)), width = 10, height = 7)
  print(ggplot(data = res.go[1:20,]) +
    geom_bar(aes(x = Term, y = FC, fill = Term), stat="identity", width=0.5) +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, color= 'black', hjust = 1),
          legend.position = 'none')+
    scale_fill_manual(values = colorRampPalette(brewer.pal(8,'Set1'))(20)))
  dev.off()
}

af <- list.files(ddir, pattern = 'csv')
af <- af[grepl('all', af)]
for (f in af){
  print(f)
  res <- read.csv(paste0(ddir, f), as.is = TRUE)
  gs <- res[res[,5] < 0.05, 6]
  res.go <- myGO(gs, res[,6])
  res.go <- res.go[res.go$FDR<0.05, ]
  res.go <- res.go[order(res.go$FDR, -res.go$FC), ]
  print(res.go[1:min(20, nrow(res.go)),c(2,7,8)])
  write.csv(res.go, paste0(rdir, 'GO_basedOnFilteredGene_', f))
}

