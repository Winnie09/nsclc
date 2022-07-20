sourcepath <- '/home/whou10/scratch16/whou10/trajectory_variability/package/Lamian/R/'
af <- list.files(sourcepath)
for (f in af){
  source(paste0(sourcepath, f))
}

input = readRDS('/home/whou10/scratch4/whou10/pardoll/nsclc/pseudotime_treg/mypt/data/C0_to_C7_lamian_input_R3_NR3.rds')

library(splines)
library(matrixcalc)
library(circlize)
Res <-
  cellPropTest(
    cellanno = input$cellanno,
    pseudotime = input$pseudotime,
    design = input$design,
    ncores = 1,
    test.type = 'Variable'
  )
rdir2 <- '/home/whou10/scratch4/whou10/pardoll/nsclc/pseudotime_treg/mypt/res_xcd/R3_NR3/pm/'
dir.create(rdir2, recursive = T, showWarnings = F)
saveRDS(Res, paste0(rdir2, 'res_tcd.rds'))

library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(reshape2)
library(scattermore)
library(viridis)
pdir <- '/home/whou10/scratch4/whou10/pardoll/nsclc/pseudotime_treg/mypt/plot_xcd/R3_NR3/pm/'
dir.create(pdir, recursive = T, showWarnings = F)
pdf(paste0(pdir, 'cell_density.pdf'), width = 3.2, height = 2.2)
plotGene(testobj = Res,
           gene = 'prop',
           variable = 'response',
           cellProp = T,
           x.lab = 'Pseudotime bin',
           y.lab = 'Cell density') +
  ggtitle('')
dev.off()

write.csv(Res[[1]], paste0(pdir, 'statistics.csv'))


pdf(paste0(pdir, 'cell_density_individual.pdf'), width = 2.9, height = 2.2)
plotGene(testobj = Res,
           gene = 'prop',
           variable = NULL,
           cellProp = T,
           x.lab = 'Pseudotime bin',
           y.lab = 'Cell density') +
  ggtitle('')
dev.off()

