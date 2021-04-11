library(here)
setwd(here())
ddir <- 'pseudotime_lung/data/'
rdir <- 'pseudotime_lung/result/path2/'
dir.create(rdir, showWarnings = F, recursive = T)
pt <-  readRDS(paste0(ddir, 'path2_pseudotime.rds'))
expr <- readRDS(paste0(ddir, 'path2_logsaver.rds'))
cellanno <- readRDS(paste0(ddir, 'path2_cellanno.rds'))
design <- readRDS(paste0(ddir, 'path2_design.rds'))

expr <- expr[rowMeans(expr > 0.01) > 0.01, , drop=FALSE]
source('/home-4/whou10@jhu.edu/scratch/Wenpin/trajectory_variability/function/01_function.R')
res <- testpt(expr, cellanno = cellanno, pseudotime = pt, design=design, ncores=24, test.type='Time')
saveRDS(res, paste0(rdir, '/testtime_res.rds'))
