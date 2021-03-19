library(here)
here()
library(Seurat)
dfs <- readRDS(here('pseudotime_lung', 'data', 'diffusion.antigen.rds'))
dfs_meta <- readRDS(here('pseudotime_lung', 'data', 'diffusion.coordinate.meta.rds'))
    
dfs_meta$resi_tumor = as.factor(dfs_meta$resi_tumor)
library(ggplot2)
library(RColorBrewer)
ggplot(data = dfs_meta) + 
  geom_point(aes(x = DC1, y = DC2, color = sample), size = 0.3, alpha = 0.8) +
  theme_classic() 
expr = as.matrix(dfs@assays$RNA@counts)

cellanno = data.frame(cell = rownames(dfs_meta),
                      sample = dfs_meta$patient,
                      stringsAsFactors = FALSE)
data = dfs_meta[dfs_meta$DC1 > 0.003, ]
table(data$patient)

