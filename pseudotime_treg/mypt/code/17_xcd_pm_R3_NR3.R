sourcepath <- '/home/whou10/scratch16/whou10/trajectory_variability/package/Lamian/R/'
af <- list.files(sourcepath)
for (f in af){
  source(paste0(sourcepath, f))
}

input = readRDS('/home/whou10/scratch4/whou10/pardoll/nsclc/pseudotime_treg/mypt/data/C0_to_C7_lamian_input_R3_NR3.rds')

Res <-
  cellPropTest(
    cellanno = input$cellanno,
    pseudotime = input$pseudotime,
    design = input$design,
    ncores = 1,
    test.type = 'Variable'
  )
rdir2 <- '/home/whou10/scratch4/whou10/pardoll/nsclc/pseudotime_treg/mypt/res_xcd/N3_NR3/'
dir.create(rdir2, recursive = T, showWarnings = F)
saveRDS(Res, paste0(rdir2, 'res_tcd.rds'))
