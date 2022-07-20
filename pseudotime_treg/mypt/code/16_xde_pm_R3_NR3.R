source('/home/whou10/scratch16/whou10/trajectory_variability/function/01_function.R')
input = readRDS('/home/whou10/scratch4/whou10/pardoll/nsclc/pseudotime_treg/mypt/data/C0_to_C7_lamian_input_R3_NR3.rds')
Res <- testpt(expr = input$expr, cellanno = input$cellanno, pseudotime = input$pseudotime, design=input$design,  test.type='VARIABLE', test.method = 'permutation', ncores = 10) 
dir.create('/home/whou10/scratch4/whou10/pardoll/nsclc/pseudotime_treg/mypt/res_xde/pm', recursive = T)
saveRDS(Res, '/home/whou10/scratch4/whou10/pardoll/nsclc/pseudotime_treg/mypt/res_xde/R3_NR3/pm/res_xde_pm_R3_NR3.rds')

