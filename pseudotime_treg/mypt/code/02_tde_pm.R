source('/home/whou10/scratch16/whou10/trajectory_variability/function/01_function.R')
input = readRDS('/home/whou10/scratch4/whou10/pardoll/nsclc/pseudotime_treg/mypt/data/C0_to_C7_lamian_input.rds')
Res <- testpt(expr = input$expr, cellanno = input$cellanno, pseudotime = input$pseudotime, design=input$design,  test.type='Time', test.method = 'permutation', ncores = 1) 
saveRDS(Res, '/home/whou10/scratch4/whou10/pardoll/nsclc/pseudotime_treg/mypt/res_tde/pm/res_tde.rds')
