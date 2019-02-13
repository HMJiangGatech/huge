library(huge)
#generate data

cat("---------------------------------------------\n")
cat("--------------without scr--------------------\n")
cat("---------------------------------------------\n")
set.seed(123)
L = huge.generator(n = 100, d = 200, graph = "hub", g = 10)
out_glasso = huge(L$data, method = "mb", verbose = F)
out_glasso = huge(L$data, method = "glasso", verbose = F)

cat("---------------------------------------------\n")
cat("---------------with scr----------------------\n")
cat("---------------------------------------------\n")
huge(L$data, method = "mb", verbose = F, scr = T)
huge(L$data, method = "glasso", verbose = F, scr = T)

cat("---------------------------------------------\n")
cat("---------------use tiger---------------------\n")
cat("---------------------------------------------\n")
huge(L$data, method = "tiger", verbose = FALSE)

cat("---------------------------------------------\n")
cat("---------------inference---------------------\n")
cat("---------------------------------------------\n")
huge.inference(L$data, tail(out_glasso$icov, 1)[[1]], L$theta)
