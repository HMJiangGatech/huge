library(huge)
trialN = 5
n = 200

cat("---------------------------------------------\n")
cat("-------------Graphical Inference-------------\n")
set.seed(123)
L1 = huge.generator(n = 200, d = 200, graph = "band")
error_b_1 <- c()
error_b_2 <- c()
for(i in 1:trialN){
  cat(".")
  out_glasso = huge(L1$data, method = "glasso", verbose = FALSE)
  error_b_1 <- c(huge.inference(L1$data, tail(out_glasso$icov, 1)[[1]], L1$theta)$error,error_b_1)
  error_b_2 <- c(huge.inference(L1$data, tail(out_glasso$icov, 1)[[1]], L1$theta, 0.1)$error,error_b_2)
}
cat("\n")
cat("n=200\t d=200\t alpha=0.05\t mean type I error of band: ", mean(error_b_1))
cat("n=200\t d=200\t alpha=0.1\t mean type I error of band: ", mean(error_b_2))
cat("---------------------------------------------\n")

set.seed(123)
L2 = huge.generator(n = 200, d = 200, graph = "scale-free")
error_s_1 <- c()
error_s_2 <- c()
for(i in 1:trialN){
  cat(".")
  out_glasso = huge(L2$data, method = "glasso", verbose = FALSE)
  error_s_1 <- c(huge.inference(L2$data, tail(out_glasso$icov, 1)[[1]], L2$theta)$error,error_s_1)
  error_s_2 <- c(huge.inference(L2$data, tail(out_glasso$icov, 1)[[1]], L2$theta, 0.1)$error,error_s_2)
}
cat("\n")
cat("n=200\t d=200\t alpha=0.05\t mean type I error of scale-free: ", mean(error_s_1))
cat("n=200\t d=200\t alpha=0.1\t mean type I error of scale-free: ", mean(error_s_2))
cat("---------------------------------------------\n")

set.seed(123)
L3 = huge.generator(n = 200, d = 200, graph = "hub")
error_h_1 <- c()
error_h_2 <- c()
for(i in 1:trialN){
  cat(".")
  out_glasso = huge(L3$data, method = "glasso", verbose = FALSE)
  error_h_1 <- c(huge.inference(L3$data, tail(out_glasso$icov, 1)[[1]], L3$theta)$error,error_h_1)
  error_h_2 <- c(huge.inference(L3$data, tail(out_glasso$icov, 1)[[1]], L3$theta, 0.1)$error,error_h_2)
}
cat("\n")
cat("n=200\t d=200\t alpha=0.05\t mean type I error of hub: ", mean(error_h_1))
cat("n=200\t d=200\t alpha=0.1\t mean type I error of hub: ", mean(error_h_2))
cat("---------------------------------------------\n")