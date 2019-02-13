library(huge)
#generate data

trialN = 10

cat("---------------------------------------------\n")
cat("--------------without scr--------------------\n")
cat("---------------------------------------------\n")
set.seed(123)
L = huge.generator(n = 100, d = 200, graph = "hub", g = 10)
t <- c()
for( i in 1:trialN){
  cat(".")
  t <- c(system.time(huge(L$data, method = "glasso", verbose = F))[1],t)
}
cat("\n")
cat("n=100\t d=200\t mean: ", mean(t), "\tvar: ", var(t), "\n")
set.seed(123)
L = huge.generator(n = 1000, d = 200, graph = "hub", g = 10)
t <- c()
for( i in 1:trialN){
  cat(".")
  t <- c(system.time(huge(L$data, method = "glasso", verbose = F))[1],t)
}
cat("\n")
cat("n=1000\t d=200\t mean: ", mean(t), "\tvar: ", var(t), "\n")
set.seed(123)
L = huge.generator(n = 100, d = 500, graph = "hub", g = 10)
t <- c()
for( i in 1:trialN){
  cat(".")
  t <- c(system.time(huge(L$data, method = "glasso", verbose = F))[1],t)
}
cat("\n")
cat("n=100\t d=500\t mean: ", mean(t), "\tvar: ", var(t), "\n")
set.seed(123)
L = huge.generator(n = 1000, d = 500, graph = "hub", g = 10)
t <- c()
for( i in 1:trialN){
  cat(".")
  t <- c(system.time(huge(L$data, method = "glasso", verbose = F))[1],t)
}
cat("\n")
cat("n=1000\t d=500\t mean: ", mean(t), "\tvar: ", var(t), "\n")

cat("---------------------------------------------\n")
cat("---------------with scr----------------------\n")
cat("---------------------------------------------\n")
set.seed(123)
L = huge.generator(n = 100, d = 200, graph = "hub", g = 10)
t <- c()
for( i in 1:trialN){
  cat(".")
  t <- c(system.time(huge(L$data, method = "glasso", verbose = F, scr = T))[1],t)
}
cat("\n")
cat("n=100\t d=200\t mean: ", mean(t), "\tvar: ", var(t), "\n")
set.seed(123)
L = huge.generator(n = 1000, d = 200, graph = "hub", g = 10)
t <- c()
for( i in 1:trialN){
  cat(".")
  t <- c(system.time(huge(L$data, method = "glasso", verbose = F, scr = T))[1],t)
}
cat("\n")
cat("n=1000\t d=200\t mean: ", mean(t), "\tvar: ", var(t), "\n")
set.seed(123)
L = huge.generator(n = 100, d = 500, graph = "hub", g = 10)
t <- c()
for( i in 1:trialN){
  cat(".")
  t <- c(system.time(huge(L$data, method = "glasso", verbose = F, scr = T))[1],t)
}
cat("\n")
cat("n=100\t d=500\t mean: ", mean(t), "\tvar: ", var(t), "\n")
set.seed(123)
L = huge.generator(n = 1000, d = 500, graph = "hub", g = 10)
t <- c()
for( i in 1:trialN){
  cat(".")
  t <- c(system.time(huge(L$data, method = "glasso", verbose = F, scr = T))[1],t)
}
cat("\n")
cat("n=1000\t d=500\t mean: ", mean(t), "\tvar: ", var(t), "\n")

