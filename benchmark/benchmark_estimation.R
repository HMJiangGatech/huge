library(huge)
trialN = 5
n = 200

cat("---------------------------------------------\n")
cat("-------------Graphical Estimation------------\n")
set.seed(123)
L = huge.generator(n = 200, d = 200, graph = "hub", g = 4)
S = cor(scale(L$data))
t_huge_glasso <- c()
out_huge_glasso <- c()
objectvalue_huge_glasso <- c()
for( i in 1:trialN){
  cat(".")
  t_huge_glasso <- c(system.time(out_huge_glasso[[i]] <- huge(L$data, method = "glasso", verbose = FALSE))[1],t_huge_glasso)
  T <- tail(out_huge_glasso[[i]]$icov,1)[[1]]
  lambda <- tail(out_huge_glasso[[i]]$lambda,1)
  objectvalue_huge_glasso <- c(log(det(T))-sum(diag(S%*%T)) - lambda*norm(T,"1"), objectvalue_huge_glasso)
}
cat("\n")
cat("n=200\t d=200\t mean time of huge: ", mean(t_huge_glasso))
cat("n=200\t d=200\t mean object value of huge glasso: ", mean(objectvalue_huge_glasso))

cat("---------------------------------------------\n")

t_huge_tiger <- c()
out_huge_tiger <- c()
objectvalue_huge_tiger <- c()
for( i in 1:trialN){
  cat(".")
  t_huge_tiger <- c(system.time(out_huge_tiger[[i]] <- huge(L$data, method = "tiger", verbose = FALSE))[1],t_huge_tiger)
  T <- tail(out_huge_tiger[[i]]$icov,1)[[1]]
  lambda <- tail(out_huge_tiger[[i]]$lambda,1)
  objectvalue_huge_tiger <- c(log(det(T))-sum(diag(S%*%T)) - lambda*norm(T,"1"), objectvalue_huge_tiger)
}
cat("\n")
cat("n=200\t d=200\t mean time of huge: ", mean(t_huge_tiger))
cat("n=200\t d=200\t mean object value of huge tiger: ", mean(objectvalue_huge_tiger))

cat("---------------------------------------------\n")
library(QUIC)
t_quic <- c()
out_quic <- c()
objectvalue_quic <- c()
for(i in 1:trialN){
  cat(".")
  t_quic <- c(system.time(out_quic[[i]] <- QUIC(S, out_huge_glasso[[i]]$lambda[[1]], out_huge_glasso[[i]]$lambda))[1],t_quic)
  T <- out_quic[[i]]$X[,,dim(out_quic[[i]]$X)[[3]]]
  lambda <- tail(out_huge_glasso[[i]]$lambda,1)
  objectvalue_quic <- c(log(det(T))-sum(diag(S%*%T)) - lambda*norm(T,"1"), objectvalue_quic)
}
cat("\n")
cat("n=200\t d=200\t mean time of QUIC: ", mean(t_quic))
cat("n=200\t d=200\t mean object value of QUIC: ", mean(objectvalue_quic))

cat("---------------------------------------------\n")
library(clime)
t_clime <- c()
out_clime <- c()
objectvalue_clime <- c()
for(i in 1:trialN){
  cat(".")
  t_clime <- c(system.time(out_clime[[i]] <- clime(L$data, out_huge_glasso[[i]]$lambda))[1],t_clime)
  T <- tail(out_clime[[i]]$Omegalist,1)[[1]]
  lambda <- tail(out_huge_glasso[[i]]$lambda,1)
  objectvalue_clime <- c(log(det(T))-sum(diag(S%*%T)) - lambda*norm(T,"1"), objectvalue_clime)
}
cat("\n")
cat("n=200\t d=200\t mean time of clime: ", mean(t_clime))
cat("n=200\t d=200\t mean object value of clime: ", mean(objectvalue_clime))


