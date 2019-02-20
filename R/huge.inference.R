#-----------------------------------------------------------------------#
# Package: High-dimensional Undirected Graph Estimation                 #
# huge.inference(): graph inference                                     #
#-----------------------------------------------------------------------#

#' Graph inference
#'
#' Implements the inference for high dimensional graphical models, including Gaussian and Nonparanormal graphical models
#' We consider the problems of testing the presence of a single edge and the hypothesis is that the edge is absent.
#'
#' For Nonparanormal graphical model we provide Score test method and Wald Test. However it is really slow for inferencing on Nonparanormal model, especially for large data.
#'
#'
#' @param data The input \code{n} by \code{d} data matrix(\code{n} is the sample size and \code{d} is the dimension).
#' @param T The estimated inverse of correlation matrix of the data.
#' @param adj The adjacency matrix corresponding to the graph.
#' @param alpha The significance level of hypothesis.The default value is \code{0.05}.
#' @param type The type of input data. There are 2 options: \code{"Gaussian"} and \code{"Nonparanormal"}. The default value is \code{"Gaussian"}.
#' @param method When using nonparanormal graphical model. Test method with 2 options: \code{"score"} and \code{"wald"}. The default value is \code{"score"}.
#' @seealso \code{\link{huge}}, and \code{\link{huge-package}}.
#' @return
#' An object is returned:
#' \item{data}{
#'   The \code{n} by \code{d} data matrix from the input.
#' }
#' \item{p}{
#'   The \code{d} by \code{d} p-value matrix of hypothesis.
#' }
#' \item{error}{
#'   The type I error of hypothesis at alpha significance level.
#' }
#' @examples
#' #generate data
#' L = huge.generator(n = 50, d = 12, graph = "hub", g = 4)
#'
#' #graph path estimation using glasso
#' est = huge(L$data, method = "glasso")
#'
#' #inference of Gaussian graphical model at 0.05 significance level
#' T = tail(est$icov, 1)[[1]]
#' out1 = huge.inference(L$data, T, L$theta)
#'
#' #inference of Nonparanormal graphical model using score test at 0.05 significance level
#' T = tail(est$icov, 1)[[1]]
#' out2 = huge.inference(L$data, T, L$theta, type = "Nonparanormal")
#'
#' #inference of Nonparanormal graphical model using wald test at 0.05 significance level
#' T = tail(est$icov, 1)[[1]]
#' out3 = huge.inference(L$data, T, L$theta, type = "Nonparanormal", method = "wald")
#'
#' #inference of Nonparanormal graphical model using wald test at 0.1 significance level
#' T = tail(est$icov, 1)[[1]]
#' out4 = huge.inference(L$data, T, L$theta, 0.1, type = "Nonparanormal", method = "wald")
#' @references
#' 1.Q Gu, Y Cao, Y Ning, H Liu. Local and global inference for high dimensional nonparanormal graphical models.\cr
#' 2.J Jankova, S Van De Geer. Confidence intervals for high-dimensional inverse covariance estimation. \emph{Electronic Journal of Statistics}, 2015.\cr
#' @export
huge.inference = function(data, T, adj, alpha = 0.05, type = "Gaussian", method = "score"){
  d = ncol(data)
  n = nrow(data)

  if(type == "Gaussian")
  {
    x=scale(data)
    U=cor(x)
    W<-2*T - T%*%U%*%T
    W_n<-matrix(0,d,d)
    sigma<-matrix(0,d,d)
    for(j in 1:d)
    {
      for (k in 1:d)
      {
        sigma[j, k] = sqrt(T[j, j]*T[k, k]+T[j, k]^2)
        W_n[j, k] = sqrt(n)*W[j, k]/sigma[j, k]
      }
    }
    p<-matrix(0,d,d)
    p<-apply(W_n, 1, function(x) 2*(1 - pnorm(abs(x))))
    rm(sigma,W_n)
  }
  if(type == "Nonparanormal")
  {
    x=data
    U<-matrix(0,d,d)
    G=list()
    Test<-matrix(0,d,d)
    for(i in 1:n)
      G[[i]]<-matrix(0,d,d)
    Temp_jk<-matrix(0,n,n)
    for(j in 1:d)
    {
      for(k in 1:d)
      {
        if(j==k)
        {
          U[j, k] = 1
          next
        }

        for(i1 in 1:n)
        {
          for(i2 in 1:n)
          {
            Temp_jk[i1, i2] = sign((x[i1, j] - x[i2, j])*(x[i1, k] - x[i2, k]))
            G[[i1]][j, k] = G[[i1]][j, k] - pi/2*Temp_jk[i1, i2]

          }
        }
        U[j, k] = sin(pi/2*sum(Temp_jk)/((n-1)*n))
        for(i in 1:n)
          G[[i]][j, k] = G[[i]][j, k]/(n-1) + asin(U[j, k])
      }
    }
    #F
    F<-apply(U,1,function(x) sqrt(1-x^2))
    #R
    R<-matrix(0,d^2,d^2)
    for (i in 1:n)
      R<-R + as.matrix(as.vector(F*G[[i]]))%*%as.vector(F*G[[i]])
    R<-R/n

    #kronecker product of T
    T_k<-kronecker(T, T)

    if(method == "score")
    {
      S<-matrix(0,d,d)
      sigma<-matrix(0,d,d)
      #ST_n
      ST_n<-matrix(0,d,d)
      for(j in 1:d)
      {
        for(k in 1:d)
        {
          #S
          ej<-matrix(0,d,1)
          ek<-matrix(0,d,1)
          ej[j] = 1
          ek[k] = 1
          T_h<-T
          T_h[j, k]=0
          S[j, k] = t(ej)%*%t(T_h)%*%U%*%T_h%*%ek/(T[j, j]*T[k, k])

          #w
          temp1<-T_k[,j*k]
          temp1<-temp1[-j*k]
          w<-as.matrix((-temp1)/T_k[j*k, j*k])

          #sigma
          temp2<-R[j*k,]
          temp2<-as.matrix(temp2[-j*k])
          sigma[j, k] = sqrt(R[j*k, j*k] - 2*t(temp2)%*%w + t(w)%*%R[-j*k, -j*k]%*%w)

          ST_n[j, k] = S[j, k]*sqrt(n)/(2*sigma[j, k])

        }
      }

      #p-value
      p<-matrix(0,d,d)
      p<-apply(ST_n, 1, function(x) 2*(1 - pnorm(abs(x))))
      rm(temp1,temp2,ST_n,S)
    }

    if(method == "wald")
    {
      T_W<-matrix(0,d,d)
      sigma<-matrix(0,d,d)
      #W_n
      W_n<-matrix(0,d,d)
      temp1<-T%*%U
      temp2<-U%*%T
      for(j in 1:d)
      {
        for(k in 1:d)
        {
          #w
          temp3<-T_k[,j*k]
          temp3<-temp3[-j*k]
          w<-as.matrix((-temp3)/T_k[j*k,j*k])

          #sigma
          temp4<-R[j*k,]
          temp4<-as.matrix(temp4[-j*k])
          sigma[j, k] = sqrt(R[j*k, j*k] - 2*t(temp4)%*%w + t(w)%*%R[-j*k, -j*k]%*%w)

          #T_W
          ej<-matrix(0,d,1)
          ek<-matrix(0,d,1)
          ej[j] = 1
          ek[k] = 1
          T_W[j ,k] = (T[j,k]*t(ej)%*%temp1%*%ek + T[j, k]*t(ej)%*%temp2%*%ek - t(ej)%*%t(T)%*%temp2%*%ek)/(t(ej)%*%temp1%*%ek + t(ej)%*%temp2%*%ek - 1)
          W_n[j, k] = T_W[j, k]*sqrt(n)/(2*sigma[j, k]*T[j, j]*T[k, k])

        }
      }


      #p-value
      p<-matrix(0,d,d)
      p<-apply(W_n, 1, function(x) 2*(1 - pnorm(abs(x))))
      rm(temp1,temp2,temp3,temp4,T_W,W_n)
    }
    rm(R,F,G,w,T_k)
  }

  error=0
  for(j in 1:d)
  {
    for(k in 1:d)
    {
      if(j==k)
        next
      if(p[j, k]<alpha && adj[j, k]==0){
        error=error+1
      }
    }
  }
  error=error/d^2


  inf = list()
  inf$data = data
  inf$p = p
  inf$error = error

  rm(U,p)
  return(inf)
}
