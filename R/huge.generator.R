#-----------------------------------------------------------------------#
# Package: High-dimensional Undirected Graph Estimation                 #
# huge.generator(): Data generator                                      #
# Authors: Tuo Zhao and Han Liu                                         #
# Emails: <tzhao5@jhu.edu> and <hanliu@cs.jhu.edu>                      #
# Date: July 15th 2011                                                  #
# Version: 1.1.0                                                        #
#-----------------------------------------------------------------------#

#' Data generator
#'
#' Implements the data generation from multivariate normal distributions with different graph structures, including \code{"random"}, \code{"hub"}, \code{"cluster"}, \code{"band"} and \code{"scale-free"}.
#'
#' @param n The number of observations (sample size). The default value is \code{200}.
#' @param d The number of variables (dimension). The default value is \code{50}.
#' @param graph The graph structure with 4 options: \code{"random"}, \code{"hub"}, \code{"cluster"}, \code{"band"} and \code{"scale-free"}.
#' @param v The off-diagonal elements of the precision matrix, controlling the magnitude of partial correlations with \code{u}. The default value is \code{0.3}.
#' @param u A positive number being added to the diagonal elements of the precision matrix, to control the magnitude of partial correlations. The default value is \code{0.1}.
#' @param g For \code{"cluster"} or \code{"hub"} graph, \code{g} is the number of hubs or clusters in the graph. The default value is about \code{d/20} if \code{d >= 40} and \code{2} if \code{d < 40}. For \code{"band"} graph, \code{g} is the bandwidth and the default value is \code{1}. NOT applicable to \code{"random"} graph.
#' @param prob For \code{"random"} graph, it is the probability that a pair of nodes has an edge. The default value is \code{3/d}. For \code{"cluster"} graph, it is the probability that a pair of nodes has an edge in each cluster. The default value is \code{6*g/d} if \code{d/g <= 30} and \code{0.3} if \code{d/g > 30}. NOT applicable to \code{"hub"} or \code{"band"} graphs.
#' @param vis Visualize the adjacency matrix of the true graph structure, the graph pattern, the covariance matrix and the empirical covariance matrix. The default value is \code{FALSE}
#' @param verbose If \code{verbose = FALSE}, tracing information printing is disabled. The default value is \code{TRUE}.
#' @details
#' Given the adjacency matrix \code{theta}, the graph patterns are generated as below:\cr\cr
#' (I) \code{"random"}: Each pair of off-diagonal elements are randomly set \code{theta[i,j]=theta[j,i]=1} for \code{i!=j} with probability \code{prob}, and \code{0} other wise. It results in about \code{d*(d-1)*prob/2} edges in the graph.\cr\cr
#' (II)\code{"hub"}:The row/columns are evenly partitioned into \code{g} disjoint groups. Each group is associated with a "center" row \code{i} in that group. Each pair of off-diagonal elements are set \code{theta[i,j]=theta[j,i]=1} for \code{i!=j} if \code{j} also belongs to the same group as \code{i} and \code{0} otherwise. It results in \code{d - g} edges in the graph.\cr\cr
#' (III)\code{"cluster"}:The row/columns are evenly partitioned into \code{g} disjoint groups. Each pair of off-diagonal elements are set \code{theta[i,j]=theta[j,i]=1} for \code{i!=j} with the probability \code{prob}if both \code{i} and \code{j} belong to the same group, and \code{0} other wise. It results in about \code{g*(d/g)*(d/g-1)*prob/2} edges in the graph.\cr\cr
#' (IV)\code{"band"}: The off-diagonal elements are set to be \code{theta[i,j]=1} if \code{1<=|i-j|<=g} and \code{0} other wise. It results in \code{(2d-1-g)*g/2} edges in the graph.\cr\cr
#' (V) \code{"scale-free"}: The graph is generated using B-A algorithm. The initial graph has two connected nodes and each new node is connected to only one node in the existing graph with the probability proportional to the degree of the each node in the existing graph. It results in \code{d} edges in the graph.
#'
#' The adjacency matrix \code{theta} has all diagonal elements equal to \code{0}. To obtain a positive definite precision matrix, the smallest eigenvalue of \code{theta*v} (denoted by \code{e}) is computed. Then we set the precision matrix equal to \code{theta*v+(|e|+0.1+u)I}. The covariance matrix is then computed to generate multivariate normal data.
#' @return
#' An object with S3 class "sim" is returned:
#' \item{data}{
#'   The \code{n} by \code{d} matrix for the generated data
#' }
#' \item{sigma}{
#'   The covariance matrix for the generated data
#' }
#' \item{omega}{
#'   The precision matrix for the generated data
#' }
#' \item{sigmahat}{
#'   The empirical covariance matrix for the generated data
#' }
#' \item{theta}{
#'   The adjacency matrix of true graph structure (in sparse matrix representation) for the generated data
#' }
#' @seealso \code{\link{huge}} and \code{\link{huge-package}}
#' @examples
#' ## band graph with bandwidth 3
#' L = huge.generator(graph = "band", g = 3)
#' plot(L)
#'
#' ## random sparse graph
#' L = huge.generator(vis = TRUE)
#'
#' ## random dense graph
#' L = huge.generator(prob = 0.5, vis = TRUE)
#'
#' ## hub graph with 6 hubs
#' L = huge.generator(graph = "hub", g = 6, vis = TRUE)
#'
#' ## hub graph with 8 clusters
#' L = huge.generator(graph = "cluster", g = 8, vis = TRUE)
#'
#' ## scale-free graphs
#' L = huge.generator(graph="scale-free", vis = TRUE)
#' @export
huge.generator = function(n = 200, d = 50, graph = "random", v = NULL, u = NULL, g = NULL, prob = NULL, vis = FALSE, verbose = TRUE){
  gcinfo(FALSE)
  if(verbose) cat("Generating data from the multivariate normal distribution with the", graph,"graph structure....")
  if(is.null(g)){
    g = 1
    if(graph == "hub" || graph == "cluster"){
      if(d > 40)  g = ceiling(d/20)
      if(d <= 40) g = 2
    }
  }

  if(graph == "random"){
    if(is.null(prob))  prob = min(1, 3/d)
    prob = sqrt(prob/2)*(prob<0.5)+(1-sqrt(0.5-0.5*prob))*(prob>=0.5)
  }

  if(graph == "cluster"){
    if(is.null(prob)){
      if(d/g > 30)  prob = 0.3
      if(d/g <= 30)  prob = min(1,6*g/d)
    }
    prob = sqrt(prob/2)*(prob<0.5)+(1-sqrt(0.5-0.5*prob))*(prob>=0.5)
  }


  # parition variables into groups
  g.large = d%%g
  g.small = g - g.large
  n.small = floor(d/g)
  n.large = n.small+1
  g.list = c(rep(n.small,g.small),rep(n.large,g.large))
  g.ind = rep(c(1:g),g.list)
  rm(g.large,g.small,n.small,n.large,g.list)
  gc()

  # build the graph structure
  theta = matrix(0,d,d);
  if(graph == "band"){
    if(is.null(u)) u = 0.1
    if(is.null(v)) v = 0.3
    for(i in 1:g){
      diag(theta[1:(d-i),(1+i):d]) = 1
      diag(theta[(1+i):d,1:(d-1)]) = 1
    }
  }
  if(graph == "cluster"){
    if(is.null(u)) u = 0.1
    if(is.null(v)) v = 0.3
    for(i in 1:g){
       tmp = which(g.ind==i)
       tmp2 = matrix(runif(length(tmp)^2,0,0.5),length(tmp),length(tmp))
       tmp2 = tmp2 + t(tmp2)
       theta[tmp,tmp][tmp2<prob] = 1
       rm(tmp,tmp2)
       gc()
    }
  }
  if(graph == "hub"){
  if(is.null(u)) u = 0.1
  if(is.null(v)) v = 0.3
  for(i in 1:g){
     tmp = which(g.ind==i)
     theta[tmp[1],tmp] = 1
     theta[tmp,tmp[1]] = 1
     rm(tmp)
     gc()
  }
  }
  if(graph == "random"){
    if(is.null(u)) u = 0.1
    if(is.null(v)) v = 0.3
  
    tmp = matrix(runif(d^2,0,0.5),d,d)
    tmp = tmp + t(tmp)
    theta[tmp < prob] = 1
    #theta[tmp >= tprob] = 0
    rm(tmp)
    gc()
  }

  if(graph == "scale-free"){
  if(is.null(u)) u = 0.1
  if(is.null(v)) v = 0.3
  out = .Call("_huge_SFGen", 2, d)
  theta = matrix(as.numeric(out$G),d,d)
  }
  diag(theta) = 0
  omega = theta*v

  # make omega positive definite and standardized
  diag(omega) = abs(min(eigen(omega)$values)) + 0.1 + u
  sigma = cov2cor(solve(omega))
  omega = solve(sigma)

  # generate multivariate normal data
  x = mvrnorm(n,rep(0,d),sigma)
  sigmahat = cor(x)

  # graph and covariance visulization
  if(vis == TRUE){
  fullfig = par(mfrow = c(2, 2), pty = "s", omi=c(0.3,0.3,0.3,0.3), mai = c(0.3,0.3,0.3,0.3))
  fullfig[1] = image(theta, col = gray.colors(256),  main = "Adjacency Matrix")

  fullfig[2] = image(sigma, col = gray.colors(256), main = "Covariance Matrix")
  g = graph.adjacency(theta, mode="undirected", diag=FALSE)
  layout.grid = layout.fruchterman.reingold(g)

  fullfig[3] = plot(g, layout=layout.grid, edge.color='gray50',vertex.color="red", vertex.size=3, vertex.label=NA,main = "Graph Pattern")

  fullfig[4] = image(sigmahat, col = gray.colors(256), main = "Empirical Matrix")
  rm(fullfig,g,layout.grid)
  gc()
  }
  if(verbose) cat("done.\n")
  rm(vis,verbose)
  gc()

  sim = list(data = x, sigma = sigma, sigmahat = sigmahat, omega = omega, theta = Matrix(theta,sparse = TRUE), sparsity= sum(theta)/(d*(d-1)), graph.type=graph)
  class(sim) = "sim"
  return(sim)
}

#' Print function for S3 class "sim"
#'
#' Print the information about the sample size, the dimension, the pattern and sparsity of the true graph strcture.
#'
#' @param x An object with S3 class \code{"sim"}.
#' @param \dots System reserved (No specific usage)
#' @seealso \code{\link{huge.generator}}
#' @export
print.sim = function(x, ...){
  cat("Simulated data generated by huge.generator()\n")
  cat("Sample size: n =", nrow(x$data), "\n")
  cat("Dimension: d =", ncol(x$data), "\n")
    cat("Graph type = ", x$graph.type, "\n")
    cat("Sparsity level:", sum(x$theta)/ncol(x$data)/(ncol(x$data)-1),"\n")
}

#' Plot function for S3 class "sim"
#'
#' Visualize the covariance matrix, the empirical covariance matrix, the adjacency matrix and the graph pattern of the true graph structure.
#'
#' @param x An object with S3 class \code{"sim"}
#' @param \dots System reserved (No specific usage)
#' @seealso \code{\link{huge.generator}} and \code{\link{huge}}
#' @export
plot.sim = function(x, ...){
  gcinfo(FALSE)
     par = par(mfrow = c(2, 2), pty = "s", omi=c(0.3,0.3,0.3,0.3), mai = c(0.3,0.3,0.3,0.3))
     image(as.matrix(x$theta), col = gray.colors(256),  main = "Adjacency Matrix")
  image(x$sigma, col = gray.colors(256), main = "Covariance Matrix")
  g = graph.adjacency(x$theta, mode="undirected", diag=FALSE)
  layout.grid = layout.fruchterman.reingold(g)

  plot(g, layout=layout.grid, edge.color='gray50',vertex.color="red", vertex.size=3, vertex.label=NA,main = "Graph Pattern")
  rm(g, layout.grid)
  gc()
  image(x$sigmahat, col = gray.colors(256), main = "Empirical Covariance Matrix")
}
