#-----------------------------------------------------------------------#
# Package: High-dimensional Undirecte d Graph Estimation                #
# huge(): The user interface for huge.mb(), huge.glasso() and huge.ct() #
# Authors: Tuo Zhao and Han Liu                                         #
# Emails: <tzhao5@jhu.edu> and <hanliu@cs.jhu.edu>                      #
# Date: Jul 15th 2011                                                   #
# Version: 1.1.0                                                        #
#-----------------------------------------------------------------------#

## Main function
#' High-dimensional undirected graph estimation
#' 
#' The main function for high-dimensional undirected graph estimation. Three graph estimation methods, including (1) Meinshausen-Buhlmann graph estimation (\code{mb}) (2) graphical lasso (\code{glasso}) and (3) correlation thresholding graph estimation (\code{ct}), are available for data analysis.
#' 
#' The graph structure is estimated by Meinshausen-Buhlmann graph estimation or the graphical lasso, and both methods can be further accelerated via the lossy screening rule by preselecting the neighborhood of each variable by correlation thresholding. We target on high-dimensional data analysis usually d >> n, and the computation is memory-optimized using the sparse matrix output. We also provide a highly computationally efficient approaches correlation thresholding graph estimation.
#' 
#' @param x There are 2 options: (1) \code{x} is an \code{n} by \code{d} data matrix (2) a \code{d} by \code{d} sample covariance matrix. The program automatically identifies the input matrix by checking the symmetry. (\code{n} is the sample size and \code{d} is the dimension).
#' @param lambda A sequence of decresing positive numbers to control the regularization when \code{method = "mb"} or \code{"glasso"}, or the thresholding in \code{method = "ct"}. Typical usage is to leave the input \code{lambda = NULL} and have the program compute its own \code{lambda} sequence based on \code{nlambda} and \code{lambda.min.ratio}. Users can also specify a sequence to override this. When \code{method = "mb"} or \code{"glasso"}, use with care - it is better to supply a decreasing sequence values than a single (small) value.
#' @param nlambda The number of regularization/thresholding paramters. The default value is \code{30} for \code{method = "ct"} and \code{10} for \code{method = "mb"} or \code{"glasso"}.
#' @param lambda.min.ratio If \code{method = "mb"} or \code{"glasso"}, it is the smallest value for \code{lambda}, as a fraction of the uppperbound (\code{MAX}) of the regularization/thresholding parameter which makes all estimates equal to \code{0}. The program can automatically generate \code{lambda} as a sequence of length = \code{nlambda} starting from \code{MAX} to \code{lambda.min.ratio*MAX} in log scale. If \code{method = "ct"}, it is the largest sparsity level for estimated graphs. The program can automatically generate \code{lambda} as a sequence of length = \code{nlambda}, which makes the sparsity level of the graph path increases from \code{0} to \code{lambda.min.ratio} evenly.The default value is \code{0.1} when \code{method = "mb"} or \code{"glasso"}, and 0.05 \code{method = "ct"}.
#' @param method Graph estimation methods with 3 options: \code{"mb"}, \code{"ct"} and \code{"glasso"}. The defaulty value is \code{"mb"}. 
#' @param scr If \code{scr = TRUE}, the lossy screening rule is applied to preselect the neighborhood before the graph estimation. The default value is  \code{FALSE}. NOT applicable when \code{method = "ct"}.
#' @param scr.num The neighborhood size after the lossy screening rule (the number of remaining neighbors per node). ONLY applicable when \code{scr = TRUE}. The default value is \code{n-1}. An alternative value is \code{n/log(n)}. ONLY applicable when \code{scr = TRUE} and \code{method = "mb"}.
#' @param cov.output If \code{cov.output = TRUE}, the output will inlcude a path of estimated covariance matrices. ONLY applicable when \code{method = "glasso"}. Since the estimated covariance matrices are generally not sparse, please use it with care, or it may take much memory under high-dimensional setting. The default value is \code{FALSE}.
#' @param sym Symmetrize the output graphs. If \code{sym = "and"}, the edge between node \code{i} and node \code{j} is selected ONLY when both node \code{i} and node \code{j} are selected as neighbors for each other. If \code{sym = "or"}, the edge is selected when either node \code{i} or node \code{j} is selected as the neighbor for each other. The default value is \code{"or"}. ONLY applicable when \code{method = "mb"}.
#' @param verbose If \code{verbose = FALSE}, tracing information printing is disabled. The default value is \code{TRUE}.
#' @return 
#' An object with S3 class \code{"huge"} is returned:  
#'   \item{data}{
#'     The \code{n} by \code{d} data matrix or \code{d} by \code{d} sample covariance matrix from the input
#'   }
#' \item{cov.input}{
#'   An indicator of the sample covariance. 
#' }
#' \item{ind.mat}{
#'   The \code{scr.num} by \code{k} matrix with each column correspondsing to a variable in \code{ind.group} and contains the indices of the remaining neighbors after the GSS. ONLY applicable when \code{scr = TRUE} and \code{approx = FALSE}
#' }
#' \item{lambda}{
#'   The sequence of regularization parameters used in mb or thresholding parameters in ct.
#' }
#' \item{sym}{
#'   The \code{sym} from the input. ONLY applicable when \code{method = "mb"}.
#' }
#' \item{scr}{
#'   The \code{scr} from the input. ONLY applicable when \code{method = "mb"} or {"glasso"}.
#' }
#' \item{path}{
#'   A list of \code{k} by \code{k} adjacency matrices of estimated graphs as a graph path corresponding to \code{lambda}.
#' }
#' \item{sparsity}{
#'   The sparsity levels of the graph path.
#' }
#' \item{icov}{
#'   A list of \code{d} by \code{d} precision matrices as an alternative graph path (numerical path) corresponding to \code{lambda}. ONLY applicable when {method = "glasso"}
#' }
#' \item{cov}{
#'   A list of \code{d} by \code{d} estimated covariance matrices corresponding to \code{lambda}. ONLY applicable when \code{cov.output = TRUE} and {method = "glasso"}
#' }
#' \item{method}{
#'   The method used in the graph estimation stage.
#' }
#' \item{df}{
#'   If \code{method = "mb"}, it is a \code{k} by \code{nlambda} matrix. Each row contains the number of nonzero coefficients along the lasso solution path. If \code{method = "glasso"}, it is a \code{nlambda} dimensional vector containing the number of nonzero coefficients along the graph path \code{icov}.
#' }
#' \item{loglik}{
#'   A \code{nlambda} dimensional vector containing the likelihood scores along the graph path (\code{icov}). ONLY applicable when \code{method = "glasso"}. For an estimated inverse convariance Z, the program only calculates log(det(Z)) - trace(SZ) where S is the empirical covariance matrix. For the likelihood for n observations, please multiply by n/2.
#' }
#' @param data The \code{n} by \code{d} data matrix or \code{d} by \code{d} sample covariance matrix from the input
#' @examples
#' #generate data
#' L = huge.generator(n = 50, d = 12, graph = "hub", g = 4)
#' 
#' #graph path estimation using mb
#' out1 = huge(L$data)
#' out1
#' plot(out1)				 #Not aligned	
#' plot(out1, align = TRUE) #Aligned
#' huge.plot(out1$path[[3]])
#' 
#' #graph path estimation using the sample covariance matrix as the input.
#' #out1 = huge(cor(L$data))
#' #out1
#' #plot(out1)				 #Not aligned	
#' #plot(out1, align = TRUE) #Aligned
#' #huge.plot(out1$path[[3]])
#' 
#' #graph path estimation using ct
#' #out2 = huge(L$data,method = "ct")
#' #out2
#' #plot(out2)
#' 
#' #graph path estimation using glasso
#' #out3 = huge(L$data, method = "glasso")
#' #out3
#' #plot(out3)
#' @export
huge = function(x, lambda = NULL, nlambda = NULL, lambda.min.ratio = NULL, method = "mb", scr = NULL, scr.num = NULL, cov.output = FALSE, sym = "or", verbose = TRUE)
{	
	gcinfo(FALSE)
	est = list()
	est$method = method	
		
	if(method == "ct")
	{
		fit = huge.ct(x, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, lambda = lambda, verbose = verbose)
		est$path = fit$path
		est$lambda = fit$lambda
		est$sparsity = fit$sparsity
		est$cov.input = fit$cov.input
		rm(fit)
		gc()
	}
	
	if(method == "mb")
	{	
		fit = huge.mb(x, lambda = lambda, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, scr = scr, scr.num = scr.num, sym = sym, verbose = verbose)
		est$path = fit$path
		est$lambda = fit$lambda
		est$sparsity = fit$sparsity
		est$df = fit$df
		est$idx_mat = fit$idx_mat
		est$sym = sym
		est$scr = fit$scr
		est$cov.input = fit$cov.input
		rm(fit,sym)
		gc()
	}
	
	
	if(method == "glasso")
	{
		fit = huge.glasso(x, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, lambda = lambda, scr = scr, cov.output = cov.output, verbose = verbose)
		est$path = fit$path
		est$lambda = fit$lambda
		est$icov = fit$icov
		est$df = fit$df
		est$sparsity = fit$sparsity
		est$loglik = fit$loglik
		if(cov.output)
			est$cov = fit$cov
		est$cov.input = fit$cov.input
		est$cov.output = fit$cov.output	
		est$scr = fit$scr
		rm(fit)
		gc()
	}			
			
	est$data = x
	
	rm(x,scr,lambda,lambda.min.ratio,nlambda,cov.output,verbose)
	gc()
	class(est) = "huge"
	return(est)
}

print.huge = function(x, ...)
{	
	if(x$method == "ct")
		cat("Model: graph estimation via correlation thresholding (ct)\n")
	if(x$method == "glasso")
		cat("Model: graphical lasso (glasso)\n")
	if(x$method == "mb")
		cat("Model: Meinshausen & Buhlmann graph estimation (mb)\n")
	
	if((x$method != "ct")&&(x$scr)) cat("lossy screening: on\n")

	if(x$cov.input) cat("Input: The Covariance Matrix\n")
	if(!x$cov.input) cat("Input: The Data Matrix\n")
	
	cat("Path length:",length(x$lambda),"\n")
	cat("Graph dimension:",ncol(x$data),"\n")
	cat("Sparsity level:",min(x$sparsity),"----->",max(x$sparsity),"\n")
}

plot.huge = function(x, align = FALSE, ...){
  gcinfo(FALSE)
	
	if(length(x$lambda) == 1)	par(mfrow = c(1, 2), pty = "s", omi=c(0.3,0.3,0.3,0.3), mai = c(0.3,0.3,0.3,0.3))
	if(length(x$lambda) == 2)	par(mfrow = c(1, 3), pty = "s", omi=c(0.3,0.3,0.3,0.3), mai = c(0.3,0.3,0.3,0.3))
	if(length(x$lambda) >= 3)	par(mfrow = c(1, 4), pty = "s", omi=c(0.3,0.3,0.3,0.3), mai = c(0.3,0.3,0.3,0.3))
	
	if(length(x$lambda) <= 3)	z.final = 1:length(x$lambda)
	
	if(length(x$lambda) >=4){
		z.max = max(x$sparsity)
		z.min = min(x$sparsity)
		z = z.max - z.min
		z.unique = unique(c(which(x$sparsity>=(z.min + 0.03*z))[1],which(x$sparsity>=(z.min + 0.07*z))[1],which(x$sparsity>=(z.min + 0.15*z))[1]))

		
		if(length(z.unique) == 1){
			if(z.unique<(length(x$lambda)-1))	z.final = c(z.unique,z.unique+1,z.unique+2)
			if(z.unique==(length(x$lambda)-1)) z.final = c(z.unique-1,z.unique,z.unique+1)
			if(z.unique==length(x$lambda)) 	z.final = c(z.unique-2,z.unique-1,z.unique)
		}
		
		if(length(z.unique) == 2){
			if(diff(z.unique)==1){
				if(z.unique[2]<length(x$lambda)) z.final = c(z.unique,z.unique[2]+1) 
				if(z.unique[2]==length(x$lambda)) z.final = c(z.unique[1]-1,z.unique)
			}
			if(diff(z.unique)>1) z.final = c(z.unique[1],z.unique[1]+1,z.unique[2])
		}
		
		if(length(z.unique) == 3) z.final = z.unique
		
		rm(z.max,z.min,z,z.unique)
		gc()
		
	}
	plot(x$lambda, x$sparsity, log = "x", xlab = "Regularization Parameter", ylab = "Sparsity Level", type = "l",xlim = rev(range(x$lambda)), main = "Sparsity vs. Regularization")
	
	lines(x$lambda[z.final],x$sparsity[z.final],type = "p")
	
	if(align){
		layout.grid = layout.fruchterman.reingold(graph.adjacency(as.matrix(x$path[[z.final[length(z.final)]]]), mode="undirected", diag=FALSE))
		for(i in z.final){
			g = graph.adjacency(as.matrix(x$path[[i]]), mode="undirected", diag=FALSE)
			plot(g, layout=layout.grid, edge.color='gray50',vertex.color="red", vertex.size=3, vertex.label=NA, main = paste("lambda = ",as.character(round(x$lambda[i],3)),sep = ""))
			rm(g)
			gc()
		}
	rm(layout.grid)
	}
	if(!align){
		for(i in z.final){
			g = graph.adjacency(as.matrix(x$path[[i]]), mode="undirected", diag=FALSE)
			layout.grid = layout.fruchterman.reingold(g)
			plot(g, layout=layout.grid, edge.color='gray50',vertex.color="red", vertex.size=3, vertex.label=NA, main = paste("lambda = ",as.character(round(x$lambda[i],3)),sep = ""))
			rm(g,layout.grid)
			gc()
		}
	}
	if(align) cat("Three plotted graphs are aligned according to the third graph\n")
}