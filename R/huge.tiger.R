#-----------------------------------------------------------------------#
# Package: High-dimensional Undirected Graph Estimation                 #
# huge.tiger(): Tuning-insensitive graph estimation                     #
#-----------------------------------------------------------------------#

#' Tuning-insensitive graph estimation
#' 
#' See more details in \code{\link{huge}}
#' @param x There are 2 options: (1) \code{x} is an \code{n} by \code{d} data matrix (2) a \code{d} by \code{d} sample covariance matrix. The program automatically identifies the input matrix by checking the symmetry. (\code{n} is the sample size and \code{d} is the dimension).
#' @param lambda A sequence of decreasing positive numbers to control the regularization when \code{method = "mb"}, \code{"glasso"} or \code{"tiger"}, or the thresholding in \code{method = "ct"}. Typical usage is to leave the input \code{lambda = NULL} and have the program compute its own \code{lambda} sequence based on \code{nlambda} and \code{lambda.min.ratio}. Users can also specify a sequence to override this. When \code{method = "mb"}, \code{"glasso"} or \code{"tiger"}, use with care - it is better to supply a decreasing sequence values than a single (small) value.
#' @param nlambda The number of regularization/thresholding parameters. The default value is \code{30} for \code{method = "ct"} and \code{10} for \code{method = "mb"}, \code{"glasso"} or \code{"tiger"}.
#' @param lambda.min.ratio If \code{method = "mb"}, \code{"glasso"} or \code{"tiger"}, it is the smallest value for \code{lambda}, as a fraction of the upperbound (\code{MAX}) of the regularization/thresholding parameter which makes all estimates equal to \code{0}. The program can automatically generate \code{lambda} as a sequence of length = \code{nlambda} starting from \code{MAX} to \code{lambda.min.ratio*MAX} in log scale. If \code{method = "ct"}, it is the largest sparsity level for estimated graphs. The program can automatically generate \code{lambda} as a sequence of length = \code{nlambda}, which makes the sparsity level of the graph path increases from \code{0} to \code{lambda.min.ratio} evenly.The default value is \code{0.1} when \code{method = "mb"}, \code{"glasso"} or \code{"tiger"}, and 0.05 \code{method = "ct"}.
#' @param sym Symmetrize the output graphs. If \code{sym = "and"}, the edge between node \code{i} and node \code{j} is selected ONLY when both node \code{i} and node \code{j} are selected as neighbors for each other. If \code{sym = "or"}, the edge is selected when either node \code{i} or node \code{j} is selected as the neighbor for each other. The default value is \code{"or"}. ONLY applicable when \code{method = "mb"} or {"tiger"}.
#' @param verbose If \code{verbose = FALSE}, tracing information printing is disabled. The default value is \code{TRUE}.
#' @seealso \code{\link{huge}}, and \code{\link{huge-package}}.
#' @export
huge.tiger = function(x, lambda = NULL, nlambda = NULL, lambda.min.ratio = NULL, sym = "or", verbose = TRUE)
{
	gcinfo(FALSE)
	n = nrow(x);
	d = ncol(x);
	fit = list()
	fit$cov.input = isSymmetric(x);
	if(fit$cov.input)
	{
		if(verbose) cat("The input is identified as the covariance matrix.\n")
		S = cov2cor(x);
	}
	if(!fit$cov.input)
	{
		x = scale(x)
		S = cor(x)
	}

	gc()

	if(!is.null(lambda)) nlambda = length(lambda)
	if(is.null(lambda))
	{
		if(is.null(nlambda))
			nlambda = 10
		if(is.null(lambda.min.ratio))
			lambda.min.ratio = 0.1
		lambda.max = max(max(S-diag(d)),-min(S-diag(d)))
		lambda.min = lambda.min.ratio*lambda.max
		lambda = exp(seq(log(lambda.max), log(lambda.min), length = nlambda))
		rm(lambda.max,lambda.min,lambda.min.ratio)
		gc()
	}
	maxdf = min(d,n);

	if(verbose)
	{
	  cat("Conducting graph estimation through a tunning-insensitive approach(tiger)....")
	  flush.console()
	}
	fit$idx_mat = NULL
	out = .Call("_huge_SPMBgraphsqrt", x, lambda, nlambda, d)
	for(i in 1:d)
	{
	  if(out$col_cnz[i+1]>out$col_cnz[i])
	  {
	    idx.tmp = (out$col_cnz[i]+1):out$col_cnz[i+1]
	    ord = order(out$row_idx[idx.tmp])
	    out$row_idx[idx.tmp] = out$row_idx[ord + out$col_cnz[i]]
	    out$x[idx.tmp] = out$x[ord + out$col_cnz[i]]
	  }
	}


  G = new("dgCMatrix", Dim = as.integer(c(d*nlambda,d)), x = as.vector(out$x[1:out$col_cnz[d+1]]),p = as.integer(out$col_cnz), i = as.integer(out$row_idx[1:out$col_cnz[d+1]]))

	fit$beta = list()
	fit$path = list()
	fit$df = matrix(0,d,nlambda)
	fit$rss = matrix(0,d,nlambda)
	fit$sparsity = rep(0,nlambda)
	for(i in 1:nlambda)
	{
		fit$beta[[i]] = G[((i-1)*d+1):(i*d),]
		fit$path[[i]] = abs(fit$beta[[i]])
		fit$df[,i] = apply(sign(fit$path[[i]]),2,sum)

		if(sym == "or")
			fit$path[[i]] = sign(fit$path[[i]] + t(as.matrix(fit$path[[i]])))
		if(sym == "and")
			fit$path[[i]] = sign(fit$path[[i]] * t(fit$path[[i]]))
		fit$sparsity[i] = sum(fit$path[[i]])/d/(d-1)
	}
	fit$icov = out$icov
 	rm(x, G, out)

 	fit$lambda = lambda

	if(verbose)
 	{
 		cat("done\n")
      flush.console()
  }

 	rm(verbose,nlambda)
 	gc()
 	return(fit)
}
