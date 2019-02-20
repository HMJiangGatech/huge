#-------------------------------------------------------------------------#
# Package: High-dimensional Undirected Graph Estimation                   #
# huge.select(): Model selection using:                                   #
#                1.rotation information criterion (ric)                   #
#                2.stability approach to regularization selection (stars) #
#                3.extended Bayesian informaition criterion (ebic)        #
#-------------------------------------------------------------------------#

#' Model selection for high-dimensional undirected graph estimation
#'
#' Implements the regularization parameter selection for high dimensional undirected graph estimation. The optional approaches are rotation information criterion (ric), stability approach to regularization selection (stars) and extended Bayesian information criterion (ebic).
#'
#' Stability approach to regularization selection (stars) is a natural way to select optimal regularization parameter for all three estimation methods. It selects the optimal graph by variability of subsamplings and tends to overselect edges in Gaussian graphical models. Besides selecting the regularization parameters, stars can also provide an additional estimated graph by merging the corresponding subsampled graphs using the frequency counts. The subsampling procedure in stars may NOT be very efficient, we also provide the recent developed highly efficient, rotation information criterion approach (ric). Instead of tuning over a grid by cross-validation or subsampling, we directly estimate the optimal regularization parameter based on random Rotations. However, ric usually has very good empirical performances but suffers from underselections sometimes. Therefore, we suggest if user are sensitive of false negative rates, they should either consider increasing \code{r.num} or applying the stars to model selection. Extended Bayesian information criterion (ebic) is another competitive approach, but the \code{ebic.gamma} can only be tuned by experience.
#'
#' @param est An object with S3 class \code{"huge"}.
#' @param criterion Model selection criterion. \code{"ric"} and \code{"stars"} are available for all 3 graph estimation methods. \code{ebic} is only applicable when \code{est$method = "glasso"} in \code{huge()}. The default value is \code{"ric"}.
#' @param ebic.gamma The tuning parameter for ebic. The default value is 0.5. Only applicable when \code{est$method = "glasso"} and \code{criterion = "ebic"}.
#' @param stars.thresh The variability threshold in stars. The default value is \code{0.1}. An alternative value is \code{0.05}. Only applicable when \code{criterion = "stars"}.
#' @param stars.subsample.ratio The subsampling ratio. The default value is \code{10*sqrt(n)/n} when \code{n>144} and \code{0.8} when \code{n<=144}, where \code{n} is the sample size. Only applicable when \code{criterion = "stars"}.
#' @param rep.num The number of subsamplings when \code{criterion = "stars"} or rotations when \code{criterion = "ric"}. The default value is \code{20}. NOT applicable when \code{criterion = "ebic"}.
#' @param verbose If \code{verbose = FALSE}, tracing information printing is disabled. The default value is \code{TRUE}.
#' @return
#' An object with S3 class "select" is returned:
#'   \item{refit}{
#'     The optimal graph selected from the graph path
#'   }
#' \item{opt.icov}{
#'   The optimal precision matrix from the path only applicable when \code{method = "glasso"}
#' }
#' \item{opt.cov}{
#'   The optimal covariance matrix from the path only applicable when \code{method = "glasso"} and \code{est$cov} is available.
#' }
#' \item{merge}{
#'   The graph path estimated by merging the subsampling paths. Only applicable when the input \code{criterion = "stars"}.
#' }
#' \item{variability}{
#'   The variability along the subsampling paths. Only applicable when the input \code{criterion = "stars"}.
#' }
#' \item{ebic.scores}{
#'   Extended BIC scores for regularization parameter selection. Only applicable when \code{criterion = "ebic"}.
#' }
#' \item{opt.index}{
#'   The index of the selected regularization parameter. NOT applicable when the input \code{criterion = "ric"}
#' }
#' \item{opt.lambda}{
#'   The selected regularization/thresholding parameter.
#' }
#' \item{opt.sparsity}{
#'   The sparsity level of \code{"refit"}.
#' }
#'
#' and anything else included in the input \code{est}
#'
#' @note The model selection is NOT available when the data input is the sample covariance matrix.
#' @seealso \code{\link{huge}} and \code{\link{huge-package}}.
#' @examples
#' #generate data
#' L = huge.generator(d = 20, graph="hub")
#' out.mb = huge(L$data)
#' out.ct = huge(L$data, method = "ct")
#' out.glasso = huge(L$data, method = "glasso")
#'
#' #model selection using ric
#' out.select = huge.select(out.mb)
#' plot(out.select)
#'
#' #model selection using stars
#' #out.select = huge.select(out.ct, criterion = "stars", stars.thresh = 0.05,rep.num=10)
#' #plot(out.select)
#'
#' #model selection using ebic
#' out.select = huge.select(out.glasso,criterion = "ebic")
#' plot(out.select)
#' @export
huge.select = function(est, criterion = NULL, ebic.gamma = 0.5, stars.thresh = 0.1, stars.subsample.ratio = NULL, rep.num = 20, verbose = TRUE){

  gcinfo(FALSE)

  if(est$cov.input){
    cat("Model selection is not available when using the covariance matrix as input.")
    class(est) = "select"
      return(est)
  }
  if(!est$cov.input)
  {
    if(est$method == "mb"&&is.null(criterion))
      criterion = "ric"
    if(est$method == "ct"&&is.null(criterion))
      criterion = "stars"
    if(est$method == "glasso"&&is.null(criterion))
      criterion = "ebic"

    n = nrow(est$data)
    d = ncol(est$data)
    nlambda = length(est$lambda)

    if(criterion == "ric")
    {
      if(verbose)
      {
        cat("Conducting rotation information criterion (ric) selection....")
            flush.console()
          }

      if(n>rep.num){
        nr = rep.num
        r = sample(n,rep.num)
      }
      if(n<=rep.num){
        nr = n
        r = 1:n
      }

      est$opt.lambda = .Call("_huge_RIC", est$data,d,n,r,nr)*1.0/n
      gc()

      if(verbose){
        cat("done\n")
        flush.console()
      }

      if(verbose)
      {
        cat("Computing the optimal graph....")
        flush.console()
      }

      if(est$opt.lambda>max(cor(est$data)))
        est$refit = Matrix(0,d,d)
      else{

      if(est$method == "mb")
        est$refit = huge.mb(est$data, lambda = est$opt.lambda, sym = est$sym, idx.mat = est$idx.mat, verbose = FALSE)$path[[1]]
      if(est$method == "glasso")
      {
        if(!is.null(est$cov))
        {
          tmp = huge.glasso(est$data, lambda = est$opt.lambda, scr = est$scr, cov.output = TRUE, verbose = FALSE)
          est$opt.cov = tmp$cov[[1]]
        }
        if(is.null(est$cov))
          tmp = huge.glasso(est$data, lambda = est$opt.lambda, verbose = FALSE)

        est$refit = tmp$path[[1]]
        est$opt.icov = tmp$icov[[1]]
        rm(tmp)
        gc()
      }
      if(est$method == "ct")
        est$refit = huge.ct(est$data, lambda = est$opt.lambda, verbose = FALSE)$path[[1]]
      }
      est$opt.sparsity=sum(est$refit)/d/(d-1)

      if(verbose){
        cat("done\n")
        flush.console()
      }
    }

    if(criterion == "ebic"&&est$method == "glasso")
    {
      if(verbose)
      {
        cat("Conducting extended Bayesian information criterion (ebic) selection....")
            flush.console()
      }
      est$ebic.score = -n*est$loglik + log(n)*est$df + 4*ebic.gamma*log(d)*est$df
      est$opt.index = which.min(est$ebic.score)
      est$refit = est$path[[est$opt.index]]
      est$opt.icov = est$icov[[est$opt.index]]
      if(est$cov.output)
        est$opt.cov = est$cov[[est$opt.index]]
        est$opt.lambda = est$lambda[est$opt.index]
        est$opt.sparsity = est$sparsity[est$opt.index]
        if(verbose){
        cat("done\n")
        flush.console()
      }
    }

    if(criterion == "stars"){
      if(is.null(stars.subsample.ratio))
      {
        if(n>144) stars.subsample.ratio = 10*sqrt(n)/n
        if(n<=144) stars.subsample.ratio = 0.8
      }

      est$merge = list()
      for(i in 1:nlambda) est$merge[[i]] = Matrix(0,d,d)

        for(i in 1:rep.num)
        {
          if(verbose)
          {
          mes <- paste(c("Conducting Subsampling....in progress:", floor(100*i/rep.num), "%"), collapse="")
             cat(mes, "\r")
                flush.console()
        }
          ind.sample = sample(c(1:n), floor(n*stars.subsample.ratio), replace=FALSE)

          if(est$method == "mb")
            tmp = huge.mb(est$data[ind.sample,],lambda = est$lambda, scr = est$scr, idx.mat = est$idx.mat, sym = est$sym, verbose = FALSE)$path
             if(est$method == "ct")
            tmp = huge.ct(est$data[ind.sample,], lambda = est$lambda,verbose = FALSE)$path
          if(est$method == "glasso")
            tmp = huge.glasso(est$data[ind.sample,], lambda = est$lambda, scr = est$scr, verbose = FALSE)$path

          for(i in 1:nlambda)
            est$merge[[i]] = est$merge[[i]] + tmp[[i]]

          rm(ind.sample,tmp)
           gc()
      }

      if(verbose){
        mes = "Conducting Subsampling....done.                 "
          cat(mes, "\r")
          cat("\n")
          flush.console()
        }

        est$variability = rep(0,nlambda)
      for(i in 1:nlambda){
        est$merge[[i]] = est$merge[[i]]/rep.num
          est$variability[i] = 4*sum(est$merge[[i]]*(1-est$merge[[i]]))/(d*(d-1))
        }

        est$opt.index = max(which.max(est$variability >= stars.thresh)[1]-1,1)
         est$refit = est$path[[est$opt.index]]
        est$opt.lambda = est$lambda[est$opt.index]
        est$opt.sparsity = est$sparsity[est$opt.index]
        if(est$method == "glasso")
        {
          est$opt.icov = est$icov[[est$opt.index]]
        if(!is.null(est$cov))
          est$opt.cov = est$cov[[est$opt.index]]
        }
    }
      est$criterion = criterion
      class(est) = "select"
      return(est)
    }
}

#-----------------------------------------------------------------------#
# default printing function for class "select"                          #
#-----------------------------------------------------------------------#

#' Print function for S3 class "select"
#'
#' Print the information about the model usage, graph dimension, model selection criterion, sparsity level of the optimal graph.
#'
#' @param x An object with S3 class \code{"select"}.
#' @param \dots System reserved (No specific usage)
#' @seealso \code{\link{huge.select}}
#' @export
print.select = function(x, ...)
{
  if(x$cov.input){
    cat("Model selection is not available when using the covariance matrix as input.")
  }
  if(!x$cov.input)
  {
    if(x$method == "ct")
      cat("Model: graph gstimation via correlation thresholding(ct)\n")
    if(x$method == "glasso")
      cat("Model: graphical lasso (glasso)\n")
    if(x$method == "mb")
      cat("Model: Meinshausen & Buhlmann Graph Estimation (mb)\n")

    cat("selection criterion:",x$criterion,"\n")
    if((x$method != "ct")&&x$scr)
      cat("lossy screening: on\n")
    cat("Graph dimension:",ncol(x$data),"\n")
    cat("sparsity level", x$opt.sparsity,"\n")
  }
}

#' Plot function for S3 class "select"
#'
#' Plot the optimal graph by model selection.
#'
#' @param x An object with S3 class \code{"select"}
#' @param \dots System reserved (No specific usage)
#' @seealso \code{\link{huge.select}}
#' @export
plot.select = function(x, ...){
  if(x$cov.input){
    cat("Model selection is not available when using the covariance matrix as input.")
  }
  if(!x$cov.input)
  {
    par(mfrow=c(1,2), pty = "s", omi=c(0.3,0.3,0.3,0.3), mai = c(0.3,0.3,0.3,0.3))

      g = graph.adjacency(as.matrix(x$refit), mode="undirected", diag=FALSE)
    layout.grid = layout.fruchterman.reingold(g)

    plot(g, layout=layout.grid, edge.color='gray50',vertex.color="red", vertex.size=3, vertex.label=NA)
      plot(x$lambda, x$sparsity, log = "x", xlab = "Regularization Parameter", ylab = "Sparsity Level", type = "l",xlim = rev(range(x$lambda)), main = "Solution path sparsity levels")
      lines(x$opt.lambda,x$opt.sparsity,type = "p")
    }
}
