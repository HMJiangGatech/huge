#' @useDynLib huge, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @importFrom stats cov2cor cor runif qnorm sd pnorm
#' @importFrom graphics par image plot lines
#' @importFrom grDevices postscript gray.colors dev.off
#' @importFrom Matrix Matrix
#' @importFrom utils flush.console
#' @importFrom MASS mvrnorm
#' @importFrom igraph graph.adjacency layout.fruchterman.reingold
#' @importFrom methods new
NULL

#' High-Dimensional Undirected Graph Estimation
#'
#' A package for high-dimensional undirected graph estimation
#'
#' \tabular{ll}{
#'   Package: \tab huge\cr
#'   Type: \tab Package\cr
#'   Version: \tab 1.2.7\cr
#'   Date: \tab 2015-09-14\cr
#'   License: \tab GPL-2\cr
#'   LazyLoad: \tab yes\cr
#' }
#' The package "huge" provides 8 main functions:\cr
#' (1) the data generator creates random samples from multivariate normal distributions with different graph structures. Please refer to \code{\link{huge.generator}}.\cr
#' (2) the nonparanormal (npn) transformation helps relax the normality assumption. Please refer to \code{\link{huge.npn}}.\cr
#' (3) The correlation thresholding graph estimation. Please refer to \code{\link{huge}}.\cr
#' (4) The Meinshausen-Buhlmann graph estimation. Please refer to \code{\link{huge}}.\cr
#' (5) The graphical Lasso algorithm using lossless screening rule. Please refer and \code{\link{huge}}.\cr
#' \cr
#' **Both (4) and (5) can be further accelerated by the lossy screening rule preselecting the neighborhood of each node via thresholding sample correlation.
#' \cr
#' (6) The model selection using the stability approach to regularization selection. Please refer to \code{\link{huge.select}}.\cr
#' (7) The model selection using the rotation information criterion. Please refer to \code{\link{huge.select}}.\cr
#' (8) The model selection using the extended Bayesian information criterion. Please refer to \code{\link{huge.select}}.\cr
#' @docType package
#' @aliases huge-package
#' @author Tuo Zhao, Han Liu, Haoming Jiang, Kathryn Roeder, John Lafferty, and Larry Wasserman \cr
#' Maintainers: Haoming Jiang<hjiang98@gatech.edu>;
#' @references
#' 1.  T. Zhao and H. Liu. The huge Package for High-dimensional Undirected Graph Estimation in R. \emph{Journal of Machine Learning Research}, 2012\cr
#' 2.  H. Liu, F. Han, M. Yuan, J. Lafferty and L. Wasserman. High Dimensional Semiparametric Gaussian Copula Graphical Models. \emph{Annals of Statistics},2012 \cr
#' 3.  D. Witten and J. Friedman. New insights and faster computations for the graphical lasso. \emph{Journal of Computational and Graphical Statistics}, to appear, 2011.
#' 4.  Han Liu, Kathryn Roeder and Larry Wasserman. Stability Approach to Regularization Selection (StARS) for High Dimensional Graphical Models. \emph{Advances in Neural Information Processing Systems}, 2010.\cr
#' 5.  R. Foygel and M. Drton. Extended bayesian information criteria for gaussian graphical models. \emph{Advances in Neural Information Processing Systems}, 2010.\cr
#' 6.  H. Liu, J. Lafferty and L. Wasserman. The Nonparanormal: Semiparametric Estimation of High Dimensional Undirected Graphs. \emph{Journal of Machine Learning Research}, 2009 \cr
#' 7.  J. Fan and J. Lv. Sure independence screening for ultra-high dimensional feature space (with discussion). \emph{Journal of Royal Statistical Society B}, 2008.\cr
#' 8.  O. Banerjee, L. E. Ghaoui, A. d'Aspremont: Model Selection Through Sparse Maximum Likelihood Estimation for Multivariate Gaussian or Binary Data. \emph{Journal of Machine Learning Research}, 2008.\cr
#' 9.  J. Friedman, T. Hastie and R. Tibshirani. Regularization Paths for Generalized Linear Models via Coordinate Descent. \emph{Journal of Statistical Software}, 2008. \cr
#' 10. J. Friedman, T. Hastie and R. Tibshirani. Sparse inverse covariance estimation with the lasso, \emph{Biostatistics}, 2007.\cr
#' 11. N. Meinshausen and P. Buhlmann. High-dimensional Graphs and Variable Selection with the Lasso. \emph{The Annals of Statistics}, 2006.\cr
#' @seealso \code{\link{huge.generator}}, \code{\link{huge.npn}}, \code{\link{huge}}, \code{\link{huge.plot}} and \code{\link{huge.roc}}
"_PACKAGE"
#> [1] "_PACKAGE"

#' Stock price of S&P 500 companies from 2003 to 2008
#'
#' This data set consists of stock price and company information.
#'
#' This data set can be used to perform high-dimensional graph estimation to analyze the relationships between S&P 500 companies.
#'
#' @usage data(stockdata)
#' @format
#' The format is a list containing contains two matrices.
#' 1. data - 1258x452, represents the 452 stocks' close prices for 1258 trading days.
#'   2. info - 452x3:
#'   The 1st column: the query symbol for each company.
#'   The 2nd column: the category for each company.
#'   The 3rd column: the full name of each company.
#' @source It was publicly available at http://ichart.finance.yahoo.com, which is now out of date
#' @examples
#' data(stockdata)
#' image(stockdata$data)
#' stockdata$info
#' @keywords datasets
#' @docType data
"stockdata"
#> [1] "stockdata"
