#-----------------------------------------------------------------------#
# Package: High-dimensional Undirected pathraph Estimation              #
# huge(): Draw ROC Curve for a solution path                            #
#         The ground truth is required                                  #
#-----------------------------------------------------------------------#

#' Draw ROC Curve for a graph path
#'
#' Draws ROC curve for a graph path according to the true graph structure.
#'
#' To avoid the horizontal oscillation, false positive rates is automatically sorted in the ascent order and true positive rates also follow the same order.
#'
#' @param path A graph path.
#' @param theta The true graph structure.
#' @param verbose If \code{verbose = FALSE}, tracing information printing is disabled. The default value is \code{TRUE}.
#' @note For a lasso regression,  the number of nonzero coefficients is at most \code{n-1}. If \code{d>>n}, even when regularization parameter is very small, the estimated graph may still be sparse. In this case, the AUC may not be a good choice to evaluate the performance.
#' @return
#' An object with S3 class "roc" is returned:
#'   \item{F1}{
#'     The F1 scores along the graph path.
#'   }
#' \item{tp}{
#'   The true positive rates along the graph path
#' }
#' \item{fp}{
#'   The false positive rates along the graph paths
#' }
#' \item{AUC}{
#'   Area under the ROC curve
#' }
#' @seealso \code{\link{huge}} and \code{\link{huge-package}}.
#' @examples
#' #generate data
#' L = huge.generator(d = 200, graph = "cluster", prob = 0.3)
#' out1 = huge(L$data)
#'
#' #draw ROC curve
#' Z1 = huge.roc(out1$path,L$theta)
#'
#' #Maximum F1 score
#' max(Z1$F1)
#' @export
huge.roc = function(path, theta, verbose = TRUE){
  gcinfo(verbose = FALSE)
  ROC = list()

  theta = as.matrix(theta)
  d = ncol(theta)
  pos.total = sum(theta!=0)
  neg.total = d*(d-1) - pos.total

  if(verbose) cat("Computing F1 scores, false positive rates and true positive rates....")
  ROC$tp = rep(0,length(path))
     ROC$fp = rep(0,length(path))
     ROC$F1 = rep(0,length(path))
     for (r in 1:length(path)){
       tmp = as.matrix(path[[r]])
       tp.all = (theta!=0)*(tmp!=0)
       diag(tp.all) = 0
    ROC$tp[r] <- sum(tp.all!=0)/pos.total
    fp.all = (theta==0)*(tmp!=0)
    diag(fp.all) = 0
    ROC$fp[r] <- sum(fp.all!=0)/neg.total

    fn = 1 - ROC$tp[r]
    precision = ROC$tp[r]/(ROC$tp[r]+ROC$fp[r])
    recall = ROC$tp[r]/(ROC$tp[r]+fn)
    ROC$F1[r] = 2*precision*recall/(precision+recall)
    if(is.na(ROC$F1[r]))  ROC$F1[r] = 0
  }
  if(verbose) cat("done.\n")

  rm(precision,recall,tp.all,fp.all,path,theta,fn)
     gc()

  ord.fp = order(ROC$fp)

  tmp1 = ROC$fp[ord.fp]
  tmp2 = ROC$tp[ord.fp]
  par(mfrow = c(1,1))
  plot(tmp1,tmp2,type="b",main = "ROC Curve", xlab = "False Postive Rate", ylab = "True Postive Rate",ylim = c(0,1))
  ROC$AUC = sum(diff(tmp1)*(tmp2[-1]+tmp2[-length(tmp2)]))/2

  rm(ord.fp, tmp1, tmp2)
  gc()
  class(ROC) = "roc"
  return(ROC)
}

#' Print function for S3 class "roc"
#'
#' Print the information about true positive rates, false positive rates, the area under curve and maximum F1 score.
#'
#' @param x An object with S3 class \code{"roc"}.
#' @param \dots System reserved (No specific usage)
#' @seealso \code{\link{huge.roc}}
#' @export
print.roc = function(x, ...){
  cat("True Postive Rate: from",min(x$tp),"to",max(x$tp),"\n")
  cat("False Positive Rate: from",min(x$fp),"to",max(x$fp),"\n")
  cat("Area under Curve:",x$AUC,"\n")
  cat("Maximum F1 Score:",max(x$F1),"\n")
}

#' Plot function for S3 class "roc"
#'
#' Plot the ROC curve for an object with S3 class \code{"roc"}.
#'
#' @param x An object with S3 class \code{"roc"}
#' @param \dots System reserved (No specific usage)
#' @seealso \code{\link{huge.roc}}
#' @export
plot.roc = function(x, ...){
  ord.fp = order(x$fp)
  par(mfrow = c(1,1))
  plot(x$fp[ord.fp],x$tp[ord.fp],type="b",main = "ROC Curve", xlab = "False Postive Rate", ylab = "True Postive Rate",ylim = c(0,1))
}
