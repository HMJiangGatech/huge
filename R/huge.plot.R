#-----------------------------------------------------------------------#
# Package: High-dimensional Undirected Graph Estimation                 #
# huge.plot(): graph visualization                                      #
#-----------------------------------------------------------------------#

#' Graph visualization
#'
#' Implements the graph visualization using adjacency matrix. It can automatic organize 2D embedding layout.
#'
#' The user can change \code{cur.num} to plot several figures and select the best one. The implementation is based on the popular package "igraph".
#'
#' @param G The adjacency matrix corresponding to the graph.
#' @param epsflag If \code{epsflag = TRUE}, save the plot as an eps file in the target directory. The default value is \code{FALSE}.
#' @param graph.name The name of the output eps files. The default value is "default".
#' @param cur.num The number of plots saved as eps files. Only applicale when \code{epsflag = TRUE}. The default value is 1.
#' @param location Target directory. The default value is the current working directory.
#' @seealso \code{\link{huge}} and \code{\link{huge-package}}.
#' @examples
#' The following code are commented out for passing the CRAN check
#' ## visualize the hub graph
#' L = huge.generator(graph = "hub")
#' huge.plot(L$theta)
#'
#' ## visualize the band graph
#' L = huge.generator(graph = "band",g=5)
#' huge.plot(L$theta)
#'
#' ## visualize the cluster graph
#' L = huge.generator(graph = "cluster")
#' huge.plot(L$theta)
#'
#' #plot 5 graphs and save the plots as eps files in the tempdir()
#' huge.plot(L$theta, epsflag = TRUE, cur.num = 5, location = tempdir())
#' @export
huge.plot = function(G, epsflag = FALSE, graph.name = "default", cur.num = 1, location=NULL){
  gcinfo(FALSE)
  if(missing(location))  location = tempdir()
  oldlocation = getwd()
  setwd(location)
  g = graph.adjacency(as.matrix(G!=0), mode="undirected", diag=FALSE)
  layout.grid = layout.fruchterman.reingold(g)

  if(epsflag == TRUE)  postscript(paste(paste(graph.name, cur.num, sep=""), "eps", sep="."), width = 8.0, height = 8.0)
    par(mfrow = c(1,1))
  plot(g, layout=layout.grid, edge.color='gray50',vertex.color="red", vertex.size=2, vertex.label=NA)
  rm(g,location)
  gc()
  if(epsflag == TRUE) dev.off()
  setwd(oldlocation)
}
