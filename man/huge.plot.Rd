% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/huge.plot.R
\name{huge.plot}
\alias{huge.plot}
\title{Graph visualization}
\usage{
huge.plot(
  G,
  epsflag = FALSE,
  graph.name = "default",
  cur.num = 1,
  location = NULL
)
}
\arguments{
\item{G}{The adjacency matrix corresponding to the graph.}

\item{epsflag}{If \code{epsflag = TRUE}, save the plot as an eps file in the target directory. The default value is \code{FALSE}.}

\item{graph.name}{The name of the output eps files. The default value is "default".}

\item{cur.num}{The number of plots saved as eps files. Only applicale when \code{epsflag = TRUE}. The default value is 1.}

\item{location}{Target directory. The default value is the current working directory.}
}
\description{
Implements the graph visualization using adjacency matrix. It can automatic organize 2D embedding layout.
}
\details{
The user can change \code{cur.num} to plot several figures and select the best one. The implementation is based on the popular package "igraph".
}
\examples{
## visualize the hub graph
L = huge.generator(graph = "hub")
huge.plot(L$theta)

## visualize the band graph
L = huge.generator(graph = "band",g=5)
huge.plot(L$theta)

## visualize the cluster graph
L = huge.generator(graph = "cluster")
huge.plot(L$theta)

## plot 5 graphs and save the plots as eps files in the tempdir()
huge.plot(L$theta, epsflag = TRUE, cur.num = 5, location = tempdir())
}
\seealso{
\code{\link{huge}} and \code{\link{huge-package}}.
}
