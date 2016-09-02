\name{unit.distances}
\alias{unit.distances}
\title{Calculate distances between units in a SOM}
\description{
Calculate distances between units in a SOM. 
}
\usage{
unit.distances(grid, toroidal)
}
\arguments{
  \item{grid}{an object of class \code{somgrid}.}
  \item{toroidal}{if true, edges of the map are joined so that the
    topology is that of a torus.}
}
\value{
  Returns a (symmetrical) matrix containing distances. When
  \code{grid$n.hood} equals "circular", Euclidean distances are used; for
  \code{grid$n.hood} is "square" maximum distances. If \code{toroidal}
  equals \code{TRUE}, maps are joined at the edges and distances are
  calculated for the shortest path.
}
\author{Ron Wehrens}

\keyword{classif}
