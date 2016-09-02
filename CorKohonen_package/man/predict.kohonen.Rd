\name{predict.kohonen}
\alias{predict.kohonen}
\title{Predict properties using a trained Kohonen map}
\description{
  Map objects to a trained Kohonen map, and return for each object the
  property associated with the corresponding winning unit. For
  \code{som} and \code{supersom} maps, the unit properties are
  calculated using explicit arguments \code{trainX} and \code{trainY};
  for \code{xyf} and \code{bdk} maps, the predicted properties are the
  Y-codebookvectors. Note that in the latter case only the X-space is
  used for prediction.
}
\usage{
\method{predict}{kohonen}(object, newdata, trainX, trainY, unit.predictions,
                 threshold = 0, whatmap = NULL, weights = 1, ...)
}
\arguments{
  \item{object}{Trained network.}
  \item{newdata}{Data matrix for which predictions are to be made. If
    not given, defaults to the training data (when available).}
  \item{trainX}{Training data for obtaining predictions for unsupervised
    maps; necessary for \code{som} maps trained with the
    \code{keep.data = FALSE} option.}
  \item{trainY}{Values for the dependent variable for the training data;
    necessary for \code{som} and \code{supersom} maps.}
  \item{unit.predictions}{Possible override of the predictions for each
    unit.}
  \item{threshold}{Used in class predictions; see
    \code{\link{classmat2classvec}}.}
  \item{whatmap}{For \code{supersom} maps: what layers to use in the
    mapping.}
  \item{weights}{For \code{supersom} maps: weights of layers uses in the
    mapping.}
  \item{\dots}{Currently not used.} 
}
\value{
  Returns a list with components
  \item{prediction}{predicted values for the properties of
    interest. When multiple values are predicted, this element is a
    list, otherwise a vector or a matrix.}
  \item{unit.classif}{unit numbers to which objects in the data matrix are
    mapped.}
  \item{unit.predictions}{mean values associated with map units. Again,
    when multiple properties are predicted, this is a list.}
}
\author{Ron Wehrens}
\seealso{\code{\link{som}},\code{\link{xyf}},\code{\link{bdk}},
  \code{\link{supersom}}, \code{\link{map}}}

\keyword{classif}
