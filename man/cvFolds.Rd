\name{cvFolds}
\alias{cvFolds}
\alias{print.cvFolds}
\title{Cross-validation folds}
\usage{
  cvFolds(n, K = 5, R = 1,
    type = c("random", "consecutive", "interleaved"))
}
\arguments{
  \item{n}{an integer giving the number of observations to
  be split into groups.}

  \item{K}{an integer giving the number of groups into
  which the observations should be split (the default is
  five).  Setting \code{K} equal to \code{n} yields
  leave-one-out cross-validation.}

  \item{R}{an integer giving the number of replications for
  repeated \eqn{K}-fold cross-validation.  This is ignored
  for for leave-one-out cross-validation and other
  non-random splits of the data.}

  \item{type}{a character string specifying the type of
  folds to be generated.  Possible values are
  \code{"random"} (the default), \code{"consecutive"} or
  \code{"interleaved"}.}
}
\value{
  An object of class \code{"cvFolds"} with the following
  components:

  \item{n}{an integer giving the number of observations.}

  \item{K}{an integer giving the number of folds.}

  \item{R}{an integer giving the number of replications.}

  \item{subsets}{an integer matrix in which each column
  contains a permutation of the indices.}

  \item{which}{an integer vector giving the fold for each
  permuted observation.}
}
\description{
  Split \eqn{n} observations into \eqn{K} groups to be used
  for (repeated) \eqn{K}-fold cross-validation.  \eqn{K}
  should thereby be chosen such that all groups are of
  approximately equal size.
}
\examples{
set.seed(1234)  # set seed for reproducibility
cvFolds(20, K = 5, type = "random")
cvFolds(20, K = 5, type = "consecutive")
cvFolds(20, K = 5, type = "interleaved")
cvFolds(20, K = 5, R = 10)
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link{cvFit}}, \code{\link{cvSelect}},
  \code{\link{cvTuning}}
}
\keyword{utilities}

