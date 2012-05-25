\name{cvFolds}
\alias{cvFolds}
\alias{print.cvFolds}
\title{Cross-validation folds}
\usage{
  cvFolds(n, K = 5, R = 1,
    type = c("random", "consecutive", "interleaved"),
    grouping = NULL)
}
\arguments{
  \item{n}{an integer giving the number of observations to
  be split into folds.  This is ignored if \code{grouping}
  is supplied in order to split groups of observations into
  folds.}

  \item{K}{an integer giving the number of folds into which
  the observations should be split (the default is five).
  Setting \code{K} equal to the number of observations or
  groups yields leave-one-out cross-validation.}

  \item{R}{an integer giving the number of replications for
  repeated \eqn{K}-fold cross-validation.  This is ignored
  for for leave-one-out cross-validation and other
  non-random splits of the data.}

  \item{type}{a character string specifying the type of
  folds to be generated.  Possible values are
  \code{"random"} (the default), \code{"consecutive"} or
  \code{"interleaved"}.}

  \item{grouping}{a factor specifying groups of
  observations.  If supplied, the data are split according
  to the groups rather than individual observations such
  that all observations within a group belong to the same
  fold.}
}
\value{
  An object of class \code{"cvFolds"} with the following
  components:

  \item{n}{an integer giving the number of observations or
  groups.}

  \item{K}{an integer giving the number of folds.}

  \item{R}{an integer giving the number of replications.}

  \item{subsets}{an integer matrix in which each column
  contains a permutation of the indices of the observations
  or groups.}

  \item{which}{an integer vector giving the fold for each
  permuted observation or group.}

  \item{grouping}{a list giving the indices of the
  observations belonging to each group.  This is only
  returned if a grouping factor has been supplied.}
}
\description{
  Split observations or groups of observations into \eqn{K}
  folds to be used for (repeated) \eqn{K}-fold
  cross-validation.  \eqn{K} should thereby be chosen such
  that all folds are of approximately equal size.
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

