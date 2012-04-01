\name{subset.cv}
\alias{subset.cv}
\alias{subset.cvSelect}
\title{Subsetting cross-validation results}
\usage{
  \method{subset}{cv} (x, select = NULL, ...)

  \method{subset}{cvSelect} (x, subset = NULL,
    select = NULL, ...)
}
\arguments{
  \item{x}{an object inheriting from class \code{"cv"} or
  \code{"cvSelect"} that contains cross-validation
  results.}

  \item{subset}{a character, integer or logical vector
  indicating the subset of models for which to keep the
  cross-validation results.}

  \item{select}{a character, integer or logical vector
  indicating the columns of cross-validation results to be
  extracted.}

  \item{\dots}{currently ignored.}
}
\value{
  An object similar to \code{x} containing just the
  selected results.
}
\description{
  Extract subsets of results from (repeated) \eqn{K}-fold
  cross-validation.
}
\examples{
library("robustbase")
data("coleman")
set.seed(1234)  # set seed for reproducibility

## set up folds for cross-validation
folds <- cvFolds(nrow(coleman), K = 5, R = 10)


## compare raw and reweighted LTS estimators for 
## 50\% and 75\% subsets

# 50\% subsets
fitLts50 <- ltsReg(Y ~ ., data = coleman, alpha = 0.5)
cvFitLts50 <- cvLts(fitLts50, cost = rtmspe, folds = folds, 
    fit = "both", trim = 0.1)

# 75\% subsets
fitLts75 <- ltsReg(Y ~ ., data = coleman, alpha = 0.75)
cvFitLts75 <- cvLts(fitLts75, cost = rtmspe, folds = folds, 
    fit = "both", trim = 0.1)

# combine results into one object
cvFitsLts <- cvSelect("0.5" = cvFitLts50, "0.75" = cvFitLts75)
cvFitsLts

# extract reweighted LTS results with 50\% subsets
subset(cvFitLts50, select = "reweighted")
subset(cvFitsLts, subset = c(TRUE, FALSE), select = "reweighted")
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link{cvFit}}, \code{\link{cvSelect}},
  \code{\link{cvTuning}}, \code{\link{subset}}
}
\keyword{utilities}

