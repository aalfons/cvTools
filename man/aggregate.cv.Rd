\name{aggregate.cv}
\alias{aggregate.cv}
\alias{aggregate.cvSelect}
\alias{aggregate.cvTuning}
\title{Aggregate cross-validation results}
\usage{
  \method{aggregate}{cv} (x, FUN = mean, select = NULL,
    ...)

  \method{aggregate}{cvSelect} (x, FUN = mean,
    select = NULL, ...)

  \method{aggregate}{cvTuning} (x, ...)
}
\arguments{
  \item{x}{an object inheriting from class \code{"cv"} or
  \code{"cvSelect"} that contains cross-validation results
  (note that the latter includes objects of class
  \code{"cvTuning"}).}

  \item{FUN}{a function to compute the summary statistics.}

  \item{select}{a character, integer or logical vector
  indicating the columns of cross-validation results for
  which to compute the summary statistics.}

  \item{\dots}{for the \code{"cvTuning"} method, additional
  arguments to be passed to the \code{"cvSelect"} method.
  Otherwise additional arguments to be passed to
  \code{FUN}.}
}
\value{
  The \code{"cv"} method returns a vector or matrix of
  aggregated cross-validation results, depending on whether
  \code{FUN} returns a single value or a vector.

  For the other methods, a data frame containing the
  aggregated cross-validation results for each model is
  returned.  In the case of the \code{"cvTuning"} method,
  the data frame contains the combinations of tuning
  parameters rather than a column describing the models.
}
\description{
  Compute summary statistics of results from repeated
  \eqn{K}-fold cross-validation.
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

# summary of the results with the 50\% subsets
aggregate(cvFitLts50, summary)
# summary of the combined results
aggregate(cvFitsLts, summary)


## evaluate MM regression models tuned for 
## 80\%, 85\%, 90\% and 95\% efficiency
tuning <- list(tuning.psi=c(3.14, 3.44, 3.88, 4.68))

# set up function call
call <- call("lmrob", formula = Y ~ .)
# perform cross-validation
cvFitsLmrob <- cvTuning(call, data = coleman, 
    y = coleman$Y, tuning = tuning, cost = rtmspe, 
    folds = folds, costArgs = list(trim = 0.1))
cvFitsLmrob

# summary of results
aggregate(cvFitsLmrob, summary)
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link{cvFit}}, \code{\link{cvSelect}},
  \code{\link{cvTuning}}, \code{\link[stats]{aggregate}}
}
\keyword{utilities}

