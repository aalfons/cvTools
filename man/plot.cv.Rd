\name{plot.cv}
\alias{plot.cv}
\alias{plot.cvSelect}
\title{Plot cross-validation results}
\usage{
  \method{plot}{cv} (x,
    method = c("bwplot", "densityplot", "xyplot", "dotplot"),
    ...)

  \method{plot}{cvSelect} (x,
    method = c("bwplot", "densityplot", "xyplot", "dotplot"),
    ...)
}
\arguments{
  \item{x}{an object inheriting from class \code{"cv"} or
  \code{"cvSelect"} that contains cross-validation
  results.}

  \item{method}{a character string specifying the type of
  plot.  For the \code{"cv"} method, possible values are
  \code{"bwplot"} to create a box-and-whisker plot via
  \code{\link[=bwplot.cv]{bwplot}} (the default), or
  \code{"densityplot"} to create a kernel density plot via
  \code{\link[=densityplot.cv]{densityplot}}.  Note that
  those plots are only meaningful for results from repeated
  cross-validation.  For the \code{"cvSelect"} method,
  there are two additional possibilities: \code{"xyplot"}
  to plot the (average) results for each model via
  \code{\link[=xyplot.cvSelect]{xyplot}}, or
  \code{"dotplot"} to create a similar dot plot via
  \code{\link[=dotplot.cvSelect]{dotplot}}.  The default is
  to use \code{"bwplot"} for results from repeated
  cross-validation and \code{"xyplot"} otherwise.  In any
  case, partial string matching allows supply abbreviations
  of the accepted values.}

  \item{\dots}{additional arguments to be passed down.}
}
\value{
  An object of class \code{"trellis"} is returned
  invisibly.  The
  \code{\link[lattice:update.trellis]{update}} method can
  be used to update components of the object and the
  \code{\link[lattice:print.trellis]{print}} method
  (usually called by default) will plot it on an
  appropriate plotting device.
}
\description{
  Plot results from (repeated) \eqn{K}-fold
  cross-validation.
}
\details{
  For objects with multiple columns of cross-validation
  results, conditional plots are produced.
}
\examples{
library("robustbase")
data("coleman")
set.seed(1234)  # set seed for reproducibility

# set up folds for cross-validation
folds <- cvFolds(nrow(coleman), K = 5, R = 50)


## compare LS, MM and LTS regression

# perform cross-validation for an LS regression model
fitLm <- lm(Y ~ ., data = coleman)
cvFitLm <- cvLm(fitLm, cost = rtmspe, 
    folds = folds, trim = 0.1)

# perform cross-validation for an MM regression model
fitLmrob <- lmrob(Y ~ ., data = coleman, k.max = 500)
cvFitLmrob <- cvLmrob(fitLmrob, cost = rtmspe, 
    folds = folds, trim = 0.1)

# perform cross-validation for an LTS regression model
fitLts <- ltsReg(Y ~ ., data = coleman)
cvFitLts <- cvLts(fitLts, cost = rtmspe, 
    folds = folds,  trim = 0.1)

# combine results into one object
cvFits <- cvSelect(LS = cvFitLm, MM = cvFitLmrob, LTS = cvFitLts)
cvFits

# plot results for the MM regression model
plot(cvFitLmrob, method = "bw")
plot(cvFitLmrob, method = "density")

# plot combined results
plot(cvFits, method = "bw")
plot(cvFits, method = "density")
plot(cvFits, method = "xy")
plot(cvFits, method = "dot")


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

# plot combined results
plot(cvFitsLts, method = "bw")
plot(cvFitsLts, method = "density")
plot(cvFitsLts, method = "xy")
plot(cvFitsLts, method = "dot")
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link{cvFit}}, \code{\link{cvSelect}},
  \code{\link{cvTuning}}, \code{\link[=bwplot.cv]{bwplot}},
  \code{\link[=densityplot.cv]{densityplot}},
  \code{\link[=xyplot.cvSelect]{xyplot}},
  \code{\link[=dotplot.cvSelect]{dotplot}}
}
\keyword{hplot}

