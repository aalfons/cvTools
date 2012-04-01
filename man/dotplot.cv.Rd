\name{dotplot.cv}
\alias{dotplot.cv}
\alias{dotplot.cvSelect}
\title{Dot plots of cross-validation results}
\usage{
  \method{dotplot}{cv} (x, data, select = NULL,
    seFactor = NA, ...)

  \method{dotplot}{cvSelect} (x, data, subset = NULL,
    select = NULL, seFactor = x$seFactor, ...)
}
\arguments{
  \item{x}{an object inheriting from class
  \code{"cvSelect"} that contains cross-validation
  results.}

  \item{data}{currently ignored.}

  \item{subset}{a character, integer or logical vector
  indicating the subset of models for which to plot the
  cross-validation results.}

  \item{select}{a character, integer or logical vector
  indicating the columns of cross-validation results to be
  plotted.}

  \item{seFactor}{a numeric value giving the multiplication
  factor of the standard error for displaying error bars.
  Error bars can be suppressed by setting this to
  \code{NA}.}

  \item{\dots}{additional arguments to be passed to the
  \code{"formula"} method of
  \code{\link[lattice:xyplot]{dotplot}}.}
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
  Produce dot plots of (average) results from (repeated)
  \eqn{K}-fold cross-validation.
}
\details{
  For objects with multiple columns of repeated
  cross-validation results, conditional dot plots are
  produced.
}
\examples{
library("robustbase")
data("coleman")
set.seed(1234)  # set seed for reproducibility

## set up folds for cross-validation
folds <- cvFolds(nrow(coleman), K = 5, R = 10)


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
    folds = folds, trim = 0.1)

# combine and plot results
cvFits <- cvSelect(LS = cvFitLm, MM = cvFitLmrob, LTS = cvFitLts)
cvFits
dotplot(cvFits)


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

# combine and plot results
cvFitsLts <- cvSelect("0.5" = cvFitLts50, "0.75" = cvFitLts75)
cvFitsLts
dotplot(cvFitsLts)
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link{cvFit}}, \code{\link{cvSelect}},
  \code{\link{cvTuning}}, \code{\link[=plot.cv]{plot}},
  \code{\link[=xyplot.cvSelect]{xyplot}},
  \code{\link[=bwplot.cv]{bwplot}},
  \code{\link[=densityplot.cv]{densityplot}}
}
\keyword{hplot}

