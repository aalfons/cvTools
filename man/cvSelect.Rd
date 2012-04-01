\name{cvSelect}
\alias{cvSelect}
\alias{print.cvSelect}
\title{Model selection based on cross-validation}
\usage{
  cvSelect(..., .reshape = FALSE,
    .selectBest = c("min", "hastie"), .seFactor = 1)
}
\arguments{
  \item{\dots}{objects inheriting from class \code{"cv"} or
  \code{"cvSelect"} that contain cross-validation results.}

  \item{.reshape}{a logical indicating whether objects with
  more than one column of cross-validation results should
  be reshaped to have only one column (see
  \dQuote{Details}).}

  \item{.selectBest}{a character string specifying a
  criterion for selecting the best model.  Possible values
  are \code{"min"} (the default) or \code{"hastie"}.  The
  former selects the model with the smallest prediction
  error.  The latter is useful for nested models or for
  models with a tuning parameter controlling the complexity
  of the model (e.g., penalized regression).  It selects
  the most parsimonious model whose prediction error is no
  larger than \code{.seFactor} standard errors above the
  prediction error of the best overall model.  Note that
  the models are thereby assumed to be ordered from the
  most parsimonious one to the most complex one.  In
  particular a one-standard-error rule is frequently
  applied.}

  \item{.seFactor}{a numeric value giving a multiplication
  factor of the standard error for the selection of the
  best model.  This is ignored if \code{.selectBest} is
  \code{"min"}.}
}
\value{
  An object of class \code{"cvSelect"} with the following
  components:

  \item{n}{an integer giving the number of observations.}

  \item{K}{an integer vector giving the number of folds
  used in cross-validation for the respective model.}

  \item{R}{an integer vector giving the number of
  replications used in cross-validation for the respective
  model.}

  \item{best}{an integer vector giving the indices of the
  models with the best prediction performance.}

  \item{cv}{a data frame containing the estimated
  prediction errors for the models.  For models for which
  repeated cross-validation was performed, those are
  average values over all replications.}

  \item{se}{a data frame containing the estimated standard
  errors of the prediction loss for the models.}

  \item{selectBest}{a character string specifying the
  criterion used for selecting the best model.}

  \item{seFactor}{a numeric value giving the multiplication
  factor of the standard error used for the selection of
  the best model.}

  \item{reps}{a data frame containing the estimated
  prediction errors from all replications for those models
  for which repeated cross-validation was performed.  This
  is only returned if repeated cross-validation was
  performed for at least one of the models.}
}
\description{
  Combine cross-validation results for various models into
  one object and select the model with the best prediction
  performance.
}
\details{
  Keep in mind that objects inheriting from class
  \code{"cv"} or \code{"cvSelect"} may contain multiple
  columns of cross-validation results.  This is the case if
  the response is univariate but the
  \code{\link[stats]{predict}} method of the fitted model
  returns a matrix.

  The \code{.reshape} argument determines how to handle
  such objects.  If \code{.reshape} is \code{FALSE}, all
  objects are required to have the same number of columns
  and the best model for each column is selected.  A
  typical use case for this behavior would be if the
  investigated models contain cross-validation results for
  a raw and a reweighted fit.  It might then be of interest
  to researchers to compare the best model for the raw
  estimators with the best model for the reweighted
  estimators.

  If \code{.reshape} is \code{TRUE}, objects with more than
  one column of results are first transformed with
  \code{\link{cvReshape}} to have only one column.  Then
  the best overall model is selected.

  It should also be noted that the argument names of
  \code{.reshape}, \code{.selectBest} and \code{.seFacor}
  start with a dot to avoid conflicts with the argument
  names used for the objects containing cross-validation
  results.
}
\note{
  Even though the function allows to compare
  cross-validation results obtained with a different number
  of folds or a different number of replications, such
  comparisons should be made with care.  Hence warnings are
  issued in those cases.  For maximum comparability, the
  same data folds should be used in cross-validation for
  all models to be compared.
}
\examples{
library("robustbase")
data("coleman")
set.seed(1234)  # set seed for reproducibility

# set up folds for cross-validation
folds <- cvFolds(nrow(coleman), K = 5, R = 10)


## compare LS, MM and LTS regression

# perform cross-validation for an LS regression model
fitLm <- lm(Y ~ ., data = coleman)
cvFitLm <- cvLm(fitLm, cost = rtmspe, 
    folds = folds, trim = 0.1)

# perform cross-validation for an MM regression model
fitLmrob <- lmrob(Y ~ ., data = coleman)
cvFitLmrob <- cvLmrob(fitLmrob, cost = rtmspe, 
    folds = folds, trim = 0.1)

# perform cross-validation for an LTS regression model
fitLts <- ltsReg(Y ~ ., data = coleman)
cvFitLts <- cvLts(fitLts, cost = rtmspe, 
    folds = folds, trim = 0.1)

# compare cross-validation results
cvSelect(LS = cvFitLm, MM = cvFitLmrob, LTS = cvFitLts)


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
cvSelect("0.5" = cvFitLts50, "0.75" = cvFitLts75)
}
\author{
  Andreas Alfons
}
\references{
  Hastie, T., Tibshirani, R. and Friedman, J. (2009)
  \emph{The Elements of Statistical Learning: Data Mining,
  Inference, and Prediction}.  Springer, 2nd edition.
}
\seealso{
  \code{\link{cvFit}}, \code{\link{cvTuning}}
}
\keyword{utilities}

