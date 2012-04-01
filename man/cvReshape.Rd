\name{cvReshape}
\alias{cvReshape}
\alias{cvReshape.cv}
\alias{cvReshape.cvSelect}
\title{Reshape cross-validation results}
\usage{
  cvReshape(x, ...)

  \method{cvReshape}{cv} (x,
    selectBest = c("min", "hastie"), seFactor = 1, ...)

  \method{cvReshape}{cvSelect} (x,
    selectBest = c("min", "hastie"), seFactor = 1, ...)
}
\arguments{
  \item{x}{an object inheriting from class \code{"cv"} or
  \code{"cvSelect"} that contains cross-validation
  results.}

  \item{selectBest}{a character string specifying a
  criterion for selecting the best model.  Possible values
  are \code{"min"} (the default) or \code{"hastie"}.  The
  former selects the model with the smallest prediction
  error.  The latter is useful for nested models or for
  models with a tuning parameter controlling the complexity
  of the model (e.g., penalized regression).  It selects
  the most parsimonious model whose prediction error is no
  larger than \code{seFactor} standard errors above the
  prediction error of the best overall model.  Note that
  the models are thereby assumed to be ordered from the
  most parsimonious one to the most complex one.  In
  particular a one-standard-error rule is frequently
  applied.}

  \item{seFactor}{a numeric value giving a multiplication
  factor of the standard error for the selection of the
  best model.  This is ignored if \code{selectBest} is
  \code{"min"}.}

  \item{\dots}{additional arguments to be passed down.}
}
\value{
  An object of class \code{"cvSelect"} with the following
  components:

  \item{n}{an integer giving the number of observations.}

  \item{K}{an integer giving the number of folds used in
  cross-validation.}

  \item{R}{an integer giving the number of replications
  used in cross-validation.}

  \item{best}{an integer giving the index of the model with
  the best prediction performance.}

  \item{cv}{a data frame containing the estimated
  prediction errors for the models.  For repeated
  cross-validation, those are average values over all
  replications.}

  \item{se}{a data frame containing the estimated standard
  errors of the prediction loss for the models.}

  \item{selectBest}{a character string specifying the
  criterion used for selecting the best model.}

  \item{seFactor}{a numeric value giving the multiplication
  factor of the standard error used for the selection of
  the best model.}

  \item{reps}{a data frame containing the estimated
  prediction errors for the models from all replications.
  This is only returned if repeated cross-validation was
  performed.}
}
\description{
  Reshape cross-validation results into an object of class
  \code{"cvSelect"} with only one column of results.
}
\examples{
library("robustbase")
data("coleman")

# perform cross-validation for an LTS regression model
fitLts <- ltsReg(Y ~ ., data = coleman)
cvFitLts <- cvLts(fitLts, cost = rtmspe, K = 5, R = 10, 
    fit = "both", trim = 0.1, seed = 1234)
# compare original and reshaped object
cvFitLts
cvReshape(cvFitLts)
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
  \code{\link{cvFit}}, \code{\link{cvSelect}},
  \code{\link{cvTuning}}
}
\keyword{utilities}

