\name{cvFit}
\alias{cvFit}
\alias{cvFit.call}
\alias{cvFit.default}
\alias{cvFit.function}
\alias{print.cv}
\title{Cross-validation for model evaluation}
\usage{
  cvFit(object, ...)

  \method{cvFit}{default} (object, data = NULL, x = NULL,
    y, cost = rmspe, K = 5, R = 1,
    foldType = c("random", "consecutive", "interleaved"),
    folds = NULL, names = NULL, predictArgs = list(),
    costArgs = list(), envir = parent.frame(), seed = NULL,
    ...)

  \method{cvFit}{function} (object, formula, data = NULL,
    x = NULL, y, args = list(), cost = rmspe, K = 5, R = 1,
    foldType = c("random", "consecutive", "interleaved"),
    folds = NULL, names = NULL, predictArgs = list(),
    costArgs = list(), envir = parent.frame(), seed = NULL,
    ...)

  \method{cvFit}{call} (object, data = NULL, x = NULL, y,
    cost = rmspe, K = 5, R = 1,
    foldType = c("random", "consecutive", "interleaved"),
    folds = NULL, names = NULL, predictArgs = list(),
    costArgs = list(), envir = parent.frame(), seed = NULL,
    ...)
}
\arguments{
  \item{object}{the fitted model for which to estimate the
  prediction error, a function for fitting a model, or an
  unevaluated function call for fitting a model (see
  \code{\link{call}} for the latter).  In the case of a
  fitted model, the object is required to contain a
  component \code{call} that stores the function call used
  to fit the model, which is typically the case for objects
  returned by model fitting functions.}

  \item{formula}{a \code{\link[stats]{formula}} describing
  the model.}

  \item{data}{a data frame containing the variables
  required for fitting the models.  This is typically used
  if the model in the function call is described by a
  \code{\link[stats]{formula}}.}

  \item{x}{a numeric matrix containing the predictor
  variables.  This is typically used if the function call
  for fitting the models requires the predictor matrix and
  the response to be supplied as separate arguments.}

  \item{y}{a numeric vector or matrix containing the
  response.}

  \item{args}{a list of additional arguments to be passed
  to the model fitting function.}

  \item{cost}{a cost function measuring prediction loss.
  It should expect the observed values of the response to
  be passed as the first argument and the predicted values
  as the second argument, and must return either a
  non-negative scalar value, or a list with the first
  component containing the prediction error and the second
  component containing the standard error.  The default is
  to use the root mean squared prediction error (see
  \code{\link{cost}}).}

  \item{K}{an integer giving the number of groups into
  which the data should be split (the default is five).
  Keep in mind that this should be chosen such that all
  groups are of approximately equal size.  Setting \code{K}
  equal to \code{n} yields leave-one-out cross-validation.}

  \item{R}{an integer giving the number of replications for
  repeated \eqn{K}-fold cross-validation.  This is ignored
  for for leave-one-out cross-validation and other
  non-random splits of the data.}

  \item{foldType}{a character string specifying the type of
  folds to be generated.  Possible values are
  \code{"random"} (the default), \code{"consecutive"} or
  \code{"interleaved"}.}

  \item{folds}{an object of class \code{"cvFolds"} giving
  the folds of the data for cross-validation (as returned
  by \code{\link{cvFolds}}).  If supplied, this is
  preferred over \code{K} and \code{R}.}

  \item{names}{an optional character vector giving names
  for the arguments containing the data to be used in the
  function call (see \dQuote{Details}).}

  \item{predictArgs}{a list of additional arguments to be
  passed to the \code{\link[stats]{predict}} method of the
  fitted models.}

  \item{costArgs}{a list of additional arguments to be
  passed to the prediction loss function \code{cost}.}

  \item{envir}{the \code{\link{environment}} in which to
  evaluate the function call for fitting the models (see
  \code{\link{eval}}).}

  \item{seed}{optional initial seed for the random number
  generator (see \code{\link{.Random.seed}}).}

  \item{\dots}{additional arguments to be passed down.}
}
\value{
  An object of class \code{"cv"} with the following
  components:

  \item{n}{an integer giving the number of observations.}

  \item{K}{an integer giving the number of folds.}

  \item{R}{an integer giving the number of replications.}

  \item{cv}{a numeric vector containing the respective
  estimated prediction errors.  For repeated
  cross-validation, those are average values over all
  replications.}

  \item{se}{a numeric vector containing the respective
  estimated standard errors of the prediction loss.}

  \item{reps}{a numeric matrix in which each column
  contains the respective estimated prediction errors from
  all replications.  This is only returned for repeated
  cross-validation.}

  \item{seed}{the seed of the random number generator
  before cross-validation was performed.}

  \item{call}{the matched function call.}
}
\description{
  Estimate the prediction error of a model via (repeated)
  \eqn{K}-fold cross-validation.  It is thereby possible to
  supply an object returned by a model fitting function, a
  model fitting function itself, or an unevaluated function
  call to a model fitting function.
}
\details{
  (Repeated) \eqn{K}-fold cross-validation is performed in
  the following way.  The data are first split into \eqn{K}
  previously obtained blocks of approximately equal size.
  Each of the \eqn{K} data blocks is left out once to fit
  the model, and predictions are computed for the
  observations in the left-out block with the
  \code{\link[stats]{predict}} method of the fitted model.
  Thus a prediction is obtained for each observation.

  The response variable and the obtained predictions for
  all observations are then passed to the prediction loss
  function \code{cost} to estimate the prediction error.
  For repeated cross-validation, this process is replicated
  and the estimated prediction errors from all replications
  as well as their average are included in the returned
  object.

  Furthermore, if the response is a vector but the
  \code{\link[stats]{predict}} method of the fitted models
  returns a matrix, the prediction error is computed for
  each column.  A typical use case for this behavior would
  be if the \code{\link[stats]{predict}} method returns
  predictions from an initial model fit and stepwise
  improvements thereof.

  If \code{formula} or \code{data} are supplied, all
  variables required for fitting the models are added as
  one argument to the function call, which is the typical
  behavior of model fitting functions with a
  \code{\link[stats]{formula}} interface.  In this case,
  the accepted values for \code{names} depend on the
  method.  For the \code{function} method, a character
  vector of length two should supplied, with the first
  element specifying the argument name for the formula and
  the second element specifying the argument name for the
  data (the default is to use \code{c("formula", "data")}).
  Note that names for both arguments should be supplied
  even if only one is actually used.  For the other
  methods, which do not have a \code{formula} argument, a
  character string specifying the argument name for the
  data should be supplied (the default is to use
  \code{"data"}).

  If \code{x} is supplied, on the other hand, the predictor
  matrix and the response are added as separate arguments
  to the function call.  In this case, \code{names} should
  be a character vector of length two, with the first
  element specifying the argument name for the predictor
  matrix and the second element specifying the argument
  name for the response (the default is to use \code{c("x",
  "y")}).  It should be noted that the \code{formula} or
  \code{data} arguments take precedence over \code{x}.
}
\examples{
library("robustbase")
data("coleman")

## via model fit
# fit an MM regression model
fit <- lmrob(Y ~ ., data=coleman)
# perform cross-validation
cvFit(fit, data = coleman, y = coleman$Y, cost = rtmspe, 
    K = 5, R = 10, costArgs = list(trim = 0.1), seed = 1234)

## via model fitting function
# perform cross-validation
# note that the response is extracted from 'data' in 
# this example and does not have to be supplied
cvFit(lmrob, formula = Y ~ ., data = coleman, cost = rtmspe, 
    K = 5, R = 10, costArgs = list(trim = 0.1), seed = 1234)

## via function call
# set up function call
call <- call("lmrob", formula = Y ~ .)
# perform cross-validation
cvFit(call, data = coleman, y = coleman$Y, cost = rtmspe, 
    K = 5, R = 10, costArgs = list(trim = 0.1), seed = 1234)
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link{cvTool}}, \code{\link{cvSelect}},
  \code{\link{cvTuning}}, \code{\link{cvFolds}},
  \code{\link{cost}}
}
\keyword{utilities}

