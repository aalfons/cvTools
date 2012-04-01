\name{cvTuning}
\alias{cvTuning}
\alias{cvTuning.call}
\alias{cvTuning.function}
\alias{print.cvTuning}
\title{Cross-validation for tuning parameter selection}
\usage{
  cvTuning(object, ...)

  \method{cvTuning}{function} (object, formula,
    data = NULL, x = NULL, y, tuning = list(),
    args = list(), cost = rmspe, K = 5, R = 1,
    foldType = c("random", "consecutive", "interleaved"),
    folds = NULL, names = NULL, predictArgs = list(),
    costArgs = list(), selectBest = c("min", "hastie"),
    seFactor = 1, envir = parent.frame(), seed = NULL, ...)

  \method{cvTuning}{call} (object, data = NULL, x = NULL,
    y, tuning = list(), cost = rmspe, K = 5, R = 1,
    foldType = c("random", "consecutive", "interleaved"),
    folds = NULL, names = NULL, predictArgs = list(),
    costArgs = list(), selectBest = c("min", "hastie"),
    seFactor = 1, envir = parent.frame(), seed = NULL, ...)
}
\arguments{
  \item{object}{a function or an unevaluated function call
  for fitting a model (see \code{\link{call}} for the
  latter).}

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

  \item{tuning}{a list of arguments giving the tuning
  parameter values to be evaluated.  The names of the list
  components should thereby correspond to the argument
  names of the tuning parameters.  For each tuning
  parameter, a vector of values can be supplied.
  Cross-validation is then applied over the grid of all
  possible combinations of tuning parameter values.}

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

  \item{selectBest}{a character string specifying a
  criterion for selecting the best model.  Possible values
  are \code{"min"} (the default) or \code{"hastie"}.  The
  former selects the model with the smallest prediction
  error.  The latter is useful for models with a tuning
  parameter controlling the complexity of the model (e.g.,
  penalized regression).  It selects the most parsimonious
  model whose prediction error is no larger than
  \code{seFactor} standard errors above the prediction
  error of the best overall model.  Note that the models
  are thereby assumed to be ordered from the most
  parsimonious one to the most complex one.  In particular
  a one-standard-error rule is frequently applied.}

  \item{seFactor}{a numeric value giving a multiplication
  factor of the standard error for the selection of the
  best model.  This is ignored if \code{selectBest} is
  \code{"min"}.}

  \item{envir}{the \code{\link{environment}} in which to
  evaluate the function call for fitting the models (see
  \code{\link{eval}}).}

  \item{seed}{optional initial seed for the random number
  generator (see \code{\link{.Random.seed}}).}

  \item{\dots}{additional arguments to be passed down.}
}
\value{
  If \code{tuning} is an empty list, \code{\link{cvFit}} is
  called to return an object of class \code{"cv"}.

  Otherwise an object of class \code{"cvTuning"} (which
  inherits from class \code{"cvSelect"}) with the following
  components is returned:

  \item{n}{an integer giving the number of observations.}

  \item{K}{an integer giving the number of folds.}

  \item{R}{an integer giving the number of replications.}

  \item{tuning}{a data frame containing the grid of tuning
  parameter values for which the prediction error was
  estimated.}

  \item{best}{an integer vector giving the indices of the
  optimal combinations of tuning parameters.}

  \item{cv}{a data frame containing the estimated
  prediction errors for all combinations of tuning
  parameter values.  For repeated cross-validation, those
  are average values over all replications.}

  \item{se}{a data frame containing the estimated standard
  errors of the prediction loss for all combinations of
  tuning parameter values.}

  \item{selectBest}{a character string specifying the
  criterion used for selecting the best model.}

  \item{seFactor}{a numeric value giving the multiplication
  factor of the standard error used for the selection of
  the best model.}

  \item{reps}{a data frame containing the estimated
  prediction errors from all replications for all
  combinations of tuning parameter values.  This is only
  returned for repeated cross-validation.}

  \item{seed}{the seed of the random number generator
  before cross-validation was performed.}

  \item{call}{the matched function call.}
}
\description{
  Select tuning parameters of a model by estimating the
  respective prediction errors via (repeated) \eqn{K}-fold
  cross-validation.  It is thereby possible to supply a
  model fitting function or an unevaluated function call to
  a model fitting function.
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
  even if only one is actually used.  For the \code{call}
  method, which does not have a \code{formula} argument, a
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
\note{
  The same cross-validation folds are used for all
  combinations of tuning parameter values for maximum
  comparability.
}
\examples{
library("robustbase")
data("coleman")

## evaluate MM regression models tuned for 85\% and 95\% efficiency
tuning <- list(tuning.psi = c(3.443689, 4.685061))

## via model fitting function
# perform cross-validation
# note that the response is extracted from 'data' in 
# this example and does not have to be supplied
cvTuning(lmrob, formula = Y ~ ., data = coleman, tuning = tuning, 
    cost = rtmspe, K = 5, R = 10, costArgs = list(trim = 0.1), 
    seed = 1234)

## via function call
# set up function call
call <- call("lmrob", formula = Y ~ .)
# perform cross-validation
cvTuning(call, data = coleman, y = coleman$Y, tuning = tuning, 
    cost = rtmspe, K = 5, R = 10, costArgs = list(trim = 0.1), 
    seed = 1234)
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
  \code{\link{cvTool}}, \code{\link{cvFit}},
  \code{\link{cvSelect}}, \code{\link{cvFolds}},
  \code{\link{cost}}
}
\keyword{utilities}

