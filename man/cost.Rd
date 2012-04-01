\name{cost}
\alias{cost}
\alias{mape}
\alias{mspe}
\alias{rmspe}
\alias{rtmspe}
\alias{tmspe}
\title{Prediction loss}
\usage{
  mspe(y, yHat, includeSE = FALSE)

  rmspe(y, yHat, includeSE = FALSE)

  mape(y, yHat, includeSE = FALSE)

  tmspe(y, yHat, trim = 0.25, includeSE = FALSE)

  rtmspe(y, yHat, trim = 0.25, includeSE = FALSE)
}
\arguments{
  \item{y}{a numeric vector or matrix giving the observed
  values.}

  \item{yHat}{a numeric vector or matrix of the same
  dimensions as \code{y} giving the fitted values.}

  \item{trim}{a numeric value giving the trimming
  proportion (the default is 0.25).}

  \item{includeSE}{a logical indicating whether standard
  errors should be computed as well.}
}
\value{
  If standard errors are not requested, a numeric value
  giving the prediction loss is returned.

  Otherwise a list is returned, with the first component
  containing the prediction loss and the second component
  the corresponding standard error.
}
\description{
  Compute the prediction loss of a model.
}
\details{
  \code{mspe} and \code{rmspe} compute the mean squared
  prediction error and the root mean squared prediction
  error, respectively.  In addition, \code{mape} returns
  the mean absolute prediction error, which is somewhat
  more robust.

  Robust prediction loss based on trimming is implemented
  in \code{tmspe} and \code{rtmspe}.  To be more precise,
  \code{tmspe} computes the trimmed mean squared prediction
  error and \code{rtmspe} computes the root trimmed mean
  squared prediction error.  A proportion of the largest
  squared differences of the observed and fitted values are
  thereby trimmed.

  Standard errors can be requested via the \code{includeSE}
  argument.  Note that standard errors for \code{tmspe} are
  based on a winsorized standard deviation.  Furthermore,
  standard errors for \code{rmspe} and \code{rtmspe} are
  computed from the respective standard errors of
  \code{mspe} and \code{tmspe} via the delta method.
}
\examples{
# fit an MM-regression model
data("coleman")
fit <- lmrob(Y~., data=coleman)

# compute the prediction loss from the fitted values
# (hence the prediction loss is underestimated in this simple
# example since all observations are used to fit the model)
mspe(coleman$Y, predict(fit))
rmspe(coleman$Y, predict(fit))
mape(coleman$Y, predict(fit))
tmspe(coleman$Y, predict(fit), trim = 0.1)
rtmspe(coleman$Y, predict(fit), trim = 0.1)

# include standard error
mspe(coleman$Y, predict(fit), includeSE = TRUE)
rmspe(coleman$Y, predict(fit), includeSE = TRUE)
mape(coleman$Y, predict(fit), includeSE = TRUE)
tmspe(coleman$Y, predict(fit), trim = 0.1, includeSE = TRUE)
rtmspe(coleman$Y, predict(fit), trim = 0.1, includeSE = TRUE)
}
\author{
  Andreas Alfons
}
\references{
  Tukey, J.W. and McLaughlin, D.H. (1963) Less vulnerable
  confidence and significance procedures for location based
  on a single sample: Trimming/winsorization.
  \emph{Sankhya: The Indian Journal of Statistics, Series
  A}, \bold{25}(3), 331--352

  Oehlert, G.W. (1992) A note on the delta method.
  \emph{The American Statistician}, \bold{46}(1), 27--29.
}
\seealso{
  \code{\link{cvFit}}, \code{\link{cvTuning}}
}
\keyword{utilities}

