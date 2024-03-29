% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/colour.R
\name{color.plot.data.frame}
\alias{color.plot.data.frame}
\alias{color.plot.formula}
\title{Plot cases as colored points}
\usage{
\method{color.plot}{data.frame}(x, z, zlab = NULL, xlab = NULL, ylab = NULL, labels = F, ...)
}
\arguments{
\item{labels}{If NULL, cases are plotted as points.  If T, cases are plotted
as labels, according to \code{rownames}.}

\item{...}{Extra arguments passed to \code{\link{color.plot.default}}.}

\item{data}{a data frame.}

\item{formula}{a formula specifying a response and two predictors from
\code{data}}
}
\description{
Like \code{\link{plot}} and \code{\link{text.plot}} but colors according to
the response variable.
}
\details{
Calls \code{\link{color.plot.default}} with \code{x} as the first predictor
in \code{data}, \code{y} as the second predictor, and \code{z} as the
response.  To get a different predictor/response breakdown than the default,
use \code{color.plot(formula, x, ...)}, which is shorthand for
\code{color.plot(model.frame(formula, x), ...)}.

Each case is plotted with a color determined by the response.  If the
response is a factor, each factor level is in a different color.  If the
response is numeric, then it is color-coded by assigning each quantile a
different color.
}
\examples{

data(iris)
color.plot(iris)
color.plot(Species ~ Petal.Length + Petal.Width, iris)
#color.plot(Species ~ Petal.Length, iris)
#color.plot(Species ~ Petal.Length, iris,jitter=T)
color.plot(iris, col=1)
color.plot(iris, col=c(1,2))

data(state)
x <- data.frame(state.x77)
color.plot(Murder ~ Frost + Illiteracy, x, labels=T, cex=0.5)

}
\seealso{
\code{\link{color.plot.default}}
}
\author{
Tom Minka
}
