% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/glm.R
\name{linechart}
\alias{linechart}
\title{Linechart}
\usage{
linechart(
  y,
  se = NULL,
  xlab = NULL,
  ylab,
  effects = F,
  med = F,
  xscale = c("equal", "linear", "none"),
  ...
)
}
\arguments{
\item{y}{a named vector or matrix.}

\item{se}{a vector or matrix, same size as \code{y}, of error bounds.
Alternatively, \code{y} can be \code{list(y,se)}.}

\item{effects}{If \code{TRUE}, the columns are shifted to have zero mean.}

\item{med}{If \code{TRUE} and \code{effects=TRUE}, the columns are shifted
to have zero median.}

\item{xscale}{describes how the columns should be placed.}

\item{...}{additional arguments to \code{\link{labeled.curves}}.}
}
\description{
Plot each data row as a curve.
}
\details{
If \code{xscale="linear"}, the columns are placed to make each curve as
straight as possible.  If \code{xscale="equal"}, the columns are placed
similar to \code{"linear"} but with the constraint that they must be equally
spaced.  If \code{xscale="none"}, the columns are placed in the order that
they appear in the matrix.  This is automatic if \code{y} has ordered
columns (see \code{\link{dim.ordered}}).  If \code{se != NULL}, error bars
are drawn around each point.

Linecharts are a replacement for dotcharts and mosaics.
}
\examples{

# compare to a dotchart
data(VADeaths)
dotchart(VADeaths, main = "Death Rates in Virginia - 1940")
dimOrdered(VADeaths)[2] = F
linechart(VADeaths)
linechart(t(VADeaths))

# compare to a mosaicplot
data(HairEyeColor)
x <- margin.table(HairEyeColor,c(1,2))
dimOrdered(x) = F
mosaicplot(x)
x = t(x)
col = c("brown","blue","red","green")
linechart(row.probs(x),color.pal=col)
linechart(row.probs(x,se=T),color.pal=col)
linechart(row.probs(x,se=T),jitter=0.02,color.pal=col)
mosaicplot(x)
linechart(row.probs(t(x),se=T))

data(blood)
dimOrdered(blood) = F
linechart(row.probs(blood,se=T))

data(antacids)
dimOrdered(antacids) = F
linechart(row.probs(antacids,se=T))
mosaicplot(t(antacids))

}
\seealso{
\code{\link{dotchart}}, \code{\link{mosaicplot}}
}
\author{
Tom Minka
}
