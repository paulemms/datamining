% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mining.R
\name{break_ward}
\alias{break_ward}
\alias{break_kmeans}
\alias{break_hclust}
\title{Quantize by clustering}
\usage{
break_ward(x, n = 2, plot = T)
}
\arguments{
\item{x}{a numerical vector}

\item{n}{the desired number of bins}

\item{plot}{If TRUE, a histogram with break lines is plotted
(\code{\link{hist.hclust}} or \code{\link{plot_breaks}}). For
\code{break_ward} and \code{break_hclust}, also shows a merging trace
(\code{\link{plot_hclust_trace}}).}

\item{method}{argument given to \code{\link{hclust}}}
}
\value{
A break vector, suitable for use in \code{\link{cut}} or
\code{\link{hist}}.
}
\description{
Quantize a one-dimensional variable by calling a clustering routine.
}
\details{
These are convenience routines which simply call the appropriate clustering
routine (\code{\link{ward}}, \code{\link{hclust}}, or \code{\link{kmeans}}),
convert the output to a break vector, and make plots.
}
\examples{

x <- c(rnorm(700,-2,1.5),rnorm(300,3,0.5))
break_ward(x,2)
break_hclust(x,2,method="complete")
break_kmeans(x,2)

x <- c(rnorm(700,-2,0.5),rnorm(1000,2.5,1.5),rnorm(500,7,0.1))
break_ward(x,3)
break_hclust(x,3,method="complete")
break_kmeans(x,3)

}
\author{
Tom Minka
}
