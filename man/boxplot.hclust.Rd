% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ts.R
\name{boxplot.hclust}
\alias{boxplot.hclust}
\title{Boxplot with hierarchical cluster breaks}
\usage{
\method{boxplot}{hclust}(hc, x, k = 2:5, col = "bisque", ...)
}
\arguments{
\item{x}{the list of vectors that was clustered to produce \code{h}
(typically via \code{\link{ward}})}

\item{k}{a vector of the cluster cardinalities to plot}

\item{...}{arguments passed to \code{boxplot}}

\item{h}{an \code{hclust} object}
}
\value{
A boxplot of \code{x} is shown with blue lines cutting the x-axis.
The tallest lines correspond to divisions made at the top of the hierarchy.
By reading top to bottom, you can see how each cluster is subdivided. This
can be far more illuminating than a plot of the hierarchy as a tree.
}
\description{
A representation of a hierarchical clustering of predefined groups
}
\seealso{
\code{\link{hist.hclust}}, \code{\link{ward}},
\code{\link{break_ward}}
}
\author{
Tom Minka
}
