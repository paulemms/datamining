% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mining.R
\name{ward}
\alias{ward}
\title{Create a hierarchy by Ward's method}
\usage{
ward(
  x,
  n = rep(1, length(x)),
  s = rep(1, length(x)),
  sortx = TRUE,
  same.var = T
)
}
\arguments{
\item{x}{a numerical vector, or a list of vectors.}

\item{n}{if x is a vector of cluster means, n is the size of each cluster.}

\item{s}{if x is a vector of cluster means, s is the sum of squares in each
cluster.  only needed if \code{same.var=F}.}

\item{sortx}{if \code{sortx=F}, only clusters which are adjacent in \code{x}
can be merged.  Used by \code{\link{break.ts}}.}

\item{same.var}{if \code{same.var=T}, clusters are assumed to have the same
true variance, otherwise not.  This affects the cost function for merging.}
}
\value{
The same type of object returned by \code{\link{hclust}}.
}
\description{
Produces a hierarchical clustering of one-dimensional data via Ward's method.
}
\details{
Repeatedly merges clusters in order to minimize the clustering cost.
By default, it is the same as \code{hclust(method="ward")}.
If \code{same.var=T}, the cost is the sum of squares:
\deqn{sum_c sum_{i in c} (x_i - m_c)^2}
The incremental cost of merging clusters ca and cb is
\deqn{(n_a*n_b)/(n_a+n_b)*(m_a - m_b)^2}
It prefers to merge clusters which are small and have similar means.

If \code{same.var=F}, the cost is the sum of log-variances:
\deqn{sum_c n_c*log(1/n_c*sum_{i in c} (x_i - m_c)^2)}
It prefers to merge clusters which are small, have similar means,
and have similar variances.

If \code{x} is a list of vectors, each vector is assumed to be a
cluster.  \code{n} and \code{s} are computed for each cluster and
\code{x} is replaced by the cluster means.
Thus you can say \code{ward(split(x,f))} to cluster the data for different
factors.
}
\section{Bugs}{
 Because of the adjacency constraint used in implementation,
the clustering that results
from \code{sortx=T} and \code{same.var=F} may occasionally be suboptimal.
}

\examples{
x <- c(rnorm(700,-2,1.5),rnorm(300,3,0.5))
hc <- ward(x)
opar <- par(mfrow=c(2,1))
# use dev.new() in RStudio
plot_hclust_trace(hc)
hist(hc,x)
par(opar)

x <- c(rnorm(700,-2,0.5),rnorm(1000,2.5,1.5),rnorm(500,7,0.1))
hc <- ward(x)
opar <- par(mfrow=c(2,1))
plot_hclust_trace(hc)
hist(hc,x)
par(opar)

data(OrchardSprays)
x <- OrchardSprays$decrease
f <- factor(OrchardSprays$treatment)
# shuffle levels
#lev <- levels(OrchardSprays$treatment)
#f <- factor(OrchardSprays$treatment,levels=sample(lev))
hc <- ward(split(x,f))
# is equivalent to:
#n <- tapply(x,f,length)
#m <- tapply(x,f,mean)
#s <- tapply(x,f,var)*n
#hc <- ward(m,n,s)
boxplot(hc,split(x,f))
}
\seealso{
\code{\link{hclust}},
\code{\link{plot_hclust_trace}},
\code{\link{hist.hclust}},
\code{\link{boxplot.hclust}},
\code{\link{break_ward}},
\code{\link{break.ts}},
\code{\link{merge_factor}}
}
\author{
Tom Minka
}
