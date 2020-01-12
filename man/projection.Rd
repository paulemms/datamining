\author{Tom Minka}
\name{pca}
\alias{pca}
\title{Principal Component Analysis}
\description{Computes a projection matrix that preserves spread.}
\usage{
pca(x,k=1)
}
\arguments{
  \item{x}{a data frame.}
  \item{k}{the number of dimensions to project down to.}
}
\details{
  The projection is chosen so that the projected data is as "spread out"
  as possible, measured according to the determinant of the covariance
  matrix.
  This turns out to be the top \code{k} eigenvectors of the data
  covariance matrix.

  The projection is "stabilized" so that small changes in the data do
  not cause sign flips in the projection.
}
\value{
  A matrix with named rows matching the numeric columns of \code{x} and
  columns named \code{h1}, ..., \code{hk}.
  Each column denotes a new dimension to be obtained as a linear combination
  of the numeric variables in \code{x}.
}
\seealso{\code{\link{project}}, \code{\link{projection}},
  \code{pca} in the \code{multiv} package}
\examples{
data(Housing)
w = pca(HousingT,k=2)
plot(project(HousingT,w),asp=1)
plot.axes(w)
}
