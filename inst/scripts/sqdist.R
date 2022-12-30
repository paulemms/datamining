# improved performance of sqdist2 using outer approach
# but what does A do in sqdist2?
library(microbenchmark)

set.seed(1)
A=matrix(rnorm(1000,5,50),ncol=5)
B=matrix(rnorm(10000,0,50),ncol=5)

sqdist2 <- function(x,y=x) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  xmag <- rowSums(x * x)
  ymag <- rowSums(y * y)
  #browser()
  d <- rep(ymag, rep.int(nrow(x), length(ymag))) +
    rep(xmag, times = nrow(y)) - 2*(x %*% t(y))
  if(identical(x,y)) diag(d) = 0
  # fix numeric errors
  d[d < 0] <- 0
  d
}

b=microbenchmark(times=100,
                 sqdist(A, B),
                 sqdist2(A, B),
                 apply(B,1,function(x)colSums((t(A)-x)^2)),
                 outer(rowSums(A^2),rowSums(B^2),"+")-2*tcrossprod(A,B),
                 outer(rowSums(A^2),rowSums(B^2),"+")-2*A%*%t(B)
)

a=aggregate(b$time,list(b$expr),median)
a=a[order(a[,2]),]
writeLines(paste(sprintf("%.3f",a[,2]/min(a[,2])),gsub(" ","",a[,1])))

all(sqdist(A,B) == sqdist2(A,B))
