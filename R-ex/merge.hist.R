### Name: merge.hist
### Title: Merge histogram bins
### Aliases: merge.hist


### ** Examples

x <- c(rnorm(100,-2,0.5),rnorm(100,2,0.5))
b <- seq(-4,4,by=0.25)
merge.hist(x,b,10)
# according to the merging trace, n=5 and n=11 are most interesting.

x <- runif(1000)
b <- seq(0,1,by=0.05)
merge.hist(x,b,10)
# according to the merging trace, n=6 and n=9 are most interesting.
# because the data is uniform, there should only be one bin,
# but chance deviations in density prevent this.
# a multiple comparisons correction in merge.hist may fix this.



