### Name: merge.factor
### Title: Merge factor levels
### Aliases: merge.factor


### ** Examples

n <- 20
x <- c(rnorm(n)+1, rnorm(n)+2, rnorm(n)*4+2)
f <- gl(3,n)
levels(f) <- c("a","b","c")
merge.factor(f,x,2,same.var=T)
merge.factor(f,x,2,same.var=F)

# an ordered factor
data(va.deaths)
merge.factor(va.deaths$Age,va.deaths$Rate,2)



