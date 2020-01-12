### Name: ward
### Title: Create a hierarchy by Ward's method
### Aliases: ward


### ** Examples

x <- c(rnorm(700,-2,1.5),rnorm(300,3,0.5))
hc <- ward(x)
opar <- par(mfrow=c(2,1))
plot.hclust.trace(hc)
hist.hclust(hc,x)
par(opar)

x <- c(rnorm(700,-2,0.5),rnorm(1000,2.5,1.5),rnorm(500,7,0.1))
hc <- ward(x)
opar <- par(mfrow=c(2,1))
plot.hclust.trace(hc)
hist.hclust(hc,x)
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
boxplot.hclust(hc,split(x,f))



