### Name: projection
### Title: Discriminative projection
### Aliases: projection


### ** Examples

# illustrate difference between (m,v,mv)
library(MASS)
m1 <- c(6,6)
v1 <- array(c(2,1.9,1.9,2),c(2,2))
#v1 <- array(c(1,0,0,1),c(2,2))
x1 <- mvrnorm(100,m1,v1)
m2 <- c(0,0)
v2 <- array(c(20,0,0,10),c(2,2))
x2 <- mvrnorm(300,m2,v2)
x = as.data.frame(rbind(x1,x2))
y = factor(c(rep(1,nrow(x1)),rep(2,nrow(x2))))
plot(x[,1],x[,2],col=1,xlab="",ylab="",asp=1)
points(x2[,1],x2[,2],col=2)
w = projection(x,y,type="m")
abline(0,w[2]/w[1],col=3)
w = projection(x,y,type="v")
abline(0,w[2]/w[1],col=4)
w = projection(x,y,type="mv")
abline(0,w[2]/w[1],col=5)
my.legend(1,c("m","v","mv"),col=3:5,lty=1)

# regression projection
x1 <- 2*runif(200)-1
x2 <- 2*runif(200)-1
y <- x1^2/2 + x2^2
x <- data.frame(x1,x2)
color.plot(x[,1],x[,2],y)
w = projection(x,y)
abline(0,w[2]/w[1],col=4)



