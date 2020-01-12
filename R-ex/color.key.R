### Name: color.key
### Title: Add a thin key to the top of a plot
### Aliases: color.key


### ** Examples

data(iris)
y = as.numeric(iris$Species)
plot(Sepal.Width ~ Sepal.Length, iris,col=y,pch=y)
color.key(1:3,1:3,levels(iris$Species))



