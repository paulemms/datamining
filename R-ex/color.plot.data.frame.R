### Name: color.plot.data.frame
### Title: Plot cases as colored points
### Aliases: color.plot.data.frame color.plot.formula


### ** Examples

data(iris)
color.plot(iris)
color.plot(Species ~ Petal.Length + Petal.Width, iris)
color.plot(Species ~ Petal.Length, iris)
color.plot(Species ~ Petal.Length, iris,jitter=T)
color.plot(iris, col=1)
color.plot(iris, col=c(1,2))

data(state)
x <- data.frame(state.x77)
color.plot(Murder ~ Frost + Illiteracy, x, labels=T, cex=0.5)



