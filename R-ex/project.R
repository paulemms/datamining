### Name: project
### Title: Project data into fewer dimensions
### Aliases: project


### ** Examples

data(iris)
w = projection(iris,k=2)
# w only involves the continuous attributes
# the new variables are h1 and h2
x = project(iris,w)
color.plot(x)
plot.axes(w)



