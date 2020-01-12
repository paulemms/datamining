### Name: Cars
### Title: Cars, Housing, States, and VA Deaths
### Aliases: Cars Housing States va.deaths


### ** Examples

data(States)
hist(States)
hist(StatesT)
w = pca(StatesT,2)
text.plot(project(StatesT,w),asp=1,cex=.6)
plot.axes(w)



