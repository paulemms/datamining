### Name: pca
### Title: Principal Component Analysis
### Aliases: pca


### ** Examples

data(Housing)
w = pca(HousingT,k=2)
plot(project(HousingT,w),asp=1)
plot.axes(w)



