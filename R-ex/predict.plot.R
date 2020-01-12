### Name: predict.plot
### Title: Plot predictors versus response.
### Aliases: predict.plot predict.plot.data.frame predict.plot.formula


### ** Examples

data(Cars)
predict.plot(Price~.,CarsT)
fit = lm(Price~.,CarsT)
predict.plot(Price~.,CarsT,partial=fit)
# same thing using predict.plot.lm
predict.plot(fit,partial=T)



