### Name: color.plot.loess
### Title: Contour plot of a regression surface
### Aliases: color.plot.loess


### ** Examples

data(Housing)
fit = loess(Price ~ Rooms + Low.Status, Housing)
color.plot(fit)



