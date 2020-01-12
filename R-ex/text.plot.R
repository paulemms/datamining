### Name: text.plot
### Title: Make a plot of text labels
### Aliases: text.plot text.plot.default text.plot.data.frame
###   text.plot.formula


### ** Examples

data(state)
x <- data.frame(state.x77)
text.plot(x$Frost, x$HS.Grad, rownames(x))
# same thing, using text.plot.formula
text.plot(HS.Grad ~ Frost, x)
# notice how the limits change
text.plot(HS.Grad ~ Frost, x, srt=45)
text.plot(HS.Grad ~ Frost, x, srt=90)



