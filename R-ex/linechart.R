### Name: linechart
### Title: Linechart
### Aliases: linechart


### ** Examples

# compare to a dotchart
data(VADeaths)
dotchart(VADeaths, main = "Death Rates in Virginia - 1940")
dimOrdered(VADeaths)[2] = F
linechart(VADeaths)
linechart(t(VADeaths))

# compare to a mosaicplot
data(HairEyeColor)
x <- margin.table(HairEyeColor,c(1,2))
dimOrdered(x) = F
mosaicplot(x)
x = t(x)
col = c("brown","blue","red","green")
linechart(row.probs(x),color.pal=col)
linechart(row.probs(x,se=T),color.pal=col)
linechart(row.probs(x,se=T),jitter=0.02,color.pal=col)
mosaicplot(x)
linechart(row.probs(t(x),se=T))

data(blood)
dimOrdered(blood) = F
linechart(row.probs(blood,se=T))

data(antacids)
dimOrdered(antacids) = F
linechart(row.probs(antacids,se=T))
mosaicplot(t(antacids))



