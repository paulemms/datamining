### Name: state.pop
### Title: U.S. state population 1790-1990



### ** Examples

dotchart(sort(log(state.pop["1990",])))

# distribution of growth rates over time
rates <- diff(log(state.pop))
boxplot(as.data.frame(t(rates)))



