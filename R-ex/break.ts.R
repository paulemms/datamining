### Name: break.ts
### Title: Change-point analysis by clustering
### Aliases: break.ts


### ** Examples

library(ts)
data(LakeHuron)
# single major change
break.ts(LakeHuron,2)
# merging trace suggests n=6 is also interesting:
break.ts(LakeHuron,6)
# interesting oscillation

data(treering)
break.ts(treering[1:500],9,same=T)
break.ts(treering[1:100],7,same=T)
# interesting multiscale structure

x <- c(rnorm(100),rnorm(300)*3,rnorm(200)*2)
b <- break.ts(x,3,same=F)
plot(x,type="l")
plot.breaks(b)



