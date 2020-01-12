### Name: break.ward
### Title: Quantize by clustering
### Aliases: break.ward break.kmeans break.hclust


### ** Examples

x <- c(rnorm(700,-2,1.5),rnorm(300,3,0.5))
break.ward(x,2)
break.hclust(x,2,method="complete")
break.kmeans(x,2)

x <- c(rnorm(700,-2,0.5),rnorm(1000,2.5,1.5),rnorm(500,7,0.1))
break.ward(x,3)
break.hclust(x,3,method="complete")
break.kmeans(x,3)



