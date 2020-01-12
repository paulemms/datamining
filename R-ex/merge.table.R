### Name: merge.table
### Title: Table merging
### Aliases: merge.table


### ** Examples

i <- factor(c(1,1,2,2,3,3,4,4))
j <- factor(c(3,4,3,4,1,2,1,2))
x <- table(i,j)
merge.table(x,c(2,2))

i <- factor(c(1,1,3,3,2,2,4,4))
j <- factor(c(2,4,2,4,1,3,1,3))
x <- table(i,j)
merge.table(x,c(2,2))

# one ordered dimension
data(education)
merge.table(education,c(3,2))

data(occupation)
merge.table(occupation,c(3,4))



