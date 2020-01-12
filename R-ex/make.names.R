### Name: make.names
### Title: Make Syntactically Valid Names
### Aliases: make.names
### Keywords: character

### ** Examples

make.names(c("a and b", "a_and_b"), unique=TRUE)
# "a.and.b"  "a.and.b1"

data(state)
state.name[make.names(state.name) != state.name]# those 10 with a space



