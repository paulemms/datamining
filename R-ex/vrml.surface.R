### Name: vrml.surface
### Title: View a surface in 3D
### Aliases: vrml.surface vrml.surface.default vrml.surface.loess


### ** Examples

#  The Obligatory Mathematical surface.
#  Rotated sinc function.
x <- seq(-10, 10, length= 30)
y <- x
f <- function(x,y) { r <- sqrt(x^2+y^2); 10 * sin(r)/r }
z <- outer(x, y, f)
z[is.na(z)] <- 1
vrml.surface(x,y,z)
vrml.surface(x,y,z,border=T)



