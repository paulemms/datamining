### Name: vrml.plot3d
### Title: View points in 3D
### Aliases: vrml.plot3d vrml.plot3d.default vrml.plot3d.data.frame
###   vrml.plot3d.formula


### ** Examples

data(Housing)
w = pca(HousingT,k=3)
x = project(HousingT,w)
plot.new()
vrml.plot3d(x)



