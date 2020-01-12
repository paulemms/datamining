### Name: colors
### Title: Color schemes
### Aliases: default.colors default.colors.w YR.colors YlGnBu.colors
###   OrRd.colors gray.colors RYB.colors BrBg.colors RC.colors GM.colors


### ** Examples

data(Housing)
color.plot(Price ~ Rooms + Low.Status, Housing, bg=gray(0.5),
           color.palette=YlGnBu.colors)
color.plot(Price ~ Rooms + Low.Status, Housing, bg=gray(0.5),
           color.palette=YR.colors)
color.plot(Price ~ Rooms + Low.Status, Housing, bg=gray(0.5),
           color.palette=RYB.colors,nlevels=5)

# also see examples for color.cone



