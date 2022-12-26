#utils::data(state, envir = environment()) # stops install with hard-coded path
States <- data.frame(datasets::state.x77)
names(States)[5] <- "Homocide"
States$Density <- States$Population/States$Area
States$Population <- States$Area <- NULL

StatesT = States
i <- c("Illiteracy","Density")
StatesT[,i] <- log(StatesT[,i])
StatesT <- scale(StatesT)
