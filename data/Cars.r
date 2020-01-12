library(MASS)
utils::data(Cars93, envir = environment())
Cars = Cars93
rownames(Cars) = Cars$Make
Cars$Model <- Cars$Manufacturer <- Cars$Make <- NULL
Cars$Max.Price <- Cars$Min.Price <- NULL
Cars$RPM <- Cars$Rev.per.mile <- Cars$MPG.city <- Cars$Fuel.tank.capacity <- NULL
# these have missing data
Cars$Luggage.room <- Cars$Rear.seat.room <- NULL
# non-numeric
Cars$DriveTrain <- Cars$Cylinders <- Cars$Man.trans.avail <- NULL
Cars$AirBags <- Cars$Origin <- NULL
# specify order of levels
Cars$Type <- factor(Cars$Type, c("Small","Sporty","Compact","Midsize","Large","Van"))
Cars <- Cars[order(Cars$Type),]
Cars = stats::model.frame(Type~.,Cars)
# bugfix
Cars["Pontiac Bonneville","Length"] <- 200

CarsT = Cars
i <- c("Price","MPG.highway","EngineSize")
CarsT[,i] <- log(CarsT[,i])
i <- "Horsepower"
CarsT[,i] <- sqrt(CarsT[,i])
ix <- sapply(CarsT, is.numeric)
CarsT[,ix] = lapply(CarsT[,ix], scale)
