library(MASS)
utils::data(Boston, envir = environment())
names(Boston) <- c("Crime","Zoned","Industry","River","Pollution","Rooms","Old","Distance","Highway","Tax","Student.Teacher.Ratio","White","Low.Status","Price")
Housing <- Boston
Housing$White <- NULL
Housing$Zoned <- NULL
Housing$River <- NULL

HousingT = Housing
i <- c("Crime","Distance","Pollution")
HousingT[,i] <- log(HousingT[,i])
i <- c("Low.Status")
HousingT[,i] <- sqrt(HousingT[,i])

# scales the numeric columns, leaving rest alone
scale_df <- function(x) {
  i <- sapply(x, is.numeric)
  for(j in which(i)) {
    x[,j] <- scale.default(x[,j])
  }
  x
}
HousingT <- scale_df(HousingT)
