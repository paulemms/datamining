# library(Rd2roxygen)
# rm(list = ls())
# ln <- readLines("original_man/mining.Rd")
# p <- paste(ln, collapse = "\n")
# ss <- strsplit(p, "\\eof", fixed = TRUE)[[1]]
# r <- regmatches(ss, regexec("name\\{(.*?)\\}", ss))
# names(ss) <- sapply(r, function(x) x[2])
#
# o <- lapply(names(ss), function(x) writeLines(ss[x], con = paste0("man/", x, ".Rd")))
# Rd2roxygen(".", "datamining")
