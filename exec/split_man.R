rm(list = ls())
ln <- readLines("original_man/mining.Rd")
p <- paste(ln, collapse = "\n")
ss <- strsplit(p, "\\eof", fixed = TRUE)[[1]]
r <- regexec("name(.*)", ss)
n <- substr(ss[1], r[[1]][1], r[[1]][2])
#write(ss[1], "1.md")

regmatches(ss, regexec("name\\{(.*?)\\}", ss))

