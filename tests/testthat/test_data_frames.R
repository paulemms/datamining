test_that("Data frames are sorted correctly", {
  dfr <- data.frame(a = c(NA, 1, 2, NA), b = c(NA, NA, NA, NA), c = 1:4,
                    d = c(1:3, NA), e = c(NA, 2, NA, 3))
  output_dfr <- sort.data.frame(dfr)
  expect_equal(output_dfr, dfr[c(2, 4, 1, 3),,drop=FALSE])
})
