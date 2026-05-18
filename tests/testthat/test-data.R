test_that("gold_msci_returns contains return dates", {
  data("gold_msci_returns", package = "bekkIRF", envir = environment())

  expect_equal(dim(gold_msci_returns), c(4920L, 2L))
  expect_equal(colnames(gold_msci_returns), c("msci.Ret", "gold.Ret"))
  expect_equal(head(rownames(gold_msci_returns), 1L), "2007-01-03")
  expect_equal(tail(rownames(gold_msci_returns), 1L), "2025-11-28")

  dates <- attr(gold_msci_returns, "dates")
  expect_s3_class(dates, "Date")
  expect_equal(length(dates), nrow(gold_msci_returns))
  expect_equal(range(dates), as.Date(c("2007-01-03", "2025-11-28")))
  expect_equal(attr(gold_msci_returns, "price_date_range"), as.Date(c("2007-01-02", "2025-11-28")))
  expect_equal(attr(gold_msci_returns, "return_date_range"), as.Date(c("2007-01-03", "2025-11-28")))
})
