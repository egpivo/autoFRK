tolerance <- 1e-6

test_r <- matrix(c(0.1, 0.2), 1, 2)
test_that("Wendland function", {
  expect_error(
    wendland(c(-.3, 0.1, -0.3)),
    "Invalid values:   -0.3  -0.3 "
  )
  expect_equal(sum(wendland(test_r)), 2)
})

matrix <- matrix(c(1, 2, 2, 1), 2, 2)
result <- eigenDecomposeInDecreasingOrder(matrix)
test_that("Eigen Decomposition", {
  expect_lte(norm(result$vector, "F") - 1.414214, tolerance)
  expect_equal(sum(result$value[1]), 3)
})

test_that("Convert value to bytes", {
  expect_equal(toBytes(c("1", "gb")), 1073741824)
})

test_that("Fetch system RAM", {
  expect_true(any("numeric" %in% class(fetchSystemRam(R.Version()$platform))))
})

test_that("Remove white space", {
  expect_equal(removeWhitespace(" test test "), "test test")
})
