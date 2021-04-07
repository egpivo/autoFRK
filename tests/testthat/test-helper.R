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

test_that("Is an object diagonal", {
  expect_true(isDiagonal(1))
  expect_true(isDiagonal(diag(2)))
  expect_false(isDiagonal("ss"))
  expect_false(isDiagonal(matrix))
})

grid <-seq(0, 1, l = 30)
z <- sample(30, 10)

test_that("nc for LkrigInfo", {
  expect_equal(setNC(z, grid, 1), 4)
  expect_equal(setNC(z, grid, 2), 4)
})

test_that("Sample knots", {
  expect_equal(subKnot(z, 4), c(3, 8, 19, 30))
})

a <- matrix(1)
attr(a, "UZ") <- matrix(2)
test_that("Set mrts object to matrix class", {
  expect_equal(as.matrix.mrts(a), matrix(1))
})


