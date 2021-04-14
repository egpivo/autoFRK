set.seed(1234)
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

grid <- seq(0, 1, l = 30)
z <- sample(30, 10)

test_that("nc for LkrigInfo", {
  expect_equal(setNC(z, grid, 1), 4)
  expect_equal(setNC(z, grid, 2), 4)
})

two_dim_knots <- subKnot(matrix(grid, ncol = 2), 2)
true_two_dim_knots <- matrix(c(0, 0.3448276, 0.5172414, 0.8620690), 2, 2)
test_that("Sample knots", {
  expect_equal(subKnot(z, 4), c(5, 12, 22, 30))
  expect_lte(norm(two_dim_knots - true_two_dim_knots, "F"), tolerance)
})

a <- matrix(1)
attr(a, "UZ") <- matrix(2)
test_that("Set mrts object to matrix class", {
  expect_equal(as.matrix.mrts(a), matrix(1))
})


test_that("Sparse matrix", {
  expect_error(toSparseMatrix(1), "Wrong format for toSparseMatrix")
  expect_message(toSparseMatrix(spam(0, 10, 10), TRUE), "The input is already a sparse matrix")
  expect_true(is.spam(toSparseMatrix(matrix(c(0, 0, 0, 1), 2, 2))))
})


mrts_message <- capture_output(print.mrts(mrts(1, 2)), print = TRUE)
test_that("Print mrts", {
  expect_error(print.mrts(1), "Invalid object! Please enter an `mrts` object")
  expect_equal(mrts_message, "[1]   1 NaN")
})

FRK_message <- capture_output(print.FRK(autoFRK(Data = rnorm(10), loc = 1:10, maxK = 3)))
test_that("Print FRK", {
  expect_error(print.FRK(1), "Invalid object! Please enter an `FRK` object")
  expect_equal(FRK_message, "[1] \"a 10 by 2 mrts matrix\"")
})
