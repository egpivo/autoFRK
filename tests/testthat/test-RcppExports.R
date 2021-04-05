tolerance <- 1e-6
matrix <- matrix(c(1, 2, 2, 1), 2, 2)
result <- eigenDecompose(matrix)
test_that("Eigen Decomposition", {
  expect_lte(norm(result$vector, "F") - 1.414214, tolerance)
  expect_equal(sum(result$value), 2)
})
