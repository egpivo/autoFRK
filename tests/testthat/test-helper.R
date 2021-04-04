test_r <- matrix(c(0.1, 0.2), 1, 2)
tol <- 1e-6
test_that("Wendland function", {
  expect_error(
    wendland(c(-.3, 0.1, -0.3)),
    "Invalid values:   -0.3  -0.3 "
  )
  expect_equal(sum(wendland(test_r)), 2)
})
