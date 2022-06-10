set.seed(0)
n <- 100
s <- 5
grid1 <- grid2 <- seq(0, 1, l = 30)
location <- expand.grid(grid1, grid2)
z <- rnorm(n)

# Test
test_that("NC with different levels", {
  expect_equal(setNC(z, location, 0), 9)  
  expect_equal(setNC(z, location, 1), 10)
  expect_equal(setNC(z, location, 2), 4)
})