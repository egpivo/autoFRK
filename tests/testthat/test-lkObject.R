set.seed(0)
n <- 100
s <- 5
grid1 <- grid2 <- seq(0, 1, l = 30)
location <- expand.grid(grid1, grid2)
z <- rnorm(n)

dummy_array <- array(1:n, s)

# Test
test_that("NC with different levels", {
  expect_equal(setNC(z, location, 0), 9)  
  expect_equal(setNC(z, location, 1), 10)
  expect_equal(setNC(z, location, 2), 4)
})

test_that("shiftArray with different shifts", {
  expect_equal(shiftArray(dummy_array, c(-1, 0, 1))[1], 2)  
  expect_equal(shiftArray(dummy_array, 1)[1], 5)
  expect_error(shiftArray(dummy_array, 10),
               "shift exceeds array dimensions")
})