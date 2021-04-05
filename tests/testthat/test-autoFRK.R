set.seed(0)
n <- 150
s <- 5
grid1 <- grid2 <- seq(0, 1, l = 30)
grids <- expand.grid(grid1, grid2)
fn <- matrix(0, 900, 2)
fn[, 1] <- cos(sqrt((grids[, 1] - 0)^2 + (grids[, 2] - 1)^2) * pi)
fn[, 2] <- cos(sqrt((grids[, 1] - 0.75)^2 + (grids[, 2] - 0.25)^2) * 2 * pi)

#### single realization simulation example
w <- c(rnorm(1, sd = 5), rnorm(1, sd = 3))
y <- fn %*% w
obs <- sample(900, n)
z <- y[obs] + rnorm(n) * sqrt(s)
X <- grids[obs, ]

#### Example1
single_realization_object <- autoFRK(Data = z, loc = X, maxK = 15)
yhat_example1 <- predict(single_realization_object, newloc = grids)

#### Example2
G <- mrts(X, 15)
Gpred <- predict(G, newx = grids)
user_defined_basis_object <- autoFRK(Data = z, loc = X, G = G)
yhat_example2 <- predict(user_defined_basis_object, newloc = grids, basis = Gpred)

tolerance <- 1e-4
# Test
test_that("Automatic selection and prediction", {
  expect_lte(
    abs(mean(yhat_example1$pred.value) + 3.029648), tolerance
  )
  expect_lte(abs(sum(yhat_example1$pred.value) + 2726.683), tolerance)
  expect_null(yhat_example1$se)
  expect_equal(length(yhat_example1$pred.value), 900)
})

test_that("User-specified basis functionsn", {
  expect_lte(
    abs(mean(yhat_example2$pred.value) + 3.024821), tolerance
  )
  expect_lte(abs(sum(yhat_example2$pred.value) + 2722.339), tolerance)
  expect_null(yhat_example2$se)
  expect_equal(length(yhat_example2$pred.value), 900)
})