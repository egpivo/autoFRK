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

test_that("Calculate the log determinant for the likelihood use", {
  expect_lte(abs(calculateLogDeterminant(diag(2), matrix(1:6, 2, 3), 3) - 4.75359), tolerance)
})

test_that("Convert value to bytes", {
  expect_equal(toBytes(c("1", "gb")), 1073741824)
  expect_equal(toBytes(c("1", "KBytes")), 1)
  expect_equal(toBytes(c("1", "megabytes")), 1048576)
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
  expect_true(isDiagonal(spam(0, 10, 10)))
  expect_false(isDiagonal("ss"))
  expect_false(isDiagonal(matrix))
})

set.seed(1234)
grid <- seq(0, 1, l = 30)
z <- sample(30, 10)
test_that("nc for LkrigInfo", {
  expect_equal(setNC(z, grid, 1), 4)
  expect_equal(setNC(z, grid, 2), 4)
})

two_dim_knots_example_1 <- subKnot(matrix(grid, ncol = 2), 2)
true_two_dim_knots_example_1 <- matrix(c(0, 0.3448276, 0.5172414, 0.8620690), 2, 2)
two_dim_knots_example_2 <- subKnot(matrix(grid, ncol = 2), nknot = 10, xrng = matrix(c(0, 0.1, 0, 0.1), 2))
true_sum_of_two_dim_knots_example_2 <- 13
test_that("Sample knots", {
  expect_equal(subKnot(z, 4), c(5, 12, 22, 30))
  expect_lte(norm(two_dim_knots_example_1 - true_two_dim_knots_example_1, "F"), tolerance)
  expect_equal(sum(two_dim_knots_example_2), true_sum_of_two_dim_knots_example_2)
})

a <- matrix(1)
attr(a, "UZ") <- matrix(2)
test_that("Set mrts object to matrix class", {
  expect_equal(as.matrix.mrts(a), matrix(1))
})

test_that("Fetch non-zero indeces", {
  expect_error(fetchNonZeroIndexs(1), "Wrong matrix format, but got numeric")
  expect_equal(fetchNonZeroIndexs(matrix(c(0, 1, 1, 0), 2, 2)), c(2, 3))
})

test_that("Sparse matrix", {
  expect_error(toSparseMatrix(1), "Wrong format for toSparseMatrix")
  expect_message(toSparseMatrix(spam(0, 10, 10), TRUE), "The input is already a sparse matrix")
  expect_true(is.spam(toSparseMatrix(matrix(c(0, 0, 0, 1), 2, 2))))
  expect_true(is.spam(toSparseMatrix(matrix(rnorm(100), 1, 100))))
  expect_true(is.spam(toSparseMatrix(matrix(rnorm(100), 100, 1))))
  expect_true(is.spam(toSparseMatrix(matrix(0, 100, 100))))
  expect_true(is.spam(toSparseMatrix(data.frame(1))))
})

R <- matrix(c(1, 2, 2, 1), 2)
L <- matrix(c(0.1, 0, 0, 0.1), 2)
z <- c(0, 1)
true_matrix <- matrix(c(0.6711635, -0.3389375), 1, 2)
test_that("Interanl matrix calculation function", {
  expect_lte(sum(ZinvC(R, L, z) - true_matrix), tolerance)
  expect_lte(sum(invCz(R, L, z) - t(true_matrix)), tolerance)
})

mrts_message <- capture_output(print.mrts(mrts(1, 2)), print = TRUE)
test_that("Print mrts", {
  expect_error(print.mrts(1), "Invalid object! Please enter an `mrts` object")
  expect_equal(mrts_message, "[1]   1 NaN")
  expect_equal(sum(print.mrts(mrts(matrix(1:10, 2), 6))), 2)
})

set.seed(1234)
FRK_message <- capture_output(print.FRK(autoFRK(Data = rnorm(10), loc = 1:10, maxK = 3)))
test_that("Print FRK", {
  expect_error(print.FRK(1), "Invalid object! Please enter an `FRK` object")
  expect_equal(FRK_message, "[1] \"a 10 by 2 mrts matrix\"")
})

test_coverted_matrix <- convertToPositiveDefinite(matrix(c(1, 2, 3, 4), 2, 2))
true_pd_matrix <- matrix(c(1.415476, 2.5, 2.5, 4.415476), 2, 2)
test_that("Convert a matrix to positive definite", {
  expect_lte(norm(test_coverted_matrix - true_pd_matrix, "F"), tolerance)
})

shifted_array <- shiftArray(array(1:10, 2), c(-1, 0, 0))
test_that("Shift an array", {
  expect_equal(sum(shifted_array - c(2, 1)), 0)
  expect_equal(attributes(shifted_array)$dim, 2)
  expect_error(shiftArray(array(1:10, 2), c(-100, 0, 0)), "shift exceeds array dimensions")
})

set.seed(1234)
n <- 150
s <- 5
grid1 <- grid2 <- seq(0, 1, l = 30)
grids <- expand.grid(grid1, grid2)
Fk <- matrix(0, 900, 2)
Fk[, 1] <- cos(sqrt((grids[, 1] - 0)^2 + (grids[, 2] - 1)^2) * pi)
Fk[, 2] <- cos(sqrt((grids[, 1] - 0.75)^2 + (grids[, 2] - 0.25)^2) * 2 * pi)
w <- matrix(c(rnorm(2, sd = 5), rnorm(2, sd = 3)), 2, 2)
y <- Fk %*% w
obs <- sample(900, n)
epsilon <- rexp(n) * sqrt(s)
data <- y[obs] + epsilon
data_2D <- y[obs, ] + epsilon
M <- matrix(rnorm(4), 2, 2)
M <- (M + t(M)) / 2

estimeated_log_likelihood_K_2 <- computeLikelihood(data, Fk[obs, ], M, s, diag(epsilon))
true_log_likelihood_K_2 <- 935.087343
estimeated_log_likelihood_K_1 <- computeLikelihood(data, as.matrix(Fk[obs, 1]), M[1], s, diag(epsilon))
true_log_likelihood_K_1 <- 1421.554255359496
test_that("Negative log likelihood", {
  expect_lte(estimeated_log_likelihood_K_2 - true_log_likelihood_K_2, tolerance)
  expect_lte(estimeated_log_likelihood_K_1 - true_log_likelihood_K_1, tolerance)
})

selected_basis <- selectBasis(data, grids)
selected_basis_em <- selectBasis(data, grids, method = "EM")
data_2D[sample(1:200, 15)] <- NA
selected_basis_na <- selectBasis(cbind(data_2D, NA), grids[obs, ])

test_that("Basis functions selection", {
  expect_equal(names(attributes(selected_basis)), c("dim", "UZ", "Xu", "nconst", "BBBH", "class"))
  expect_equal(class(selected_basis), "mrts")
  expect_equal(dim(selected_basis), c(900, 112))
  expect_equal(dim(attributes(selected_basis)$UZ), c(153, 112))
  expect_lte(norm(attributes(selected_basis)$UZ, "F") - 2493869, tolerance)
  expect_equal(dim(attributes(selected_basis)$Xu), c(150, 2))
  expect_lte(norm(attributes(selected_basis)$Xu, "F") - 7.206402, tolerance)
  expect_lte(sum(attributes(selected_basis)$nconst - c(0.29846350, 0.04876598)), tolerance)
  expect_equal(dim(attributes(selected_basis)$BBBH), c(3, 150))
  expect_lte(norm(attributes(selected_basis)$BBBH, "F") - 0.09683313, tolerance)
  expect_equal(names(attributes(selected_basis_em)), c("dim", "UZ", "Xu", "nconst", "BBBH", "class"))
  expect_lte(norm(attributes(selected_basis_em)$UZ, "F") - 2493869, tolerance)
  expect_equal(dim(attributes(selected_basis_em)$Xu), c(150, 2))
  expect_lte(norm(attributes(selected_basis_em)$Xu, "F") - 7.206402, tolerance)
  expect_lte(sum(attributes(selected_basis_em)$nconst - c(0.29846350, 0.04876598)), tolerance)
  expect_equal(dim(attributes(selected_basis_em)$BBBH), c(3, 150))
  expect_lte(norm(attributes(selected_basis_em)$BBBH, "F") - 0.09683313, tolerance)
  expect_equal(dim(selected_basis_na), c(150, 13))
  expect_lte(norm(attributes(selected_basis_na)$UZ, "F") - 669565.8, tolerance)
  expect_equal(dim(attributes(selected_basis_na)$Xu), c(149, 2))
  expect_lte(norm(attributes(selected_basis_na)$Xu, "F") - 10.21243, tolerance)
  expect_lte(sum(attributes(selected_basis_na)$nconst - c(0.3217056, 0.2795377)), tolerance)
  expect_equal(dim(attributes(selected_basis_na)$BBBH), c(3, 149))
  expect_lte(norm(attributes(selected_basis_na)$BBBH, "F") - 0.07779111, tolerance)
  expect_warning(selectBasis(data, grids, sequence_rank = 1:10, max_knot = 1000), "The minimum of sequence_rank can not less than 3. Too small values will be ignored.")
  expect_error(selectBasis(data, grids, sequence_rank = 1:2, max_knot = 1000), "Not valid sequence_rank!")
})

test_that("Estimate the parameter v", {
  expect_lte(abs(estimateV(1:3, 2, 100, 3) - 31.3333333), tolerance)
  expect_equal(estimateV(1:30, 2, 30, 30), 0)
})

test_that("Estimate the parameter eta", {
  expect_equal(estimateEta(1:3, 0.1, 1), c(0.0, 0.9, 1.9))
})

test_that("Estimate negative log-likelihood", {
  expect_lte(abs(neg2llik(1:3, 0.1, 1, 30, 30) - 84.32403), tolerance)
  expect_equal(neg2llik(1:3, 0, 0, 30, 30), Inf)
})

test_that("Log determinant", {
  set.seed(1234)
  expect_lte(abs(logDeterminant(matrix(rnorm(4), 2, 2))[1] - 0.9284389), tolerance)
  set.seed(1234)
  expect_lte(abs(logDeterminant(matrix(rnorm(10000), 100, 100))[1] - 176.308193), tolerance)
})

set.seed(1234)
test_matrix <- matrix(rnorm(25), 5, 5)
test_matrix <- (test_matrix + t(test_matrix)) / 2
mle1 <- cMLE(Fk,
  4,
  5,
  diag(1, 5),
  test_matrix,
  wSave = TRUE
)
mle2 <- cMLE(Fk,
  4,
  5,
  diag(1, 5),
  test_matrix,
  wSave = FALSE,
  onlylogLike = FALSE,
)
mle3 <- cMLE(Fk,
  4,
  5,
  diag(1, 5),
  test_matrix,
  wSave = TRUE,
  s = 10
)
test_that("cMLE", {
  expect_lte(abs(mle1$v - 0.003904234), tolerance)
  expect_lte(abs(norm(mle1$M, "F") - 1.373454), tolerance)
  expect_equal(mle1$s, 0)
  expect_lte(abs(mle1$negloglik - -9710.933343), tolerance)
  expect_lte(abs(norm(mle1$L, "F") - 27.275983), tolerance)
  expect_null(mle2$L)
  expect_true(all(mle3$L == 0))
})

likeilihood_result <- computeNegativeLikelihood(100, 10, 1000, 4, diag(1,10), 1)
test_that("helper function for cMLE", {
  expect_error(
    computeNegativeLikelihood(100, 10, 1000, 4, matrix(c(1, 2)), 1),
    "Please input a symmetric matrix"           
  )
  expect_error(
    computeNegativeLikelihood(100, 10, 1000, 4, diag(1,4), 1),
    "Please input the rank of a matrix larger than ncol_Fk = 10"           
  )
  expect_lte(abs(likeilihood_result$negloglik - 3498.2569382), tolerance)
})
