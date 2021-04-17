#'
#' Internal function: Remove attributes of mrts
#'
#' @keywords internal
#' @param x A mrts object
#' @param ... Not used directly
#' @return A matrix object
#'
as.matrix.mrts <- function(x, ...) {
  attr(x, "S") <- NULL
  attr(x, "UZ") <- NULL
  attr(x, "Xu") <- NULL
  attr(x, "nconst") <- NULL
  attr(x, "BBBH") <- NULL
  attr(x, "class") <- NULL
  class(x) <- "matrix"
  return(x)
}

#'
#' Internal function: eigen-decomposition in decreasing order
#'
#' @keywords internal
#' @param mat A matrix
#' @return A list
#'
eigenDecomposeInDecreasingOrder <- function(mat) {
  obj <- eigenDecompose(mat)
  obj$value <- rev(obj$value)
  obj$vector <- obj$vector[, ncol(mat):1]
  return(obj)
}

#'
#' Internal function: calculate the log determinant for the likelihood use.
#'
#' @keywords internal
#' @param R A p x p positive-definite matrix
#' @param L A p x K matrix
#' @param K A numeric
#' @return A numeric
#'
calculateLogDeterminant <- function(R, L, K) {
  first_part_determinant <- spam::determinant(
    diag(1, K) + t(L) %*% solve(R) %*% L,
    logarithm = TRUE
  )$modulus
  second_part_determinant <- spam::determinant(R, logarithm = TRUE)$modulus
  return(first_part_determinant + second_part_determinant)
}

#'
#' Internal function: compute a negative log-likelihood
#'
#' @keywords internal
#' @param data  A \emph{n} by \emph{T} data matrix (NA allowed) with
#' \eqn{z[t]} as the \emph{t}-th column.
#' @param Fk A  \emph{n} by \emph{K} matrix of basis function values with
#'  each column being a basis function taken values at \code{loc}
#' @param M A \emph{K} by \emph{K} symmetric matrix
#' @param s A scalar
#' @param Depsilion A \emph{n} by \emph{n} diagonal matrix
#' @return A numeric
#'
computeLikelihood <- function(Data, Fk, M, s, Depsilon) {
  Data <- as.matrix(Data)
  O <- as.matrix(!is.na(Data))
  TT <- NCOL(Data)
  n2loglik <- sum(O) * log(2 * pi)
  R <- toSparseMatrix(s * Depsilon)
  eg <- eigen(M)
  K <- NCOL(Fk)
  L <- Fk %*% eg$vector %*% diag(sqrt(pmax(eg$value, 0)), K) %*% t(eg$vector)
  for (tt in 1:TT) {
    zt <- Data[O[, tt], tt]
    Rt <- R[O[, tt], O[, tt]]
    Lt <- L[O[, tt], ]
    n2loglik <- n2loglik + calculateLogDeterminant(Rt, Lt, K) + sum(zt * invCz(Rt, Lt, zt))
  }
  return(n2loglik)
}

#'
#' Internal function: check if a numeric-like object is diagonal
#'
#' @keywords internal
#' @param object An R object
#' @return logical
#'
isDiagonal <- function(object) {
  if (is.numeric(object) & (length(object) == 1)) {
    return(TRUE)
  }
  if (is.matrix(object)) {
    return(sum(abs(diag(diag(object)) - object)) < .Machine$double.eps)
  }
  else if (is.spam(object)) {
    x <- diag.spam(diag.of.spam(object), NROW(object))
    return(sum(abs(x - object)) < .Machine$double.eps)
  } else {
    return(FALSE)
  }
}

#'
#' Internal function: set 'nc' for LKrigInfo
#'
#' @keywords internal
#' @param z A matrix
#' @param location A location matrix
#' @param nlelve An integer
#' @return numeric
#'
setNC <- function(z, location, nlevel) {
  location_dim <- NCOL(location)
  n <- nrow(z)
  a <- sum(2^(location_dim * (0:(nlevel - 1))))
  nc_estimate <- round((n / a)^(1 / location_dim))
  return(max(4, nc_estimate))
}

#'
#' Internal function: sampling knots
#'
#' @keywords internal
#' @param x A matrix or an array
#' @param nknot A location matrix
#' @param xrng An array including two integers
#' @return sampling knots
#'
subKnot <- function(x, nknot, xrng = NULL, nsamp = 1) {
  x <- apply(as.matrix(x), 2, sort)
  xdim <- dim(x)

  if (is.null(xrng)) {
    xrng <- apply(x, 2, range)
  }
  # To-do: move out the function based on SRP
  mysamp <- function(z_and_id) {
    z <- as.double(names(z_and_id))
    if (length(z) == 1L) {
      return(z)
    }
    else {
      set.seed(mean(z_and_id))
      return(sample(z, size = min(nsamp, length(z))))
    }
  }

  rng <- sqrt(xrng[2, ] - xrng[1, ])
  rng[rng == 0] <- min(rng[rng > 0]) / 5
  rng <- rng * 10 / min(rng)
  rng_max_index <- which.max(rng)

  nmbin <- round(exp(log(rng) * log(nknot) / sum(log(rng))))
  nmbin <- pmax(2, nmbin)
  while (prod(nmbin) < nknot) {
    nmbin[rng_max_index] <- nmbin[rng_max_index] + 1
  }

  gvec <- matrix(1, nrow = xdim[1])
  cnt <- 0
  while (length(unique(gvec)) < nknot) {
    nmbin <- nmbin + cnt
    kconst <- 1
    gvec <- matrix(1, nrow = xdim[1])
    for (kk in 1:xdim[2]) {
      grp <- pmin(
        round((nmbin[kk] - 1) * ((x[, kk] - xrng[1, kk]) / (xrng[2, kk] - xrng[1, kk]))),
        nmbin[kk] - 1L
      )
      if (length(unique(grp)) < nmbin[kk]) {
        brk <- quantile(x[, kk], seq(0, 1, l = nmbin[kk] + 1))
        brk[1] <- brk[1] - 0.1^8
        grp <- as.double(cut(x[, kk], brk))
      }
      gvec <- gvec + kconst * grp
      kconst <- kconst * nmbin[kk]
    }
    cnt <- cnt + 1
  }
  # To-do: refactor the following four lines
  gvec <- as.factor(gvec)
  gid <- as.double(as.character(gvec))
  names(gid) <- 1:xdim[1]
  index <- unlist(tapply(gid, gvec, mysamp))
  return(x[index, ])
}

#'
#' Internal function: convert to a sparse matrix
#'
#' @keywords internal
#' @param mat A matrix
#' @return An array of indeces
#'
fetchNonZeroIndexs <- function(mat) {
  if (!is.matrix(mat)) {
    stop(paste0(c("Wrong matrix format, but got ", class(mat))))
  }
  db <- tempfile()
  NR <- NROW(mat)
  NC <- NCOL(mat)
  f <- fm.create(db, NR, NC)
  f[, 1:NCOL(mat)] <- mat
  j <- sapply(1:NC, function(j) which(f[, j] != 0))
  ridx <- unlist(j)
  k <- sapply(1:NR, function(k) rbind(k, which(f[k, ] != 0)))
  kk <- matrix(unlist(k), ncol = 2, byrow = T)
  cidx <- sort(kk[, 2])
  where <- (cidx - 1) * NR + ridx
  closeAndDeleteFiles(f)
  return(where)
}

#'
#' Internal function: convert to a sparse matrix
#'
#' @keywords internal
#' @param mat A matrix or a dataframe
#' @param verbose A boolean
#' @return sparse matrix
#'
toSparseMatrix <- function(mat, verbose = FALSE) {
  if (is.spam(mat)) {
    if (verbose) message("The input is already a sparse matrix")
    return(mat)
  }
  if (!(is.data.frame(mat) || is.matrix(mat))) {
    stop(paste0(c("Wrong format for toSparseMatrix, but got ", class(mat))))
  }
  if (is.data.frame(mat)) {
    mat <- as.matrix(mat)
  }

  if (length(mat) > 1e8) {
    warnings("Use sparse matrix as input instead; otherwise it could take a very long time!")
    where <- fetchNonZeroIndexs(mat)
  } else {
    where <- which(mat != 0)
  }
  matrix_dim <- dim(mat)
  sparse_matrix <- spam(0, nrow = NROW(mat), ncol = NCOL(mat))
  for (j in 1:matrix_dim[2]) {
    flat_indexs <- where[where >= ((j - 1) * matrix_dim[1] + 1) & where <= (j * matrix_dim[1])]
    if (length(flat_indexs) > 0) {
      row_indexs <- flat_indexs - ((j - 1) * matrix_dim[1])
      sparse_matrix[row_indexs, j] <- mat[row_indexs, j]
    }
  }
  return(sparse_matrix)
}

#'
#' Internal function: internal matrix calcuation (will be deprecated)
#'
#' @keywords internal
#' @param R A p x p matrix
#' @param L A p x K matrix
#' @param z An array with length p or 1 x p matrix
#' @return A 1 x p matrix
#'
ZinvC <- function(R, L, z) {
  K <- NCOL(L)
  iR <- solve(R)
  ZiR <- z %*% iR
  left <- ZiR %*% L %*% solve(diag(1, K) + t(L) %*% iR %*% L) %*%
    t(L)
  return(ZiR - left %*% iR)
}

#'
#' Internal function: internal matrix calcuation
#'
#' @keywords internal
#' @param R A p x p matrix
#' @param L A p x K matrix
#' @param z An array with length p or 1 x p matrix
#' @return A 1 x p matrix
#'
invCz <- function(R, L, z) {
  K <- NCOL(L)
  iR <- solve(R)
  iRZ <- iR %*% z
  right <- L %*% solve(diag(1, K) + as.matrix(t(L) %*% iR %*% L)) %*% (t(L) %*% iRZ)

  return(iRZ - iR %*% right)
}


#'
#' Internal function: print an FRK object
#'
#' @keywords internal
#' @param x An FRK object
#' @param ... Not used directly
#'
print.FRK <- function(x, ...) {
  if (class(x) != "FRK") {
    stop("Invalid object! Please enter an `FRK` object")
  }
  attr(x, "pinfo") <- NULL
  if (!is.null(x$LKobj)) {
    x$LKobj <- x$LKobj$summary
  }
  out <- paste0("a ", NROW(x$G), " by ", NCOL(x$G), " mrts matrix")
  return(print(out))
}

#'
#' Internal function: print an mrts object
#'
#' @keywords internal
#' @param x An mrts object
#' @param ... Not used directly
#'
print.mrts <- function(x, ...) {
  if (class(x) != "mrts") {
    stop("Invalid object! Please enter an `mrts` object")
  }

  if (NCOL(x) == 1) {
    return(c(x))
  } else {
    return(x[, 1:NCOL(x)])
  }
}

#'
#' Internal function: convert a matrix to positive definite
#'
#' @keywords internal
#' @param mat A matrix or a dataframe
#' @return A positive-definite matrix
#'
convertToPositiveDefinite <- function(mat) {
  v <- try(min(eigen(mat, only.values = T)$value), silent = TRUE)
  if (is(v, "try-error") || !isSymmetric.matrix(mat)) {
    mat <- (mat + t(mat)) / 2
    v <- min(eigen(mat, only.values = T)$value)
  }
  if (v <= 0) {
    mat <- mat + diag(max(0, -v) + 0.1^7.5, NROW(mat))
  }
  return(mat)
}

extractLK <- function(obj, loc = NULL, w = NULL, pick = NULL) {
  out <- list()
  if (is.null(loc)) {
    if (is.null(obj$LKinfo.MLE$x)) {
      loc <- obj$LKinfo.MLE$call["x"][[1]]
    } else {
      loc <- obj$LKinfo.MLE$x
    }
  }
  phi <- LKrig.basis(loc, obj$LKinfo)
  Q <- LKrig.precision(obj$LKinfo)
  out$Q <- Q
  if (!is.null(w)) {
    out$weights <- w
  } else {
    out$weights <- obj$LKinfo.MLE$weights
  }
  w <- diag.spam(sqrt(out$weights))
  wX <- w %*% phi
  out$wX <- wX
  out$G <- t(wX) %*% wX + obj$lambda.MLE * Q
  out$lambda <- obj$lambda.MLE
  if (is.null(pick)) {
    pick <- 1:NROW(loc)
  }
  out$pick <- pick

  return(out)
}

#'
#' Internal function: A wrapper for LatticeKrig::LKrigSetup
#'
#' @keywords internal
#' @param x Spatial locations that define the spatial domain for prediction.
#' @param nlevel Number of levels in multi-resolution.
#' @param alpha A vector of length nlevel with the relative variances for the different multi- resolution levels.
#' @param a.wght The correlation range in the SAR model.
#' @param NC The maximum number of lattice grid points for a spatial coordinate and at the coarsest level of resolution.
#' @param lambda A smoothing parameter.
#' @param LKGeometry A text string that gives the names of the model geometry.
#' @param ... Not used directly
#' @return list
#'
setUpKrigInfo <- function(x = NULL,
                          nlevel = NULL,
                          alpha = NA,
                          a.wght = NA,
                          NC = NULL,
                          lambda = NA,
                          LKGeometry = "LKRectangle",
                          ...) {
  setupArgs <- list(...)
  LKinfo <- list(
    x = x,
    nlevel = nlevel,
    alpha = alpha,
    alphaObject = NULL,
    a.wght = a.wght,
    a.wghtObject = NULL,
    NC = NC,
    NC.buffer = NULL,
    nu = NULL,
    normalize = TRUE,
    lambda = lambda,
    sigma = NA,
    rho = NA,
    rho.object = NULL,
    LKGeometry = LKGeometry,
    distance.type = "Euclidean",
    BasisFunction = "wendland",
    overlap = 2.5,
    V = NULL,
    BasisType = "Radial",
    fixedFunction = "LKrigDefaultFixedFunction",
    fixedFunctionArgs = list(m = 2),
    collapseFixedEffect = FALSE,
    max.points = NULL,
    mean.neighbor = 50,
    choleskyMemory = NULL,
    setupArgs = setupArgs,
    dense = FALSE
  )

  LKinfo$basisInfo <- list(
    BasisType = "Radial",
    BasisFunction = "wendland",
    overlap = 2.5,
    max.points = NULL,
    mean.neighbor = 50,
    V = NULL
  )

  class(LKinfo) <- c("LKinfo", LKGeometry)
  LKinfo <- setDefaultsLKinfo(LKinfo)

  LKinfo$latticeInfo <- do.call(
    "LKrigSetupLattice",
    c(
      list(object = LKinfo, verbose = FALSE),
      setupArgs
    )
  )

  LKinfo$a.wght <- setDefaultAwght(LKinfo)

  LKinfo$call <- NULL
  LKinfoCheck(LKinfo)

  return(LKinfo)
}

#'
#' Internal function: A wrapper of LatticeKrig::LKrigSetupAwght
#'
#' @keywords internal
#' @param LKinfo LKinfo object
#' @return list
#'
setDefaultAwght <- function(LKinfo) {
  a_wght <- LKinfo$a.wght
  nlevel <- LKinfo$nlevel
  isotropic <- length(a_wght) == 1

  if (!is.list(a_wght)) {
    if (nlevel == 1) {
      a_wght <- list(a_wght)
    }
    else {
      if (length(a_wght) == 1) {
        a_wght <- rep(a_wght, nlevel)
      }
      a_wght <- as.list(c(a_wght))
    }
  }
  if (length(a_wght) != nlevel) {
    stop("length of a_wght list differs than of nlevel")
  }
  attr(a_wght, "fastNormalize") <- FALSE
  attr(a_wght, "isotropic") <- isotropic
  return(a_wght)
}

#'
#' Internal function: A modifier of LatticeKrig::LKrig.basis
#'
#' @keywords internal
#' @param x1 A matrix
#' @param LKinfo LKinfo object
#' @return A matrix
#'
calculateLatticeKrigBasis <- function(x1, LKinfo) {
  x1 <- as.matrix(x1)

  PHI <- NULL
  basis_delta <- LKinfo$latticeInfo$delta * LKinfo$basisInfo$overlap
  for (l in 1:LKinfo$nlevel) {
    centers <- LKinfo$latticeInfo$grid[[l]]
    if (LKinfo$LKGeometry == "LKBox") class(centers) <- "gridList"

    PHItemp <- Radial.basis(
      x1,
      centers,
      basis_delta[l],
      max.points = LKinfo$basisInfo$max.points,
      mean.neighbor = LKinfo$basisInfo$mean.neighbor,
      BasisFunction = get(LKinfo$basisInfo$BasisFunction),
      distance.type = LKinfo$distance.type
    )
    # Normalization
    wght <- normalizeBasis(LKinfo, level = l, PHI = PHItemp)
    indZero <- wght == 0
    if (any(indZero)) {
      warning("Some normalization weights are zero")
    }
    wght[indZero] <- 1.0
    if (nrow(x1) > 1) {
      PHItemp <- diag.spam(1 / sqrt(wght)) %*% PHItemp
    }
    else {
      PHItemp@entries <- PHItemp@entries / sqrt(wght)
    }

    if (length(LKinfo$alpha[[l]]) > 1) {
      PHItemp <- diag.spam(sqrt(LKinfo$alpha[[l]])) %*% PHItemp
    }
    else {
      PHItemp <- sqrt(LKinfo$alpha[[l]]) * PHItemp
    }
    PHI <- spam::cbind.spam(PHI, PHItemp)
  }

  return(PHI)
}

#'
#' Internal function: A modifier of LatticeKrig::LKrigNormalizeBasis
#'
#' @keywords internal
#' @param x1 A matrix
#' @param LKinfo LKinfo object
#' @return A matrix
#'
normalizeBasis <- function(LKinfo, level, PHI) {
  if (LKinfo$LKGeometry == "LKInterval") {
    tempB <- calculateSARForOneDimLocation(LKinfo, level)
  } else if (LKinfo$LKGeometry == "LKRectangle") {
    tempB <- calculateSARForTowDimLocation(LKinfo, level)
  } else {
    tempB <- calculateSARForThreeDimLocation(LKinfo, level)
  }

  # tempB is in spind format
  tempB <- spam(tempB[c("ind", "ra")], nrow = tempB$da[1], ncol = tempB$da[2])
  # quadratic form involves applying inverse precision matrix to basis function evaluted at
  # locations for evaluation
  wght <- LKrig.quadraticform(t(tempB) %*% tempB,
    PHI = PHI,
    choleskyMemory = LKinfo$choleskyMemory
  )
  return(wght)
}

#'
#' Internal function: A modifier of LatticeKrig::LKrigSAR.LKInterval
#'
#' @keywords internal
#' @param LKinfo LKinfo object
#' @param level an integer
#' @return A list
#'
calculateSARForOneDimLocation <- function(LKinfo, level) {
  m <- LKinfo$latticeInfo$mLevel[level]
  a.wght <- (LKinfo$a.wght)[[level]]
  if (length(a.wght) == 1) {
    a.wght <- rep(a.wght, m)
  }
  if (length(a.wght) != m) {
    cat("Level, m, length( a.wght): ", level, m, length(a.wght), fill = TRUE)
    stop("a.wght wrong length")
  }
  da <- c(m, m)
  ra <- c(a.wght, rep(-1, (m - 1)), rep(-1, (m - 1)))
  Bi <- c(1:m, 2:m, 1:(m - 1))
  Bj <- c(1:m, 1:(m - 1), 2:m)
  return(list(ind = cbind(Bi, Bj), ra = ra, da = da))
}

#'
#' Internal function: A modifier of LatticeKrig::LKrigSAR.LKBox
#'
#' @keywords internal
#' @param LKinfo LKinfo object
#' @param level an integer
#' @return A list
#'
calculateSARForTowDimLocation <- function(LKinfo, level) {
  m <- LKinfo$latticeInfo$mLevel[level]
  a.wght <- (LKinfo$a.wght)[[level]]
  if (length(a.wght) > 1) {
    stop("a.wght must be constant")
  }
  da <- c(m, m)

  ra <- c(rep(a.wght, m), rep(-1, m * 6))
  Bi <- c(rep(1:m, 7))
  Bindex <- array(1:m, LKinfo$latticeInfo$mx[level, ])
  Bj <- c(
    1:m,
    c(shiftArray(Bindex, c(-1, 0, 0))),
    c(shiftArray(Bindex, c(1, 0, 0))),
    c(shiftArray(Bindex, c(0, -1, 0))),
    c(shiftArray(Bindex, c(0, 1, 0))),
    c(shiftArray(Bindex, c(0, 0, -1))),
    c(shiftArray(Bindex, c(0, 0, 1)))
  )
  in_range <- !is.na(Bj)
  Bi <- Bi[in_range]
  Bj <- Bj[in_range]
  ra <- ra[in_range]
  return(list(ind = cbind(Bi, Bj), ra = ra, da = da))
}

#'
#' Internal function: A modifier of LatticeKrig::LKrigSAR.LKRectangle
#'
#' @keywords internal
#' @param LKinfo LKinfo object
#' @param level an integer
#' @return A list
#'
calculateSARForThreeDimLocation <- function(LKinfo, level) {
  mx1 <- LKinfo$latticeInfo$mx[level, 1]
  mx2 <- LKinfo$latticeInfo$mx[level, 2]
  m <- mx1 * mx2
  a.wght <- (LKinfo$a.wght)[[level]]

  stationary <- (attr(LKinfo$a.wght, "stationary"))[level]
  first.order <- attr(LKinfo$a.wght, "first.order")[level]
  isotropic <- attr(LKinfo$a.wght, "isotropic")[level]
  distance.type <- LKinfo$distance.type
  if (all(stationary & isotropic)) {
    if (any(unlist(a.wght) < 4)) {
      stop("a.wght less than 4")
    }
  }

  dim.a.wght <- dim(a.wght)
  index <- c(5, 4, 6, 2, 8, 3, 9, 1, 7)
  da <- as.integer(c(m, m))
  if (first.order) {
    ra <- array(NA, c(mx1 * mx2, 5))
    ra[, 1] <- a.wght
    ra[, 2:5] <- -1
  }
  else {
    ra <- array(NA, c(mx1 * mx2, 9))
    for (kk in 1:9) {
      if (stationary) {
        ra[, kk] <- a.wght[index[kk]]
      }
      else {
        ra[, kk] <- a.wght[, index[kk]]
      }
    }
  }

  Bi <- rep(1:m, 5)
  i.c <- matrix(1:m, nrow = mx1, ncol = mx2)
  Bj <- c(
    i.c,
    LKrig.shift.matrix(i.c, 0, -1),
    LKrig.shift.matrix(i.c, 0, 1),
    LKrig.shift.matrix(i.c, 1, 0),
    LKrig.shift.matrix(i.c, -1, 0)
  )
  if (!first.order) {
    Bi <- c(Bi, rep(1:m, 4))
    Bj <- c(
      Bj,
      LKrig.shift.matrix(i.c, 1, 1),
      LKrig.shift.matrix(i.c, -1, 1),
      LKrig.shift.matrix(i.c, 1, -1),
      LKrig.shift.matrix(i.c, -1, -1)
    )
  }

  good <- !is.na(Bj)
  Bi <- as.integer(Bi[good])
  Bj <- as.integer(Bj[good])
  ra <- c(ra)[good]

  if (!is.null(LKinfo$setupArgs$BCHook)) {
    M <- da[1]
    for (i in 1:M) {
      rowI <- which(Bi == i)
      rowNN <- rowI[-1]
      ra[rowNN] <- 4 * ra[rowNN] / length(rowNN)
    }
  }
  return(list(ind = cbind(Bi, Bj), ra = ra, da = da))
}

#'
#' Internal function: A modifier of LatticeKrig::LKArrayShift
#'
#' @keywords internal
#' @param array_object array
#' @param shift_index one-dim array
#'
shiftArray <- function(array_object, shift_index) {
  shape <- dim(array_object)
  reshaped_array <- array(NA, shape)
  if (any(abs(shift_index) > shape)) {
    stop("shift exceeds array dimensions")
  }

  index_list_source <- index_list_target <- NULL

  for (k in 1:length(shape)) {
    index_list_source <- c(index_list_source, list(c(1:shape[k])))
    temp_index <- (0:(shape[k] - 1) + shift_index[k]) %% shape[k] + 1
    temp_index[(temp_index < 1) | (temp_index > shape[k])] <- NA
    index_list_target <- c(index_list_target, list(temp_index))
  }
  index_source <- as.matrix(expand.grid(index_list_source))
  index_target <- as.matrix(expand.grid(index_list_target))
  in_range <- rowSums(is.na(index_target)) == 0
  reshaped_array[index_target[in_range, ]] <- array_object[index_source[in_range, ]]
  return(reshaped_array)
}

#'
#' Internal function: A wendland function with k = 2 and l = 4
#'
#' @keywords internal
#' @param r A matrix of 2 or 3 d locations
#' @return A matrix
#'
wendland <- function(r) {
  # Ref: https://arxiv.org/pdf/1203.5696.pdf
  if (any(r < 0)) {
    stop(paste(c("Invalid values: ", r[which(r < 0)]), " "))
  }
  return(r < 1) * (1 - r)^6 * (35 * r^2 + 18 * r + 3)
}

#'
#' Internal function: Convert value to bytes
#'
#' @keywords internal
#' @param input_array An array consisting of the value and the unit
#' @return A scalar in byte
#'
toBytes <- function(input_array) {
  num <- as.numeric(input_array[1])
  # Avoid case-sensitive
  units <- tolower(input_array[2])
  lookup <- list(
    "kb" = "kilobytes",
    "mb" = "megabytes",
    "gb" = "gigabytes",
    "tb" = "terabytes"
  )

  if (units %in% lookup) {
    power <- which(units == lookup)
    return(as.numeric(num * 1024^power))
  }
  else if (units %in% names(lookup)) {
    power <- which(units == names(lookup))
    return(num * 1024^power)
  }
  else {
    return(num)
  }
}

#'
#' Internal function: fetch the system ram
#'
#' @keywords internal
#' @param os An OS name with type
#' @return An integer
#'
fetchSystemRam <- function(os) {
  if (grepl("*linux*", os)) {
    cmd <- "awk '/MemTotal/ {print $2}' /proc/meminfo"
    ram <- system(cmd, intern = TRUE, ignore.stderr = TRUE)
    ram <- as.numeric(ram) * 1024
  }
  else if (grepl("*darwin*", os)) {
    ram <- system("system_profiler -detailLevel mini | grep \"  Memory:\"",
      intern = TRUE,
      ignore.stderr = TRUE
    )[1]
    ram <- removeWhitespace(ram)
    ram <- toBytes(unlist(strsplit(ram, " "))[2:3])
  }
  else if (grepl("*solaris*", os)) {
    cmd <- "prtconf | grep Memory"
    ram <- system(cmd, intern = TRUE, ignore.stderr = TRUE)
    ram <- removeWhitespace(ram)
    ram <- toBytes(unlist(strsplit(ram, " "))[3:4])
  }
  else {
    ram <- system("wmic MemoryChip get Capacity", intern = TRUE)[-1]
    ram <- removeWhitespace(ram)
    ram <- ram[nchar(ram) > 0]
    ram <- sum(as.numeric(ram))
  }
  return(as.double(ram))
}

#'
#' Internal function: remove white spaces for a given string
#'
#' @keywords internal
#' @param x A character
#' @return A character
#'
removeWhitespace <- function(x) gsub("(^[[:space:]]+|[[:space:]]+$)", "", x)
