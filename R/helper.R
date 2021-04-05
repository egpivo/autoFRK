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
#' @param matrix A matrix
#' @return A list
#'
eigenDecomposeInDecreasingOrder <- function(matrix) {
  obj <- eigenDecompose(matrix)
  obj$value <- rev(obj$value)
  obj$vector <- obj$vector[, ncol(matrix):1]
  return(obj)
}

getLikelihood <- function(Data, Fk, M, s, Depsilon) {
  logdet <- function(R, L, K) {
    spam::determinant(diag(1, K) + t(L) %*% solve(R) %*%
      L, logarithm = TRUE)$modulus + spam::determinant(R,
      logarithm = TRUE
    )$modulus
  }
  Data <- as.matrix(Data)
  O <- as.matrix(!is.na(Data))
  TT <- NCOL(Data)
  n2loglik <- sum(O) * log(2 * pi)
  R <- toSpMat(s * Depsilon)
  eg <- eigen(M)
  L <- Fk %*% eg$vector %*% diag(sqrt(pmax(eg$value, 0))) %*%
    t(eg$vector)
  K <- NCOL(Fk)
  for (tt in 1:TT) {
    zt <- Data[O[, tt], tt]
    Rt <- R[O[, tt], O[, tt]]
    Lt <- L[O[, tt], ]
    if (NCOL(Lt) == 1) {
      Lt <- t(Lt)
      n2loglik <- n2loglik + log(Rt + Lt %*% t(Lt))
    }
    else {
      n2loglik <- n2loglik + logdet(Rt, Lt, K) + sum(zt *
        invCz(Rt, Lt, zt))
    }
  }
  return(n2loglik)
}

ifElse <- function(cond, yes_out, no_out) {
  if (cond) {
    return(yes_out)
  } else {
    return(no_out)
  }
}

checkDiag <- function(X) {
  if (is(X, "numeric") & (length(X) == 1)) {
    return(TRUE)
  }
  if (is(X, "matrix")) {
    if (sum(abs(diag(diag(X)) - X)) < .Machine$double.eps) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
  else {
    x <- diag.spam(diag.of.spam(X), NROW(X))
    return(identical(x, X))
  }
}

setNC <- function(z, loc, nlevel) {
  Dimension <- NCOL(loc)
  N <- nrow(z)
  a <- sum(2^(Dimension * (0:(nlevel - 1))))
  NCtest <- (N / a)^(1 / Dimension)
  return(round(max(4, NCtest)))
}

subKnot <- function(x, nknot, xrng = NULL, nsamp = 1) {
  x <- as.matrix(x)
  xdim <- dim(x)
  if (xdim[2] > 1) {
    for (kk in xdim[2]:1) x <- x[order(x[, kk]), ]
  }
  else {
    x <- as.matrix(sort(x))
  }
  if (is.null(xrng)) {
    if (xdim[2] > 1) {
      xrng <- apply(x, 2, range)
    }
    else {
      xrng <- matrix(range(x), 2, 1)
    }
  }
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
  nmbin <- round(exp(log(rng) * log(nknot) / sum(log(rng))))
  nmbin <- pmax(2, nmbin)
  while (prod(nmbin) < nknot) {
    nmbin[which.max(rng)] <- nmbin[which.max(rng)] +
      1
  }
  gvec <- matrix(1, xdim[1], 1)
  cnt <- 0
  while (length(unique(gvec)) < nknot) {
    nmbin <- nmbin + cnt
    gvec <- matrix(1, xdim[1], 1)
    kconst <- 1
    for (kk in 1:xdim[2]) {
      grp <- pmin(
        round((nmbin[kk] - 1) * ((x[, kk] - xrng[1, kk]) / (xrng[2, kk] - xrng[1, kk]))),
        nmbin[kk] -1L
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
  gvec <- as.factor(gvec)
  gid <- as.double(as.character(gvec))
  names(gid) <- 1:xdim[1]
  index <- unlist(tapply(gid, gvec, mysamp))
  if (xdim[2] > 1) {
    x[index, ]
  } else {
    x[index]
  }
}

systemRam <- function(os) {
  remove_white <- function(x) {
    gsub(
      "(^[[:space:]]+|[[:space:]]+$)",
      "", x
    )
  }
  toBytes <- function(value) {
    num <- as.numeric(value[1])
    units <- value[2]
    power <- match(units, c("kB", "MB", "GB", "TB"))
    if (!is.na(power)) {
      return(num * 1024^power)
    }
    power <- match(units, c(
      "Kilobytes", "Megabytes", "Gigabytes",
      "Terabytes"
    ))
    if (!is.na(power)) {
      return(num * 1024^power)
    }
    num
  }
  if (length(grep("^linux", os))) {
    cmd <- "awk '/MemTotal/ {print $2}' /proc/meminfo"
    ram <- system(cmd, intern = TRUE)
    ram <- as.numeric(ram) * 1024
  }
  else if (length(grep("^darwin", os))) {
    ram <- system("system_profiler -detailLevel mini | grep \"  Memory:\"",
      intern = TRUE
    )[1]
    ram <- remove_white(ram)
    ram <- toBytes(unlist(strsplit(ram, " "))[2:3])
  }
  else if (length(grep("^solaris", os))) {
    cmd <- "prtconf | grep Memory"
    ram <- system(cmd, intern = TRUE)
    ram <- remove_white(ram)
    ram <- toBytes(unlist(strsplit(ram, " "))[3:4])
  }
  else {
    ram <- system("wmic MemoryChip get Capacity", intern = TRUE)[-1]
    ram <- remove_white(ram)
    ram <- ram[nchar(ram) > 0]
    sum(as.numeric(ram))
  }
  as.double(ram)
}

toSpMat <- function(mat) {
  if (is(mat, "data.frame")) {
    mat <- as.matrix(mat)
  }
  if (is(mat, "matrix")) {
    if (length(mat) > 10^8) {
      warnings("Use sparse matrix as input instead; otherwise it could take a very long time!")
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
    }
    else {
      where <- which(mat != 0)
      ridx <- row(mat)[where]
      cidx <- col(mat)[where]
    }
    nonzero <- mat[where]
    mat <- spam(0, nrow = NROW(mat), ncol = NCOL(mat))
    mat[ridx, cidx] <- nonzero
  }
  return(mat)
}

ZinvC <- function(R, L, z) {
  K <- NCOL(L)
  iR <- solve(R)
  ZiR <- z %*% iR
  left <- ZiR %*% L %*% solve(diag(1, K) + t(L) %*% iR %*% L) %*%
    t(L)
  ZiR - left %*% iR
}

print.FRK <- function(x, ...) {
  attr(x, "pinfo") <- NULL
  if (!is.null(x$LKobj)) {
    x$LKobj <- x$LKobj$summary
  }
  out <- paste("a ", NROW(x$G), " by ", NCOL(x$G), " mrts matrix",
    sep = ""
  )
  print(out)
}

print.mrts <- function(x, ...) {
  if (NCOL(x) == 1) {
    out <- c(x)
  } else {
    out <- x[, 1:NCOL(x)]
  }
  print(out)
}

mkpd <- function(M) {
  v <- try(min(eigen(M, only.values = T)$value), silent = TRUE)
  if (is(v, "try-error")) {
    M <- (M + t(M)) / 2
    v <- min(eigen(M, only.values = T)$value)
  }
  if (v <= 0) {
    M <- M + diag(max(0, -v) + 0.1^7.5, NROW(M))
  }
  return(M)
}


extractLK <- function(obj, loc = NULL, w = NULL, pick = NULL) {
  out <- list()
  if (is.null(loc)) {
    loc <- ifElse(
      is.null(obj$LKinfo.MLE$x), obj$LKinfo.MLE$call["x"][[1]],
      obj$LKinfo.MLE$x
    )
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
#' @return list
#'
LKrigSetupWrapper <- function(x = NULL,
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
