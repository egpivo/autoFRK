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

autoFRK <- function(Data, loc, mu = 0, D = diag.spam(NROW(Data)), G = NULL, isFineScale = FALSE,
                    maxit = 50, tolerance = 1e-6, maxK = NULL, Kseq = NULL, method = c("fast", "EM"),
                    n.neighbor = 3, maxknot = 5000) {
    # Validate `method`
    if (!missing(method) & length(method) != 1) 
        stop("Please enter one method at a time")
    method <- match.arg(method)
    # Coerce type of `Data`` as matrix; center `Data`
    if (!is.matrix(Data)) Data <- as.matrix(Data)
    Data <- Data - mu
    # Coerce type of `loc` as matrix
    if (!is.matrix(loc)) loc <- as.matrix(loc)
    # Determine basis functions (n x K matrix)
    if (!is.null(G)) 
        Fk <- G
    else {
        Fk <- selectBasis(Data = Data,
                          loc = loc,
                          D = D,
                          maxit = maxit,
                          avgtol = tolerance,
                          maxK = maxK,
                          Kseq = Kseq,
                          method = method,
                          n.neighbor = n.neighbor,
                          maxknot = 5000
        )
    }
    
    K <- NCOL(Fk)
    # Execute KNN for imputation
    if (method == "fast") {
        # Fill missing data by k-nearest-neighbor imputation
        Data <- apply(Data, 2, imputeByKnn , loc=loc, k=n.neighbor)
    }
    if (!isFineScale) 
        obj <- indeMLE(Data = Data,
                       Fk = Fk[, 1:K],
                       D = D,
                       maxit = maxit,
                       avgtol = tolerance, 
                       wSave = TRUE)
    else {
        nu <- .Options$LKinfoSetup$nu
        nlevel <- .Options$LKinfoSetup$nlevel
        a.wght <- .Options$LKinfoSetup$a.wght
        NC <- .Options$LKinfoSetup$NC
        if (is.null(nu)) 
            nu <- 1
        if (is.null(nlevel)) 
            nlevel <- 3
        iniobj <- LKnFRKini(Data, loc, nlevel, weights = 1/diag(D), 
                            n.neighbor, nu)
        DnLK <- LKnFRKopt(iniobj, Fk[, 1:K], nc = NC, a.wght = a.wght)
        DfromLK <- DnLK$DfromLK
        LKobj <- DnLK$LKobj
        Depsilon <- diag.spam(iniobj$weight[iniobj$pick], 
                              length(iniobj$weight[iniobj$pick]))
        obj <- indeMLE(Data = Data,
                       Fk = Fk[, 1:K],
                       D = D,
                       maxit = maxit,
                       avgtol = tolerance, 
                       wSave = TRUE,
                       DfromLK = DfromLK,
                       vfixed = DnLK$s)
    }
    
    obj$G <- Fk
    if (isFineScale) {
        obj$LKobj <- LKobj
        attr(obj, "pinfo")$loc <- loc
        attr(obj, "pinfo")$weights <- 1/diag(D)
    }
    else
        obj$LKobj <- NULL
    class(obj) <- "FRK"
    
    return(obj)
}

checkDiag <- function(X) {
    if (class(X) == "numeric")
        status <- TRUE
    else if (class(X) == "matrix") 
        status <- ifelse(sum(abs(diag(diag(X)) - X)) < .Machine$double.eps, TRUE, FALSE) 
    else
        status <- identical(diag.spam(diag.of.spam(X), NROW(X)), X)
    
    return(status)
}

cMLE <- function(Fk, TT, trS, half, JSJ = NULL, s = 0, ldet = NULL, wSave = FALSE,
                 onlylogLike = !wSave, vfixed = NULL) {
    if (is.null(ldet)) ldet <- 0
    n <- nrow(Fk)
    k <- ncol(Fk)
    eg <- eigen(JSJ)
    d <- eg$value[1:k]
    P <- eg$vector[, 1:k]
    if (is.null(vfixed)) v <- sol.v(d, s, trS, n) else v <- vfixed
    dii <- pmax(d, 0)

    if (onlylogLike)
        result_list <- list(negloglik = neg2llik(dii, s, v, trS, n) * TT + ldet * TT)
    else{
        dhat <- sol.eta(dii, s, v)
        M <- half %*% P %*% (dhat * t(P)) %*% half
        dimnames(M) <- NULL
        if (!wSave) 
            L <- NULL
        else {
            L <- Fk %*% t((sqrt(dhat) * t(P)) %*% half)
            if (all(dhat == 0)) dhat[1] <- 1e-10
            L <- as.matrix(L[, dhat > 0])
        }
    
        result_list <- list(v = v,
                            M = M,
                            s = s,
                            negloglik = neg2llik(dii, s, v, trS, n) * TT + ldet * TT,
                            L = L)
    }
    
    return(result_list)
}

cMLEimat <- function(Fk, Data, s, wSave = FALSE, S = NULL, onlylogLike = !wSave) {
    n <- nrow(Fk)
    k <- ncol(Fk)
    TT <- NCOL(Data)
    trS <- sum(rowSums(as.matrix(Data)^2))/TT
    half <- getHalf(Fk, Fk)
    ihF <- half %*% t(Fk)
    if (is.null(S))
        JSJ <- tcrossprod(ihF %*% Data)/TT
    else
        JSJ <- (ihF %*% S) %*% t(ihF)
    JSJ <- (JSJ + t(JSJ))/2
    
    eg <- eigen(JSJ)
    d <- eg$value[1:k]
    v <- sol.v(d, s, trS, n)
    dii <- pmax(d, 0)

    if (onlylogLike)
        result_list <- list(negloglik = neg2llik(dii, s, v, trS, n) * TT)
    else{
        P <- eg$vector[, 1:k]
        dhat <- sol.eta(dii, s, v)
        M <- half %*% P %*% (dhat * t(P)) %*% half
        dimnames(M) <- NULL
        if (!wSave) 
            result_list <- list(v = v,
                                M = M,
                                s = s,
                                negloglik = neg2llik(dii, s, v, trS, n) * TT)
        else {
            L <- Fk %*% t((sqrt(dhat) * t(P)) %*% half)
            if (all(dhat == 0)) dhat[1] <- 1e-10
            L <- L[, dhat > 0]
            invD <- rep(1, n)/(s + v)
            iDZ <- invD * Data
            right <- L %*% (solve(diag(1, NCOL(L)) + t(L) %*% (invD * L)) %*% (t(L) %*% iDZ))
            INVtZ <- iDZ - invD * right
            etatt <- as.matrix(M %*% t(Fk) %*% INVtZ)
            GM <- Fk %*% M
            V <- as.matrix(M - t(GM) %*% invCz((s + v) * diag.spam(n), L, GM))
            
            result_list <- list(v = v,
                                M = M, 
                                s = s, 
                                negloglik = neg2llik(dii, s, v, trS, n) * TT, 
                                w = etatt, 
                                V = V)
        }
    }
    
    return(result_list)
}

cMLElk <- function(Fk, Data, Depsilon, wSave = FALSE, DfromLK, vfixed = NULL) {
    TT <- NCOL(Data)
    N <- NROW(Data)
    lambda <- DfromLK$lambda
    pick <- DfromLK$pick
    wX <- DfromLK$wX[pick, ]
    G <- t(wX) %*% wX + lambda * DfromLK$Q
    weight <- DfromLK$weights[pick]
    wwX <- diag.spam(sqrt(weight)) %*% wX
    wXiG <- (wwX) %*% solve(G)
    iDFk <- weight * Fk - wXiG %*% (t(wwX) %*% as.matrix(Fk))
    iDZ <- weight * Data - wXiG %*% (t(wwX) %*% as.matrix(Data))
    half <- getHalf(Fk, iDFk)
    ihFiD <- half %*% t(iDFk)
    JSJ <- tcrossprod(ihFiD %*% Data)/TT
    JSJ <- (JSJ + t(JSJ))/2

    ldetD <- -nrow(DfromLK$Q) * log(lambda) + ldet(G) - ldet(DfromLK$Q) - sum(log(weight))
    trS <- sum(rowSums(as.matrix(iDZ) * Data))/TT
    out <- cMLE(Fk = Fk,
                TT = TT,
                trS = trS,
                half = half,
                JSJ = JSJ,
                ldet = as.vector(ldetD), 
                wSave = TRUE,
                vfixed = vfixed)

    out$s <- out$v
    out <- out[-which(names(out) == "v")]
    result <- out[-which(names(out) == "L")]

    if (wSave) {
        iDL <- weight * L - wXiG %*% (t(wwX) %*% L)
        itmp <- solve(diag(1, NCOL(L)) + t(L) %*% iDL/result$s)
        iiLiD <- itmp %*% t(iDL/result$s)
        MFiS11 <- out$M %*% t(iDFk)/result$s - ((result$M %*% t(iDFk/result$s)) %*% L) %*% iiLiD
        result$w <- MFiS11 %*% Data
        result$V <- MFiS11 %*% (Fk %*% result$M)
        wlk <- t(wXiG) %*% Data - t(wXiG) %*% L %*% (iiLiD %*% Data)
        ihL <- chol(itmp) %*% t(L)
        attr(result, "pinfo") <- list(wlk = wlk, pick = pick)
    }
    
    return(result)
}

cMLEsp <- function(Fk, Data, Depsilon, wSave = FALSE) {
    De <- toSpMat(Depsilon)
    iD <- solve(De)
    ldetD <- spam::determinant(De, logarithm = TRUE)$modulus
    iDFk <- iD %*% Fk
    half <- getHalf(Fk, iDFk)
    ihFiD <- half %*% t(iDFk)
    TT <- NCOL(Data)
    JSJ <- tcrossprod(ihFiD %*% Data)/TT
    JSJ <- (JSJ + t(JSJ))/2
    trS <- sum(rowSums(as.matrix(iD %*% Data) * Data))/TT
    out <- cMLE(Fk, TT, trS, half, JSJ, s = 0, ldet = ldetD, wSave)
    
    if (wSave) {
        L <- as.matrix(out$L)
        invD <- iD/(out$s + out$v)
        iDZ <- invD %*% Data
        right0 <- L %*% solve(diag(1, NCOL(L)) + t(L) %*% (invD %*% L))
        INVtZ <- iDZ - invD %*% right0 %*% (t(L) %*% iDZ)
        etatt <- out$M %*% t(Fk) %*% INVtZ
        out$w <- as.matrix(etatt)
        GM <- Fk %*% out$M
        iDGM <- invD %*% GM
        out$V <- as.matrix(out$M - t(GM) %*% (iDGM - invD %*%  right0 %*% (t(L) %*% iDGM)))
    }
    out$s <- out$v
    out <- out[-which(names(out) == "v")]
    out <- out[-which(names(out) == "L")]
    
    return(out)
}

EM0miss <- function(Fk, Data, Depsilon, maxit, avgtol, wSave = FALSE, external = FALSE, 
                    DfromLK = NULL, num.report = TRUE, vfixed = NULL) {
    
    saveOLD <- function(external) if (external) save(old, Ptt1, file = oldfile)
    
    O <- !is.na(Data)
    TT <- NCOL(Data)
    K <- NCOL(Fk)
    tmpdir <- tempfile()
    dir.create(tmpdir)
    ftmp <- paste(tmpdir, 1:TT, sep = "/")
    oldfile <- paste(tmpdir, "old_par.Rdata", sep = "/")
    ziDz <- rep(NA, TT)
    ziDB <- matrix(NA, TT, K)
    db <- list()
    D <- toSpMat(Depsilon)
    iD <- solve(D)
    diagD <- checkDiag(D)
    
    if (!is.null(DfromLK)) {
        pick <- DfromLK$pick
        if (is.null(pick)) pick <- 1:length(DfromLK$weights)
        weight <- DfromLK$weights[pick]
        DfromLK$wX <- DfromLK$wX[pick, ]
        wwX <- diag.spam(sqrt(weight)) %*% DfromLK$wX
        lQ <- DfromLK$lambda * DfromLK$Q
    }
    
    for (tt in 1:TT) {
        if (!is.null(DfromLK)) {
            iDt = NULL
            if (sum(O[, tt]) == NROW(O)) 
                wXiG <- wwX %*% solve(DfromLK$G)
            else {
                G <- t(DfromLK$wX[O[, tt], ]) %*% DfromLK$wX[O[, tt], ] + lQ
                wXiG <- wwX[O[, tt], ] %*% solve(G)
            }
            Bt <- as.matrix(Fk[O[, tt], ])
            if (NCOL(Bt) == 1) Bt <- t(Bt)
            iDBt <- as.matrix(weight[O[, tt]] * Bt - wXiG %*% (t(wwX[O[, tt], ]) %*% Bt))
            zt <- Data[O[, tt], tt]
            ziDz[tt] <- sum(zt * as.vector(weight[O[, tt]] * zt - wXiG %*% (t(wwX[O[, tt], ]) %*% zt)))
            ziDB[tt, ] <- t(zt) %*% iDBt
            BiDBt <- t(Bt) %*% iDBt
        }
        else {
            if (!diagD) 
                iDt <- solve(D[O[, tt], O[, tt]]) 
            else 
                iDt <- iD[O[, tt], O[, tt]]
            Bt <- Fk[O[, tt], ]
            if (NCOL(Bt) == 1) Bt <- t(Bt)
            iDBt <- as.matrix(iDt %*% Bt)
            zt <- Data[O[, tt], tt]
            ziDz[tt] <- sum(zt * as.vector(iDt %*% zt))
            ziDB[tt, ] <- t(zt) %*% iDBt
            BiDBt <- t(Bt) %*% iDBt
        }
        if (external)
            db[[tt]] <- dumpObjects(iDBt, zt, BiDBt, external, oldfile, dbName = ftmp[tt])
        else
            db[[tt]] <- list(iDBt = iDBt,
                             zt = zt,
                             BiDBt = BiDBt, 
                             external = external,
                             oldfile = oldfile)
    }
    # gc
    rm(iDt, Bt, iDBt, zt, BiDBt)
    gc()
    
    dif <- Inf
    cnt <- 0
    Z0 <- Data
    Z0[is.na(Z0)] <- 0
    old <- cMLEimat(Fk, Z0, s = 0, wSave = T)
    if(is.null(vfixed)) old$s <- old$v else old$s <- vfixed
    old$M <- mkpd(old$M)
    Ptt1 <- old$M
    saveOLD(external)
    inv <- MASS::ginv
    
    while ((dif > (avgtol * (100 * K^2))) && (cnt < maxit)) {
        etatt <- matrix(0, K, TT)
        sumPtt <- 0
        s1 <- rep(0, TT)
        if (external) load(oldfile)
        for (tt in 1:TT) {
            s1.eta.P <- with(db[[tt]], {
                             iP <- mkpd(MASS::ginv(mkpd(Ptt1)) + BiDBt/old$s)
                             Ptt <- solve(iP)
                             Gt <- as.matrix(Ptt %*% t(iDBt)/old$s)
                             eta <- c(0 + Gt %*% zt)
                             s1kk <- diag(BiDBt %*% (eta %*% t(eta) + Ptt))
                             rbind(s1kk, eta, Ptt)
                        })
            sumPtt <- sumPtt + s1.eta.P[-c(1:2), ]
            etatt[, tt] <- s1.eta.P[2, ]
            s1[tt] <- sum(s1.eta.P[1, ])
        }
        if (is.null(vfixed))
            new <- list(M = (etatt %*% t(etatt) + sumPtt)/TT, 
                        s = max((sum(ziDz) - 2 * sum(ziDB * t(etatt)) + sum(s1))/sum(O), 1e-8))
        else
            new <- list(M = (etatt %*% t(etatt) + sumPtt)/TT, 
                        s = vfixed)
        new$M <- (new$M + t(new$M)) / 2
        dif <- sum(abs(new$M - old$M)) + abs(new$s - old$s)
        cnt <- cnt + 1
        old <- new
        Ptt1 <- old$M
        saveOLD(external)
    }
    
    if (num.report) cat("Number of iteration: ", cnt, "\n")
    unlink(tmpdir, recursive = TRUE)
    n2loglik <- getlik(Data, Fk, new$M, new$s, Depsilon)
    
    if (!wSave) 
        out <- list(M = new$M, s = new$s, negloglik = n2loglik)
    else {
        if (!is.null(DfromLK)) {
            out <- list(M = new$M, s = new$s, negloglik = n2loglik, 
                       w = etatt, V = new$M - etatt %*% t(etatt)/TT)
            dec <- eigen(new$M)
            L <- Fk %*% dec$vector %*% diag(sqrt(pmax(dec$value, 0)))
            weight <- DfromLK$weights[pick]
            wlk <- matrix(NA, NROW(lQ), TT)
            
            for (tt in 1:TT) {
                if (sum(O[, tt]) == NROW(O)) 
                    wXiG <- wwX %*% solve(DfromLK$G)
                else {
                    G <- t(DfromLK$wX[O[, tt], ]) %*% DfromLK$wX[O[, tt], ] + lQ
                    wXiG <- wwX[O[, tt], ] %*% solve(G)
                }
                dat <- Data[O[, tt], tt]
                Lt <- L[O[, tt], ]
                iDL <- weight[O[, tt]] * Lt - wXiG %*% (t(wwX[O[, tt], ]) %*% Lt)
                itmp <- solve(diag(1, NCOL(L)) + t(Lt) %*% iDL/out$s)
                iiLiD <- itmp %*% t(iDL/out$s)
                wlk[, tt] <- t(wXiG) %*% dat - t(wXiG) %*% Lt %*% (iiLiD %*% dat)
            }
            attr(out, "pinfo") <- list(wlk = wlk, pick = pick)
            attr(out, "missing") <- list(miss = toSpMat(1 - O), 
                                         maxit = maxit, 
                                         avgtol = avgtol)
        }
        else {
            out <- list(M = as.matrix(new$M),
                        s = new$s,
                        negloglik = n2loglik, 
                        w = etatt,
                        V = new$M - etatt %*% t(etatt) / TT)
            
            attr(out, "missing") <- list(miss = toSpMat(1 - O), 
                                         maxit = maxit, 
                                         avgtol = avgtol)
        }
    }

    return(out)
}

getHalf <- function(Fk, iDFk) {
    dec <- mgcv::slanczos(t(Fk) %*% iDFk, k = NCOL(Fk))
    sroot <- sqrt(pmax(dec$value, 0))
    sroot[sroot == 0] <- Inf
    sroot <- 1/sroot

    return(dec$vector %*% (sroot * t(dec$vector)))
}

getlik <- function(Data, Fk, M, s, Depsilon) {
    
    logdet <- function(R, L, K) {
        det1 <- spam::determinant(diag(1, K) + t(L) %*% solve(R) %*% L,
                                  logarithm = TRUE)$modulus 
        det2 <- spam::determinant(R, logarithm = TRUE)$modulus
        
        return(det1 + det2)
    }
    
    Data <- as.matrix(Data)
    O <- as.matrix(!is.na(Data))
    TT <- NCOL(Data)
    n2loglik <- sum(O) * log(2 * pi)
    R <- toSpMat(s * Depsilon)
    eg <- eigen(M)
    L <- Fk %*% eg$vector %*% diag(sqrt(pmax(eg$value, 0))) %*% t(eg$vector)
    K <- NCOL(Fk)
    
    for (tt in 1:TT) {
        zt <- Data[O[, tt], tt]
        Rt <- R[O[, tt], O[, tt]]
        Lt <- L[O[, tt], ]
        if (NCOL(Lt) == 1) {
            Lt <- t(Lt)
            n2loglik <- n2loglik + log(Rt + Lt %*% t(Lt))
        }
        else 
            n2loglik <- n2loglik + logdet(Rt, Lt, K) + sum(zt * invCz(Rt, Lt, zt))
    }
    
    return(n2loglik)
}

imputeByKnn <- function(data, loc, k) {
    where <- is.na(data)
    if (sum(where) == 0) 
        next
    cidx <- which(!where)
    nnidx <- FNN::get.knnx(data = loc[cidx, ],
                           query = loc[where, ],
                           k = k)
    nnidx <- array(cidx[nnidx$nn.index], dim(nnidx$nn.index))
    nnval <- array(data[nnidx], dim(nnidx))
    data[where] <- rowMeans(nnval)
    
    return(data)
}

indeMLE <- function(Data, Fk, D = diag.spam(NROW(Data)), maxit = 50, avgtol = 1e-6, 
                    wSave = FALSE, DfromLK = NULL, vfixed = NULL, num.report = TRUE) {

    isWithNA <- sum(is.na(Data)) > 0
    if (class(Data) == "numeric") Data <- as.matrix(Data)
    TT <- NCOL(Data)
    empty <- apply(!is.na(Data), 2, sum) == 0
    notempty <- which(!empty)
    if (sum(empty) > 0) Data <- as.matrix(Data[, notempty])
    if (class(Data) == "numeric") Data <- as.matrix(Data)
    
    del <- which(rowSums(as.matrix(!is.na(Data))) == 0)
    pick <- 1:NROW(Data)
    if (!checkDiag(D)) D0 <- toSpMat(D) else D0 <- diag.spam(diag(D), NROW(Data))

    if (isWithNA && (length(del) > 0)) {
        pick <- pick[-del]
        Data <- Data[-del, ]
        Fk <- Fk[-del, ]
        if (!checkDiag(D))
            D <- D[-del, -del]
        else
            D <- diag.spam(diag(D)[-del], NROW(Data))
        isWithNA <- sum(is.na(Data)) > 0
    }
    N <- NROW(Data)
    K <- NCOL(Fk)
    Depsilon <- toSpMat(D)
    isimat <- checkDiag(D) * (sum(abs(rep(mean(diag(D)), N) -  diag(Depsilon))) < .Machine$double.eps)
    
    if (!isWithNA) {
        if (isimat & is.null(DfromLK)) {
            
            sigma <- ifelse(!is.null(.Options$sigma_FRK),
                            .Options$sigma_FRK,
                            0)
            if (NCOL(Data) == 1) Data <- as.matrix(Data)
            out <- cMLEimat(Fk, Data, s = sigma, wSave)
            
            if (!is.null(out$v)) {
                out$s <- ifelse(sigma == 0, out$v, sigma)
                out <- out[-which(names(out) == "v")]
            }
            if (wSave) {
                w <- matrix(0, K, TT)
                w[, notempty] <- out$w
                out$w <- w
                attr(out, "pinfo") <- list(D = D0, pick = pick)
            }
        }
        else {
            if (is.null(DfromLK)) {
                out <- cMLEsp(Fk, Data, Depsilon, wSave)
                if (wSave) {
                    w <- matrix(0, K, TT)
                    w[, notempty] <- out$w
                    out$w <- w
                    attr(out, "pinfo") <- list(D = D0, pick = pick)
                }
            }
            else {
                out <- cMLElk(Fk, Data, Depsilon, wSave, DfromLK, vfixed)
                if (wSave) {
                    w <- matrix(0, K, TT)
                    w[, notempty] <- out$w
                    out$w <- w
                }
            }
        }
    }
    else {
        out <- EM0miss(Fk,
                       Data,
                       Depsilon,
                       maxit,
                       avgtol,
                       wSave, 
                       external = FALSE,
                       DfromLK,
                       num.report,
                       vfixed)
        if (wSave) {
            w <- matrix(0, K, TT)
            w[, notempty] <- out$w
            out$w <- w
            if (is.null(DfromLK)) 
                attr(out, "pinfo") <- list(D = D0, pick = pick)
        }
    }
    
    return(out)
}

invCz <- function(R, L, z) {
    K <- NCOL(L)
    iR <- solve(R)
    iRZ <- iR %*% z
    right <- L %*% solve(diag(1, K) + as.matrix(t(L) %*% iR %*% L)) %*% (t(L) %*% iRZ)

    return(iRZ - iR %*% right)
}

LKextract <- function(obj, loc = NULL, w = NULL, pick = NULL) {
    out <- list()
    if (is.null(loc))
        if (is.null(obj$LKinfo.MLE$x))
            loc <- obj$LKinfo.MLE$call["x"][[1]]
        else
            loc <- obj$LKinfo.MLE$x
    
    phi <- LKrig.basis(loc, obj$LKinfo)
    Q <- LKrig.precision(obj$LKinfo)
    out$Q <- Q
    
    if (!is.null(w)) out$weights <- w else out$weights <- obj$LKinfo.MLE$weights
    
    w <- diag.spam(sqrt(out$weights))
    wX <- w %*% phi
    out$wX <- wX
    out$G <- t(wX) %*% wX + obj$lambda.MLE * Q
    out$lambda <- obj$lambda.MLE
    
    if (is.null(pick)) pick <- 1:NROW(loc)
    out$pick <- pick
    
    return(out)
}

LKnFRKini <- function(Data, loc, nlevel = 3, weights = NULL, n.neighbor = 3, nu = 1) {
    if (class(Data) == "numeric") Data <- as.matrix(Data)
    empty <- apply(!is.na(Data), 2, sum) == 0
    if (sum(empty) > 0) Data <- Data[, which(!empty)]
    if (class(Data) == "numeric") Data <- as.matrix(Data)
    
    loc <- as.matrix(loc)
    N <- NROW(Data)
    d <- NCOL(loc)
    nas <- sum(is.na(Data))
    del <- which(rowSums(as.matrix(!is.na(Data))) == 0)
    pick <- 1:N
    
    if (length(del) > 0) {
        Data <- Data[-del, ]
        loc <- loc[-del, ]
        pick <- (1:N)[-del]
    }
    
    nas <- sum(is.na(Data))
    
    if (nas > 0) {
        for (tt in 1:NCOL(Data)) {
            where <- is.na(Data[, tt])
            if (sum(where) == 0) next
            cidx <- which(!where)
            nnidx <- FNN::get.knnx(loc[cidx, ],
                                   as.matrix(loc[where, ]),
                                   k = n.neighbor)
            nnidx <- array(cidx[nnidx$nn.index], dim(nnidx$nn.index))
            nnval <- array((Data[, tt])[nnidx], dim(nnidx))
            Data[where, tt] <- rowMeans(nnval)
        }
    }
    
    x <- as.matrix(loc[pick, ])
    z <- as.matrix(Data)
    d <- NCOL(x)
    
    gtype <- ifelse(d == 1, "LKInterval", ifelse(d == 2, "LKRectangle", "LKBox"))

    thetaL <- 2^(-1 * (1:nlevel))
    alpha <- thetaL^(2 * nu)
    alpha <- alpha/sum(alpha)
    n <- NROW(x)
    
    if (is.null(weights)) weights <- rep(1, NROW(z))
    
    return(list(x = x,
                z = z,
                n = n,
                alpha = alpha,
                gtype = gtype,
                weights = weights, 
                nlevel = nlevel,
                loc = loc,
                pick = pick))
}

LKnFRKopt <- function(iniobj, Fk, nc = NULL, Ks = NCOL(Fk), a.wght = NULL) {
    x <- iniobj$x
    z <- iniobj$z
    alpha <- iniobj$alpha
    alpha <- alpha/sum(alpha)
    gtype <- iniobj$gtype
    weights <- iniobj$weights[iniobj$pick]
    nlevel <- iniobj$nlevel
    TT <- NCOL(z)
    Fk <- Fk[iniobj$pick, ]
    
    if (is.null(nc)) nc <- setNC(z, x, nlevel)
    if (is.null(a.wght)) a.wght <- 2 * NCOL(x) + 0.01
    
    info <- LKrigSetup(x = x,
                       a.wght = a.wght,
                       nlevel = nlevel,
                       NC = nc, 
                       alpha = alpha,
                       LKGeometry = gtype,
                       lambda = 1)
    
    loc <- x
    phi <- LKrig.basis(loc, info)
    w <- diag.spam(sqrt(weights))
    wX <- w %*% phi
    wwX <- w %*% wX
    XwX <- t(wX) %*% wX

    iniLike <- function(par, Data = z, full = FALSE) {
        lambda <- exp(par)
        G <- XwX + lambda * Qini
        wXiG <- (wwX) %*% solve(G)
        iDFk <- weights * Fk - wXiG %*% (t(wwX) %*% as.matrix(Fk))
        iDZ <- weights * Data - wXiG %*% (t(wwX) %*% as.matrix(Data))
        ldetD <- -nrow(Qini) * log(lambda) + ldet(G)
        ldetD <- as.vector(ldetD)
        trS <- sum(rowSums(as.matrix(iDZ) * Data))/TT
        half <- getHalf(Fk, iDFk)
        ihFiD <- half %*% t(iDFk)
        LSL <- tcrossprod(ihFiD %*% Data)/TT
        
        if (!full) 
            cMLE(Fk, TT, trS, half, LSL, s = 0, ldet = ldetD, 
                 wSave = FALSE)$negloglik
        else {
            llike <- ldetD - ldet(Qini) - sum(log(weights))
            cMLE(Fk, TT, trS, half, LSL, s = 0, ldet = llike, 
                 wSave = TRUE, onlylogLike = FALSE, vfixed = NULL)
        }
    }
    
    Qini <- LKrig.precision(info)
    sol <- optimize(iniLike, c(-16, 16), tol = .Machine$double.eps^0.025)
    lambda.MLE <- sol$minimum
    out <- iniLike(sol$minimum, z, full = TRUE)
    llike <- out$negloglik
    info.MLE <- LKrigSetup(x = x,
                           a.wght = a.wght,
                           nlevel = nlevel, 
                           NC = nc,
                           alpha = alpha,
                           LKGeometry = gtype,
                           lambda = lambda.MLE)
    info.MLE$llike <- llike
    info.MLE$time <- NA
    Q <- LKrig.precision(info.MLE)
    G <- t(wX) %*% wX + info.MLE$lambda * Q
    
    return(list(DfromLK = list(Q = Q,
                               weights = weights,
                               wX = wX, 
                               G = G,
                               lambda = info.MLE$lambda,
                               pick = iniobj$pick), 
                s = out$v,
                LKobj = list(summary = NULL,
                             par.grid = NULL, 
                             LKinfo.MLE = info.MLE,
                             lnLike.eval = NULL,
                             lambda.MLE = info.MLE$lambda, 
                             call = NA,
                             taskID = NULL)))
}

mkpd <- function(M) {
    v <- try(min(eigen(M, only.values = T)$value), silent = TRUE)
    
    if (class(v) == "trial-error") {
        M <- (M + t(M))/2
        v <- min(eigen(M, only.values = T)$value)
    }
    if (v <= 0) 
        M <- M + diag(max(0, -v) + 0.1^7.5, NROW(M))
    
    return(M)
}

mrts <- function(knot, k, x = NULL) {
    is64bit = length(grep("64", version["system"])) > 0
    rlimit <- ramSize()
    
    if ((!is64bit) & (max(NROW(x), NROW(knot)) > 20000)) 
        stop("Use 64-bit version of R for such a volume of data!")
    
    #
    # Check: xobs consists of redundant lines
    #
    if (NCOL(knot) == 1)
        xobs <- as.matrix(as.double(as.matrix(knot)))
    else
        xobs <- apply(knot, 2, as.double)
    
    Xu <- unique(cbind(xobs))
    
    if (is.null(x) & length(Xu) != length(xobs)) x <- xobs
    
    colnames(Xu) <- NULL
    n <- n.Xu <- NROW(Xu)
    ndims <- NCOL(Xu)
    if (k < (ndims + 1)) 
        stop("k-1 can not be smaller than the number of dimensions!")
    if ((n^3/rlimit) > 25) {
        bmax <- (rlimit/max(NROW(x), NROW(xobs)))^(1/3)^(1/ndims)
        if (bmax < k) {
            bmax <- k
            warnings("Set a smaller k; or it may eat up all your RAM!")
        }
        Xu <- subknot(Xu, bmax)
        Xu <- as.matrix(Xu)
        n <- NROW(Xu)
        n.Xu <- n
    }
    
    xobs_diag <- diag(sqrt(n/(n - 1))/apply(xobs, 2, sd), ndims)
    
    if (!is.null(x)) {
        if (NCOL(x) == 1)
            x <- as.matrix(as.double(as.matrix(x)))
        else
            x <- as.matrix(array(as.double(as.matrix(x)), dim(x)))
        
        if (k - ndims - 1 > 0) 
            result <- predictMrtsRcpp(Xu,
                                      xobs_diag,
                                      x,
                                      k - ndims - 1)
        else {
            X2 <- scale(Xu, scale = FALSE)
            shift <- colMeans(Xu)
            nconst <- sqrt(diag(t(X2) %*% X2))
            X2 <- cbind(1, t((t(x) - shift) / nconst) * sqrt(n))
            result <- list(X = X2[, 1:k])
            x <- NULL
        }
    }
    else {
        if (k - ndims - 1 > 0) 
            result <- mrtsRcpp(Xu, xobs_diag, k - ndims - 1)
        else {
            X2 <- scale(Xu, scale = FALSE)
            shift <- colMeans(Xu)
            nconst <- sqrt(diag(t(X2) %*% X2))
            X2 <- cbind(1, t((t(Xu) - shift)/nconst) * sqrt(n))
            result <- list(X = X2[, 1:k])
        }
    }
    obj <- result$X
    attr(obj, "UZ") <- result$UZ
    attr(obj, "Xu") <- Xu
    attr(obj, "nconst") <- result$nconst
    attr(obj, "BBBH") <- result$BBBH
    attr(obj, "class") <- c("matrix", "mrts")
    class(obj) <- "mrts"
    if (is.null(x)) 
        return(obj)
    else {
        shift <- colMeans(attr(obj, "Xu"))
        X2 <- sweep(cbind(x), 2, shift, "-")
        X2 <- cbind(1, sweep(X2, 2, attr(obj, "nconst"), "/"))
        if (k - ndims - 1 > 0)
            obj0 <- as.matrix(cbind(X2, result$X1))
        else
            obj0 <- as.matrix(X2)
        dimnames(obj) <- NULL
        aname <- names(attributes(obj))
        attributes(obj0) <- c(attributes(obj0),
                              attributes(obj)[setdiff(aname, c("dim", "dimnames"))])
        
        return(obj0)
    }
}

predict.FRK <- function(object, obsData = NULL, obsloc = NULL, mu.obs = 0, 
                        newloc = obsloc, basis = NULL, mu.new = 0, se.report = FALSE, 
                        ...) {
    
    if (!"w" %in% names(object)) 
        stop("input model (object) should use the option \"wsave=TRUE\"!")
    if (is.null(basis)) {
        if (is.null(newloc) & is.null(obsloc)) 
            basis <- object$G
        else {
            if (class(object$G) != "mrts") 
                stop("Basis matrix of new locations should be given (unless the model was fitted with mrts bases)!")
            else {
                if (is.null(newloc))
                    basis <- object$G
                else
                    basis <- predict.mrts(object$G, newx = newloc)
            }
        }
    }
    if (NROW(basis) == 1) basis <- as.matrix(t(basis))

    nobs <- ifelse(is.null(obsloc), NROW(object$G), NROW(obsloc))

    if (!is.null(obsData)) {
        obsData <- as.matrix(obsData - mu.obs)
        if (length(obsData) != nobs) 
            stop("Dimensions of obsloc and obsData are not compatible!")
    }
    if (!is.null(newloc)) {
        if (NROW(basis) != NROW(newloc)) 
            stop("Dimensions of newloc and basis are not compatible!")
    }
    else {
        if (NROW(basis) != NROW(object$G)) 
            stop("Dimensions of obsloc and basis are not compatible!")	
    }
    if (is.null(object$LKobj)) {
        if (is.null(obsloc) & is.null(obsData)) {
            miss <- attr(object, "missing")
            yhat <- basis %*% object$w
            if (se.report) {
                TT <- NCOL(object$w)
                if (is.null(miss)) {
                    se <- sqrt(pmax(0, rowSums((basis %*% object$V) * basis)))
                    se <- matrix(se, length(se), TT)
                }
                else {
                    se <- matrix(NA, NROW(basis), TT)
                    pick <- attr(object, "pinfo")$pick
                    D0 <- attr(object, "pinfo")$D[pick, pick]
                    miss <- (as.matrix(miss$miss) == 1)
                    Fk <- object$G[pick, ]
                    M <- object$M
                    dec <- eigen(M)
                    for (tt in 1:TT) {
                        G <- Fk[!miss[, tt], ]
                        GM <- G %*% M
                        De <- D0[!miss[, tt], !miss[, tt]]
                        L <- G %*% dec$vector %*% diag.spam(sqrt(pmax(dec$value, 0)), NROW(M))
                        V <- as.matrix(M - t(GM) %*% invCz(object$s * De, L, GM))
                        se[, tt] <- sqrt(pmax(0, rowSums((basis %*% V) * basis)))
                    }
                }
            }
        }
        if (!is.null(obsData)) {
            pick <- which(!is.na(obsData))
            if (is.null(obsloc)) {
                De <- attr(object, "pinfo")$D[pick, pick]
                G <- object$G[pick, ]
            }
            else {
                De <- diag.spam(length(pick))
                G <- predict.mrts(object$G, newx = as.matrix(obsloc)[pick, ])
            }
            M <- object$M
            GM <- G %*% M
            dec <- eigen(M)
            L <- G %*% dec$vector %*% diag.spam(sqrt(pmax(dec$value, 0)), NROW(M))
            yhat <- basis %*% t(GM) %*% invCz(object$s * De, L, 
                                              obsData[pick])
            if (se.report) {
                V <- as.matrix(M - t(GM) %*% invCz(object$s * De, L, GM))
                se <- sqrt(pmax(0, rowSums((basis %*% V) * basis)))
            }
        }
    }
    else {
        
        LKpeon <- function(M, s, Fk, basis, weight, phi1, phi0, 
                           Q, lambda, phi0P, L = NULL, Data = NULL, only.wlk = FALSE, 
                           only.se = FALSE) {
            wwX <- diag.spam(weight) %*% phi1
            wXiG <- (wwX) %*% solve(t(wwX) %*% phi1 + lambda * Q)
            fM <- Fk %*% M
            if (is.null(L)) {
                dec <- eigen(M)
                L <- Fk %*% dec$vector %*% diag.spam(sqrt(pmax(dec$value, 0)), NROW(M))
                L <- as.matrix(L)
            }
            iDL <- weight * L - wXiG %*% (t(wwX) %*% L)
            iDFk <- weight * Fk - wXiG %*% (t(wwX) %*% as.matrix(Fk))
            itmp <- solve(diag(1, NCOL(L)) + t(L) %*% iDL/s)
            ihL <- chol(itmp) %*% t(L)
            iiLiD <- itmp %*% t(iDL/s)
            
            if (only.wlk) {
                LiiLiDZ <- L %*% (iiLiD %*% Data)
                w <- M %*% t(iDFk) %*% Data/s - (M %*% t(iDFk/s)) %*% (LiiLiDZ)
                wlk <- t(wXiG) %*% Data - t(wXiG) %*% (LiiLiDZ)
                return(list(w = w, wlk = wlk))
            }
            
            MFiS11 <- M %*% t(iDFk)/s - ((M %*% t(iDFk/s)) %*% L) %*% iiLiD
            FMfi <- basis %*% MFiS11
            p0Pp1 <- as.matrix(phi0P %*% t(phi1))
            se0 <- rowSums((basis %*% M) * basis) + rowSums(as.matrix(phi0P * phi0))/lambda * s
            se11 <- rowSums((FMfi %*% fM) * basis)
            se12 <- rowSums(p0Pp1 * (FMfi)) * s/lambda
            se13 <- se12
            se14 <- rowSums(as.matrix(phi0 %*% t(wXiG)) * p0Pp1) * s/lambda - colSums((ihL %*% wXiG %*% t(phi0))^2)
            se <- sqrt(pmax(se0 - (se11 + se12 + se13 + se14), 0))
            if (only.se) 
                return(se)
            else {
                return(list(se = se,
                            w = MFiS11 %*% Data,
                            wlk = t(wXiG) %*% Data - t(wXiG) %*% L %*% (iiLiD %*% Data)))
            }
        }
        if (is.null(obsloc) & is.null(obsData)) {
            if (is.null(newloc)) newloc <- attr(object, "pinfo")$loc
            miss <- attr(object, "missing")
            info <- object$LKobj$LKinfo.MLE
            phi0 <- LKrig.basis(newloc, info)
            pinfo <- attr(object, "pinfo")
            yhat <- basis %*% object$w + phi0 %*% pinfo$wlk
            
            if (se.report) {
                TT <- NCOL(object$w)
                lambda <- object$LKobj$lambda.MLE
                loc <- attr(object, "pinfo")$loc
                pick <- pinfo$pick
                G <- object$G[pick, ]
                M <- object$M
                dec <- eigen(M)
                L <- G %*% dec$vector %*% diag.spam(sqrt(pmax(dec$value, 0)), NROW(M))
                L <- as.matrix(L)
                phi1 <- LKrig.basis(as.matrix(loc)[pick, ], info)
                Q <- LKrig.precision(info)
                weight <- pinfo$weights[pick]
                s <- object$s
                phi0P <- phi0 %*% solve(Q)
                if (is.null(miss)) {
                    se <- LKpeon(M = M,
                                 s = s,
                                 Fk = G,
                                 basis = basis,
                                 weight = weight,
                                 phi1 = phi1,
                                 phi0 = phi0, 
                                 Q = Q,
                                 lambda = lambda, 
                                 phi0P = phi0P,
                                 L = L,
                                 only.se = TRUE)
                    se <- matrix(se, length(se), TT)
                }
                else {
                    se <- matrix(NA, NROW(basis), TT)
                    miss <- (as.matrix(miss$miss) == 1)
                    for (tt in 1:TT) 
                        se[, tt] <- LKpeon(M = M,
                                           s = s,
                                           Fk = G[!miss[, tt], ],
                                           basis = basis,
                                           weight = weight[!miss[, tt]],
                                           phi1 = phi1[!miss[, tt], ],
                                           phi0 = phi0,
                                           Q = Q,
                                           lambda = lambda,
                                           phi0P = phi0P,
                                           L = L[!miss[, tt], ],
                                           only.se = TRUE)
                }
            }
        }
        if (!is.null(obsData)) {
            loc <- attr(object, "pinfo")$loc
            if (is.null(newloc)) 
                newloc <- loc
            pick <- which(!is.na(obsData))
            
            if (is.null(obsloc)) {
                obsloc <- loc
                De <- attr(object, "pinfo")$D[pick, pick]
                G <- object$G[pick, ]
            }
            else
                G <- predict.mrts(object$G, newx = as.matrix(obsloc)[pick, ])
            M <- object$M
            dec <- eigen(M)
            L <- G %*% dec$vector %*% diag.spam(sqrt(pmax(dec$value, 0)), NROW(M))
            L <- as.matrix(L)
            info <- object$LKobj$LKinfo.MLE
            phi1 <- LKrig.basis(as.matrix(obsloc)[pick, ], info)
            Q <- LKrig.precision(info)
            weight <- rep(1, length(pick))
            s <- object$s
            phi0 <- LKrig.basis(newloc, info)
            phi0P <- phi0 %*% solve(Q)
            lambda <- object$LKobj$lambda.MLE
            pred <- LKpeon(M = M,
                           s = s,
                           Fk = G,
                           basis = basis,
                           weight = weight,
                           phi1 = phi1,
                           phi0 = phi0, 
                           Q = Q,
                           lambda = lambda,
                           phi0P = phi0P,
                           L = L,
                           Data = obsData[pick, ], 
                           only.wlk = !se.report)

            yhat <- basis %*% pred$w + phi0 %*% pred$wlk
            if (se.report) se <- pred$se
        }
    }
    if (!se.report) 
        return(list(pred.value = yhat + mu.new, se = NULL))
    else
        return(list(pred.value = yhat + mu.new, se = as.matrix(se)))
}

predict.mrts <- function(object, newx, ...) {
    
    if (missing(newx)) 
        result <- object
    else{
        Xu <- attr(object, "Xu")
        n <- NROW(Xu)
        xobs_diag <- diag(sqrt(n/(n - 1))/apply(Xu, 2, sd), ncol(Xu))
        ndims <- NCOL(Xu)
        k <- NCOL(object)
        x0 <- matrix(as.matrix(newx), ncol = ndims)
        
        shift <- colMeans(attr(object, "Xu"))
        X2 <- sweep(cbind(x0), 2, shift, "-")
        X2 <- cbind(1, sweep(X2, 2, attr(object, "nconst"), "/"))
        
        kstar <- (k - ndims - 1)
        if (kstar <= 0) 
            result <- as.matrix(X2)
        else {
            X1 <- predictMrtsRcppWithBasis(Xu,
                                           xobs_diag,
                                           x0,
                                           attr(object,"BBBH"),
                                           attr(object, "UZ"),
                                           attr(object, "nconst"),
                                           k)$X1
            X1 <- X1[, 1:kstar]
            result <- as.matrix(cbind(X2, X1)) 
        }
    }
    
    return(result)
}

print.FRK <- function(x, ...) {
    attr(x, "pinfo") = NULL
    if (!is.null(x$LKobj)) x$LKobj = x$LKobj$summary
    out <- paste("a ", NROW(x$G), " by ", NCOL(x$G), " mrts matrix", sep = "")
    print(out)
}

print.mrts <- function(x, ...) {
    if (NCOL(x) == 1) out = c(x) else out = x[, 1:NCOL(x)]
    print(out)
}

setNC <- function(z, loc, nlevel) {
    Dimension <- NCOL(loc)
    N <- nrow(z)
    a <- sum(2^(Dimension * (0:(nlevel - 1))))
    NCtest <- (N/a)^(1/Dimension)
    
    return(round(max(4, NCtest)))
}

selectBasis <- function(Data, loc, D = diag.spam(NROW(Data)), maxit = 50, avgtol = 1e-6, 
                        maxK = NULL, Kseq = NULL, method = c("fast", "EM"), n.neighbor = 3, 
                        maxknot = 5000, DfromLK = NULL, Fk = NULL) {
    
    # Remove columnwise and rowwise NAs of Data in order
    if (class(Data) == "numeric") 
        Data <- Data[which(!is.na(Data))]
    else if (class(Data) == "matrix") 
        Data <- Data[, colSums(is.na(Data)) != nrow(Data)]
    else
        stop("Please enter a valid class of Data")
    if (class(Data) == "numeric") Data <- as.matrix(Data)
    #
    # Assume all elements of Data are not NAs
    # TODO: Add protection in `autoFRK ` for detection `length(pick) > 0`
    #
    naDataMatrix <- is.na(Data)
    isWithNA <- sum(naDataMatrix) > 0
    pick <- which(rowSums(naDataMatrix) != 0)
    N <- length(pick)
    if (N == 1) Data <- as.matrix(Data[pick, ]) else Data <- Data[pick, ]
    
    klim <-min(N, round(10 * sqrt(N)))
    if (class(loc) != "matrix") loc <- as.matrix(loc)
    if (N < maxknow)
        knot <- loc[pick, ]
    else 
        knot <- subknot(loc[pick, ], min(maxknot, klim))
    
    if (!is.null(maxK)) 
        maxK <- round(maxK)
    else
        maxK <- ifelse(!is.null(Kseq), round(max(Kseq)), klim)
    
    d <- NCOL(loc)
    TT <- NCOL(Data)
    if (!is.null(Kseq)) {
        K <- unique(round(Kseq))
        if (max(K) > maxK) stop("maximum of Kseq is larger than maxK!")
        if (any(K < (d + 1))) 
            warning("The minimum of Kseq can not less than ", 
                    d + 1, ". Too small values will be ignored.")
        K <- K[K > d]
        if (length(K) == 0) stop("Not valid Kseq!")
    }
    else {
        K <- unique(round(seq(d + 1, maxK, by = maxK^(1/3) * d)))
        if (length(K) > 30) K <- unique(round(seq(d + 1, maxK, l = 30)))
    }
    
    if (is.null(Fk)) Fk <- mrts(knot, max(K), loc)
    AIClist <- rep(Inf, length(K))
    method <- match.arg(method)
    if ((method == "EM") & (is.null(DfromLK))) {
        for (k in 1:length(K)) 
            AIClist[k] <- indeMLE(Data = Data,
                                  Fk = Fk[pick, 1:K[k]],
                                  D = D,
                                  maxit = maxit,
                                  avgtol = avgtol,
                                  num.report = FALSE)$negloglik
    }
    else {
        if (isWithNA) Data <- apply(Data, 2, imputeByKnn , loc=loc, k=n.neighbor)
        
        if (is.null(DfromLK)) {
            iD <- solve(D)
            iDFk <- iD %*% Fk[pick, ]
            iDZ <- iD %*% Data
        }
        else {
            wX <- DfromLK$wX[pick, ]
            G <- t(DfromLK$wX) %*% DfromLK$wX + DfromLK$lambda * DfromLK$Q
            weight <- DfromLK$weights[pick]
            wwX <- diag.spam(sqrt(weight)) %*% wX
            wXiG <- (wwX) %*% solve(G)
            iDFk <- weight * Fk[pick, ] - wXiG %*% (t(wwX) %*% as.matrix(Fk[pick, ]))
            iDZ <- weight * Data - wXiG %*% (t(wwX) %*% as.matrix(Data))
        }
        
        trS <- sum(rowSums(as.matrix(iDZ) * Data))/TT
        for (k in 1:length(K)) {
            half <- getHalf(Fk[pick, 1:K[k]], iDFk[, 1:K[k]])
            ihFiD <- half %*% t(iDFk[, 1:K[k]])
            JSJ <- tcrossprod(ihFiD %*% Data)/TT
            JSJ <- (JSJ + t(JSJ))/2
            AIClist[k] <- cMLE(Fk = Fk[pick, 1:K[k]],
                               TT = TT,
                               trS = trS,
                               half = half,
                               JSJ = JSJ)$negloglik
        }
    }
    
    df <- (K * (K + 1)/2 + 1) * (K <= TT) + (K * TT + 1 - TT * (TT - 1)/2) * (K > TT)
    AIClist <- AIClist + 2 * df
    Kopt <- K[which.min(AIClist)]
    out <- Fk[, 1:Kopt]
    
    dimnames(Fk) <- NULL
    aname <- names(attributes(Fk))
    attributes(out) <- c(attributes(out), attributes(Fk)[setdiff(aname, "dim")])
    
    return(out)
}

subknot <- function(x, nknot, xrng = NULL, nsamp = 1) {
    x <- as.matrix(x)
    xdim <- dim(x)
    
    if (xdim[2] > 1)
        for (kk in xdim[2]:1) 
            x <- x[order(x[, kk]), ]
    else 
        x <- as.matrix(sort(x))
    
    if (is.null(xrng))
        if (xdim[2] > 1)
            xrng <- apply(x, 2, range)
        else
            xrng <- matrix(range(x), 2, 1)
    
    mysamp <- function(zANDid) {
        z <- as.double(names(zANDid))
        if (length(z) == 1L) {
            return(z)
        }
        else {
            set.seed(mean(zANDid))
            return(sample(z, size = min(nsamp, length(z))))
        }
    }
    
    rng <- sqrt(xrng[2, ] - xrng[1, ])
    rng[rng == 0] <- min(rng[rng > 0])/5
    rng <- rng * 10/min(rng)
    nmbin <- round(exp(log(rng) * log(nknot)/sum(log(rng))))
    nmbin <- pmax(2, nmbin)
    
    while (prod(nmbin) < nknot) 
        nmbin[which.max(rng)] <- nmbin[which.max(rng)] + 1
    
    gvec <- matrix(1, xdim[1], 1)
    cnt <- 0
    
    while (length(unique(gvec)) < nknot) {
        nmbin <- nmbin + cnt
        gvec <- matrix(1, xdim[1], 1)
        kconst <- 1
        
        for (kk in 1:xdim[2]) {
            grp <- pmin(round((nmbin[kk] - 1) * ((x[, kk] - xrng[1, kk])/(xrng[2, kk] - xrng[1, kk]))),
                        nmbin[kk] - 1L)
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
    
    if (xdim[2] > 1) x[index, ] else  x[index]
}

toSpMat <- function(mat, MAX_LIMIT = 1e8) {
    if (class(mat) == "data.frame") 
        mat <- as.matrix(mat)
    
    if (class(mat) == "matrix") {
        if (length(mat) > MAX_LIMIT) {
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
        mat <- spam(0, nrow = NROW(mat), ncol = NCOL(mat))
        mat[ridx, cidx] <- mat[where]
    }
    
    return(mat)
}

ZinvC <- function(R, L, z) {
    K <- NCOL(L)
    iR <- solve(R)
    ZiR <- z %*% iR
    left <- ZiR %*% L %*% solve(diag(1, K) + t(L) %*% iR %*% L) %*% t(L)
    
    return(ZiR - left %*% iR)
}

sol.v <- function(d, s, trS, n) {
    # Assume d is in acsending order
    if (max(d) < max(trS/n, s)) 
        result <- max(trS/n - s, 0)
    else{
        estimated_ds <- (trS - cumsum(d)) / (n - 1:length(d))
        L <- max(which(d > estimated_ds))
        result <- ifelse(L == n,
                         max(estimated_ds[L-1] - s, 0), 
                         max(estimated_ds[L] - s, 0))
    }
    
    return(result)
}

sol.eta <- function(d, s, v) pmax(d - s - v, 0)

neg2llik <- function(d, s, v, trS, n) {
    eta <- sol.eta(d, s, v)
    
    if (max(eta/(s + v)) > 1e20) 
        result <- Inf
    else
        result <- n * log(2 * pi) + sum(log(eta + s + v)) + log(s + v) * 
            (n - length(d)) + 1/(s + v) * trS - 1/(s + v) * sum(d * eta/(eta + s + v))
    
    return(result)
}

