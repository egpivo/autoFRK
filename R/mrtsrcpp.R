mrtsrcpp <- function(s, xobs_diag, k) {
  computeMrtsRcpp(s, xobs_diag, k)
}

mrtsrcpp_predict0 <- function(s, xobs_diag, s_new, k) {
  predictMrtsRcpp(s, xobs_diag, s_new, k)
}

mrtsrcpp_predict <- function(s, xobs_diag, s_new, BBBH, UZ, nconst, k) {
  predictMrtsRcppWithBasis(s, xobs_diag, s_new, BBBH, UZ, nconst, k)
}
