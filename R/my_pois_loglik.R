#' Log-likelihood for poisson regression model
#' @param beta beta coefficients
#' @param x data matrix
#' @param y non-negative response variable (counts)
#' @export my_pois_loglik
my_pois_loglik = function(beta, x, y) {
  eta = x %*% beta
  mean(y * eta - exp(eta))
}
