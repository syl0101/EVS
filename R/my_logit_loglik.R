#' Log-likelihood for logistic regression model
#' @param beta beta coefficients
#' @param x data matrix
#' @param y binary response variable
#' @export my_logit_loglik
my_logit_loglik = function(beta, x, y) {
  eta = x %*% beta
  mean(y * log(logit(eta)) + (1-y) * log(1-logit(eta)))
}
