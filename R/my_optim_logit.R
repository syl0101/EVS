#' One-step optimization for logistic regression model
#'
#' @importFrom stats runif
#' @param index variable index
#' @param x data matrix
#' @param y binary response variable (zero-one)
#' @param n sample size
#' @param d number of random selection dimension
#' @param n_cand number of candidate beta coefficient for warm start
#' @export my_optim_logit
my_optim_logit = function(index, x, y, n, d, n_cand=10) {

  # generating random candidate beta for warm-start
  # non-informative candidates
  beta_init = sapply(1:n_cand, function(k) runif(d, -1, 1))

  dif_mat = logit(x[,index] %*% beta_init) - matrix(y, n, n_cand)

  id_beta_subopt = which.min(crossprod(dif_mat^2, rep(1, n)))
  beta_subopt = beta_init[,id_beta_subopt]

  # among the beta candidates, one that yielded the minimized
  # residuals = (y - hat of y)^2 are selected

  # One step update

  eta = x[, index] %*% beta_subopt
  pi = logit(eta)
  hess = c(pi * (1 - pi))
  z = eta + (y-pi)/hess

  tilde.x = x[, index] * sqrt(hess)
  tilde.z = z * sqrt(hess)
  qr.obj = qr(tilde.x)
  beta_update <- backsolve(qr.obj$qr, qr.qty(qr.obj, tilde.z))

  loglik_1_step = my_logit_loglik(beta_update, x[,index], y)
  return(loglik_1_step)
}
