#' One-step optimization for poisson regression model
#'
#' @importFrom stats runif
#' @param index variable index
#' @param x data matrix
#' @param y non-negative response vector (counts)
#' @param n sample size
#' @param d number of random selection dimension
#' @param n_cand number of candidate beta coefficient for warm start
#' @export my_optim_pois
my_optim_pois = function(index, x, y, n, d, n_cand=10) {
  # generating random candidate beta for warm-start
  # non-informative candidates

  beta_init = sapply(1:n_cand, function(k) runif(d, -1, 1))

  dif_mat = exp(x[,index] %*% beta_init) - matrix(y, n, n_cand)

  id_beta_subopt = which.min(crossprod(dif_mat^2, rep(1, n)))
  beta_subopt = beta_init[,id_beta_subopt]

  eta = x[, index] %*% beta_subopt
  mu = c(exp(eta))
  z = eta + (y-mu)/mu

  tilde.x = x[, index] * sqrt(mu)
  tilde.z = z * sqrt(mu)
  qr.obj = qr(tilde.x)
  beta_update <- backsolve(qr.obj$qr, qr.qty(qr.obj, tilde.z))

  loglik_1_step = my_pois_loglik(beta_update, x[,index], y)
  return(loglik_1_step)
}

