#' Partial log-likelihood for coxph model
#' @param beta beta coefficients
#' @param x data matrix
#' @param y 'Surv' class variable
#' @export my_coxph_ploglik
my_coxph_ploglik = function(beta, x, y) {

  index_inc=order(y[,1], decreasing = FALSE) # reorder data as time

  x=x[index_inc,]; time=y[index_inc, 1]; delta=y[index_inc, 2]

  eta=x%*%beta

  ploglik = sum(delta * (eta - log(rev(cumsum(rev(exp(eta)))))))

  return(ploglik)
}
