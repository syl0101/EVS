#' Function for parallel computation
#' @param x data matrix
#' @param y Response variable in vector or 'Surv' class object for CoxPH model
#' @param b1 Iteration for Ensemble
#' @param B2 Total number of iteration for selecting weak learners
#' @param n sample size
#' @param p Number of covariates
#' @param d Dimension for random selection
#' @param model Either 'gaussian', 'logistic', 'poisson' or 'coxph'
#' @export temp_parallel
temp_parallel = function(x, y, b1, B2, n, p, d, model) {

  set.seed(b1)

  id_ind = sapply(1:B2, function(x) sample(1:p, d, replace=FALSE))

  if(model =='gaussian') {
    temp_res = apply(id_ind, 2, my_mse_index, x, y)
  } else if(model == "poisson") {
    temp_res = apply(id_ind, 2, my_optim_pois, x, y, n, d)
  } else if(model == "logistic") {
    temp_res = apply(id_ind, 2, my_optim_logit, x, y, n, d)
  } else if(model == "coxph") {
    temp_res = apply(id_ind, 2, my_optim_coxphfit, x, y) # delta removed
  } else {
    stop("Models should be either gaussian, poisson, logistic or coxph")
  }

  sel = which.max(temp_res)
  return(id_ind[,sel])

}
