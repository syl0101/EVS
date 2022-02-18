#' One-step optimization for coxph model
#'
#' @importFrom survival coxph.fit
#' @param index variable index
#' @param x data matrix
#' @param y 'Surv' class variable
#' @param control_opt list variable for hyper-parameters in coxph.fit. Similar to coxph.control()
#' @export my_optim_coxphfit
my_optim_coxphfit = function(index, x, y,
                             control_opt = list(strata = NULL, init = NULL, offset = NULL, eps = 1e-9, maxiter=20,
                                                weights = NULL, method = 'efron', rownames = NULL)) {

  coxph.fit(x[,index], y, strata = control_opt$strata,
            offset = control_opt$offset,
            init = control_opt$init,
            control = list(eps = control_opt$eps, iter.max = control_opt$maxiter),
            weights = control_opt$weights,
            method = control_opt$method,
            rownames = control_opt$rownames)$loglik[2]
}
