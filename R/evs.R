#' Big function for Ensemble Variable Selection using random selection
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel clusterExport
#' @param x data matrix
#' @param y response vector or 'Surv' class object for CoxPH model
#' @param B1 Total number of ensemble iteration
#' @param B2 Total number of iterations for weak learners
#' @param d Random selection dimension
#' @param model Either 'gaussian', 'logistic', 'poisson' or 'coxph'
#' @param tol tolerance level vector, default setting is 0.1, 0.05, 0.01
#' @param parallel Boolean variable for parallel computation
#' @param cl 'makeCluster()' object
#' @param n_cand number of candidate beta coefficient for warm start
#' @param control_opt list variable for hyper-parameters in coxph.fit. Similar to coxph.control()
#' @export evs
evs = function(x, y, B1, B2, d, model = 'gaussian', tol = c(0.1, 0.05, 0.01),
               parallel = FALSE, cl = NULL, n_cand = 10,
               # control opt is for coxph method
               control_opt = list(strata = NULL, init = NULL, offset = NULL,
                                  eps = 1e-9, maxiter=20, weights = NULL,
                                  method = 'efron', rownames = NULL)) {
  n = dim(x)[1]
  p = ncol(x)

  if(!(model %in% c('gaussian', 'logistic', 'poisson', 'coxph'))) {
    stop("Models should be either gaussian, poisson, logistic or coxph")
  }

  if(parallel) {
    #if(is.null(ncore)) {
    #  ncore = floor(0.5 * detectCores()) # default setting = 0.25 * number of cores
    #}
    #    cl = makeCluster(ncore)
    registerDoParallel(cl)
    clusterExport(cl, c("my_mse_index", "logit", "my_logit_loglik", "coxph.fit", "my_pois_loglik",
                        "my_coxph_ploglik","my_optim_logit", "my_optim_pois",
                        "my_optim_coxphfit", "temp_parallel", "Surv"), envir=environment())

    sel.id = foreach(b1 = 1:B1, .combine = 'c') %dopar% {
      temp_parallel(x, y, b1, B2, n, p, d, model)
    }

    #    stopCluster(cl)
  } else {

    sel.id = NULL
    if(model == 'gaussian') {
      #obj_func = my_mse_index
      for(b1 in 1:B1) {
        set.seed(b1)

        id_ind = sapply(1:B2, function(x) sample(1:p, d, replace=FALSE))
        temp_res = apply(id_ind, 2, my_mse_index, x, y)
        sel = which.max(temp_res)
        sel.id = c(sel.id, id_ind[,sel])
      }
    } else if(model == 'poisson') {
      #obj_func = my_optim_pois
      for(b1 in 1:B1) {
        set.seed(b1)

        id_ind = sapply(1:B2, function(x) sample(1:p, d, replace=FALSE))
        temp_res = apply(id_ind, 2, my_optim_pois, x, y, n, d)
        sel = which.max(temp_res)
        sel.id = c(sel.id, id_ind[,sel])
      }
    } else if(model == 'logistic') {
      #obj_func = my_optim_logit
      for(b1 in 1:B1) {
        set.seed(b1)

        id_ind = sapply(1:B2, function(x) sample(1:p, d, replace=FALSE))
        temp_res = apply(id_ind, 2, my_optim_logit, x, y, n, d)
        sel = which.max(temp_res)
        sel.id = c(sel.id, id_ind[,sel])
      }
    } else if(model == 'coxph') {
      #obj_func = my_optim_coxphfit
      for(b1 in 1:B1) {
        set.seed(b1)

        id_ind = sapply(1:B2, function(x) sample(1:p, d, replace=FALSE))
        temp_res = apply(id_ind, 2, my_optim_coxphfit, x, y) # delta removed
        sel = which.max(temp_res)
        sel.id = c(sel.id, id_ind[,sel])
      }
    }
  }

  sel_vec = c(rep(0, p))
  var_sel_id=names(table(sel.id))
  sel_vec[as.integer(var_sel_id)] = table(sel.id)[var_sel_id]

  result = list(sel_vec/B1)
  for(i in 1:length(tol)) {
    result[[i+1]] = BH_method(tol[i], sel_vec)
  }

  names(result) = c("Selection_proportion", noquote(paste("Tolerance", tol[1:length(tol)])))

  return(result)
}
