#' Function for selecting variable under the tolerance level of FDR (False Discovery Rate)
#' @description Using the Idea of Benjamini and Hochberg method.
#'
#' @param tol tolerance level vector or single element
#' @param selection_total selection counts vector of size 'p'
#' @param B1 Total number of ensemble
#' @param d random seletion dimension
#' @param p total number of covariates
#' @export BH_method
BH_method = function(tol, selection_total, B1, d, p) {

  pval = get_exact_pval(selection_total, B1, d, p)
  pval_sorted = sort(pval, decreasing = FALSE)

  sel_crit = max(which(pval_sorted <= c(1:p)/p * tol))

  sel_index = which(rank(pval) <= sel_crit)
  return(sel_index)
}
