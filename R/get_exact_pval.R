#' Function for p value of binomial distribution (One(right)-sided)
#'
#' @importFrom stats dbinom pbinom
#' @param value selected counts
#' @param B1 Iteration for Ensemble
#' @param d Dimension for random selection
#' @param p Number of covariates
#' @export get_exact_pval
get_exact_pval = function(value, B1, d, p) {
  return(dbinom(value, B1, d/p)+(1-pbinom(value, B1, d/p, lower.tail = TRUE)))
}
