#' Logit function
#' @param x numerical variable
#' @export logit
logit = function(x) {1 / (1 + exp(-x))}
