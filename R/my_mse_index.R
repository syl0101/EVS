#' -1 * MSE of linear regression model
#' @param index index vector
#' @param x data matrix
#' @param y continuous response vector
#' @export my_mse_index
my_mse_index = function(index, x, y) {
  -mean((qr.resid(qr(x[,index]), y))^2)
}
