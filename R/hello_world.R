#' @title Hello world
#'
#' @useDynLib AwesomePackage
#'
#' @importFrom Rcpp evalCpp
#'
#' @description Hello world and add "AwesomePackage" to NAMESPACE.
#'
#' @return A string "Data science is fantastic!".
#'
#' @export
#'
#' @examples hello_world()
hello_world <- function ()
{
  rcpp_hello_world()
}
