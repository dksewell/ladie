#' Akaike Information Criterion for class ladie
#' 
#' AIC is given by \eqn{-2\mbox{loglik}(\hat\theta) + 2 * |\theta|}
#' 
#' @param object object of class ladie
#' @export AIC.ladie
#' @export

AIC.ladie = function(object){
  -2.0 * object$log_likelihood + 2 * nrow(object$summary)
}
