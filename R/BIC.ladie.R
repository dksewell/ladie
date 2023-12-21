#' Bayesian Information Criterion for class ladie
#' 
#' BIC is given by \eqn{-2\mbox{loglik}(\hat\theta) + \log(n) * |\theta|}
#' 
#' @param object object of class ladie
#' @export BIC.ladie
#' @export

BIC.ladie = function(object){
  -2.0 * object$log_likelihood + nrow(object$data) * nrow(object$summary)
}
