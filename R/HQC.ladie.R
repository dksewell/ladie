#' Hannon Quinn Information Criterion for class ladie
#' 
#' HQC is given by \eqn{-2\mbox{loglik}(\hat\theta) + 2 * |\theta|}
#' 
#' @param object object of class ladie
#' @export HQC.ladie
#' @export

HQC.ladie = function(object){
  -2.0 * object$log_likelihood + 2 * log(log(nrow(object$data))) * nrow(object$summary)
}

#' @export
HQC = function(object){
  UseMethod("HQC")
}