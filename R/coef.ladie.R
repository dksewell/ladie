#' Extract model coefficients
#' 
#' coef method for class "ladie"
#' 
#' @param object
#' 
#' @export coef.ladie
#' @exportS3Method ladie::coef

coef.ladie = function(object){
  ret = object$summary$`Posterior Mean`
  names(ret) = object$summary$Variable
  
  return(ret)
}
