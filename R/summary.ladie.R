#' Summarizing LADIE model fits
#' 
#' summary method for class "ladie"
#' 
#' @param object object of class "ladie"
#' @param CI_level the level of credible intervals to be provided
#' @param print logical
#' 
#' @export summary.ladie
#' @exportS3Method ladie::summary


summary.ladie = function(object,CI_level,print = TRUE){
  summ = object$summary
  if(!missing(CI_level)){
    is_quantile = function(x,w,prob = 0.025){
      w = cumsum(w[order(x)])
      x = x[order(x)]
      quant = x[min(which(w >= prob))]
      return(quant)
    }
    
    summ$Lower = 
      apply(object$importance_samples[,2:ncol(object$importance_samples) - 1],2,
            is_quantile, 
            w = object$importance_samples[,"IS_weights"], 
            prob = (1-CI_level)/2)
    summ$Upper = 
      apply(object$importance_samples[,2:ncol(object$importance_samples) - 1],2,
            is_quantile, 
            w = object$importance_samples[,"IS_weights"], 
            prob = 1 - (1-CI_level)/2)
  }
  
  if(print) print(summ)
  invisible(summ)
}
