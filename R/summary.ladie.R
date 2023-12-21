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
    alpha = 1 - CI_level
    summ$Lower = 
      qnorm(0.5 * alpha,
            summ$`Posterior Median`,
            summ$`Posterior SD`)
    summ$Upper = 
      qnorm(1 - 0.5 * alpha,
            summ$`Posterior Median`,
            summ$`Posterior SD`)
  }
  
  if(print) print(summ)
  invisible(summ)
}
