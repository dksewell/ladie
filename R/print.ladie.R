#' Print LADIE model output
#' 
#' print method for class "ladie"
#' 
#' @param object
#' 
#' @export print.ladie
#' @exportS3Method ladie::print

print.ladie = function(object){
  cat("LADIE: LAtent Dose Incidence Estimation\n\n")
  
  
  print(object$summary)
  
  cat("NOTES: 1. Intercept is confounded with dose-response model parameter.\n")
  cat("       2. Exponentiate regression coefficients to get rate of accruing dose due to a unit change in 'x'.\n")
  cat(paste0("       3. CI's are given at ",100 * object$CI_level,"%.  Use summary() to set the CI level differently.\n"))
  cat("       4. SD of sigma (and alpha, if beta-poisson model was used) is on the log scale; everything else is at the original scale.\n")
}
