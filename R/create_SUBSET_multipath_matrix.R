#' Helper function for SUBSET_multipath
#' 
#' This function is designed to help individuals create a matrix for the 
#' argument paths_by_variable in SUBSET_multipath via prompts.
#' 
#' @export


create_SUBSET_multipath_matrix = function(objects){
  
  n_pathogens = length(objects)
  
  P = 
    sapply(objects, function(X) which(X$summary$Variable == "sigma") - 1)
  if(any(diff(P) != 0)) stop("LADIE objects need to use the same covariates")
  P = P[1]
  
  variable_names = objects[[1]]$summary$Variable[1:P]
  
  paths_by_variable = matrix(0L,length(variable_names),n_pathogens)
  
  for(j in 1:length(variable_names)){
    cat("\n\n\n")
    input = 
      eval(parse(text = readline(paste0("Which pathogens should have their coefficients for",
                                        variable_names[j]," shrunk towards each other?  "))))
    if(!is.null(input)) paths_by_variable[j,input] = 1L
  }
  
  return(paths_by_variable)
}
