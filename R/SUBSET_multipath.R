#' SUBSET for multiple pathogens
#' 
#' Implement the SUBSET approach (SUBspace Shrinkage via Exponential Tilting) 
#' towards setting the regression coefficients (sans intercepts) obtained 
#' from multiple pathogens proportional to each other.
#' 
#' @param objects list of objects of class ladie
#' @param paths_by_variable binary 0/1 matrix.  The rows should correspond 1-to-1 with the 
#' covariates (including intercept) for each model.  The columns should correspond 
#' 1-to-1 with the pathogens represented by each of the objects.  If two or more 
#' entries in a row are 1, then these coefficients from the corresponding pathogens 
#' will be shrunk towards each other.  All 0 entries do not undergo any shrinkage estimates 
#' for that covariate.
#' @param nu_max positive numeric giving the maximum level of shrinkage allowed.
#' @param nu positive value giving the exact level of shrinkage.
#' 
#' @returns updated list of objects of class ladie
#' 
#' 
#' @export


SUBSET_multipath = function(objects,
                            paths_by_variable,
                            n_prior_draws = 500,
                            nu_max,
                            nu){
  # Get basic quantities
  
  n_pathogens = length(objects)
  
  P = 
    sapply(objects, function(X) which(X$summary$Variable == "sigma") - 1)
  if(any(diff(P) != 0)) stop("LADIE objects need to use the same covariates")
  P = P[1]
  
  if(missing(paths_by_variable)){
    paths_by_variable = matrix(1L,P,n_pathogens)
    paths_by_variable[1,] = 0L
    warning("No paths_by_variable provided, so shrinking all coefficients towards each other with the exception of the intercept.  For something else, use create_SUBSET_multipath_matrix().")
  }
  
  if(any(sapply(objects,function(X) !("alpha" %in% X$summary$Variable) ))) stop("Exp. DR model not supported at this time")
  Q = (P + 2) * n_pathogens
  sigma_seq = seq(P + 1, Q, by = P + 2)
  alpha_seq = seq(P + 2, Q, by = P + 2)
  
  # Get posterior mean and covariance
  Mean0 = 
    lapply(objects,
           function(X) X$summary$`Posterior Median` ) |>
    unlist()
  Mean0[sigma_seq] = log(Mean0[sigma_seq])
  Mean0[alpha_seq] = log(Mean0[alpha_seq])
  
  Sigma0_inv = Sigma0 = 
    matrix(0.0, (P + 2) * length(objects), (P + 2) * length(objects))
  for(i in 1:length(objects)){
    Sigma0[(P + 2) * (i-1) + 1:(P + 2),(P + 2) * (i-1) + 1:(P + 2)] = 
      objects[[i]]$asymptotic_covariance
    Sigma0_inv[(P + 2) * (i-1) + 1:(P + 2),(P + 2) * (i-1) + 1:(P + 2)] = 
      chol2inv(chol(objects[[i]]$asymptotic_covariance))
  }
  
  # Get prior draws of beta
  rhalft = function(n,df){
    draws = rt(n,df = df)
    
    while(min(draws) < 0){
      negative_index = which(draws < 0)
      draws[negative_index] = 
        rt(length(negative_index),df = df)
    }
    
    draws
  }
  
  prior_draws = 
    do.call(
      cbind,
      lapply(objects,
             function(X) cbind(matrix(rnorm(n_prior_draws * P,
                                            mean = c(X$prior_regr_intercept$mean,
                                                     X$prior_regr_coefs$mean),
                                            sd = c(X$prior_regr_intercept$sd,
                                                   X$prior_regr_coefs$sd)),
                                      n_prior_draws,
                                      P),
                               log(rhalft(n_prior_draws,X$prior_sigma_df)),
                               log(rexp(n_prior_draws,X$prior_alpha_rate)))
             )
      )
  
  # Create projection matrix
  ## Create projection matrix function
  Proj = function(x){
    proj_inner = function(y){
      if("matrix" %in% class(y)){
        return( tcrossprod(y %*% qr.solve(crossprod(y)),
                           y)
        )
      }else{
        return( tcrossprod(y, y) / drop(crossprod(y)) )
      }
    }
    if(isTRUE(class(x) ==  "list")){
      return(Matrix::as.matrix(Matrix::bdiag(lapply(x,proj_inner))))
    }else{
      return(proj_inner(x))
    }
  }
  
  ## Create matrix to span
  I_j_S = function(j,S){
    diag(Q)[,(P + 2) * (S-1) + j,drop=F]
  }
  v_S_e = function(S,j){
    if(length(S) > 0){
      vv = matrix(0.0,n_pathogens,1)
      ee = matrix(0.0,P + 2,1)
      vv[S] = 1.0
      ee[j] = 1.0
      return(vv %x% ee)
    }else{
      return(matrix(0.0,Q,0))
    }
  }
  
  L = matrix(0.0,Q,0)
  for(j in 1:P){
    L = 
      cbind(L,
            I_j_S(j,which(paths_by_variable[j,] != 1)),
            v_S_e(which(paths_by_variable[j,] == 1),j))
  }
  L = 
    cbind(L,
          I_j_S(P + 1, 1:n_pathogens),
          I_j_S(P + 2, 1:n_pathogens))
  
  ## Now create the projection matrix
  Proj_mat = Proj(L)
  I_m_P = diag(Q) - Proj_mat
  
  
  # Get optimal nu
  Sigma0_logdet = determinant(Sigma0)$modulus
  if(missing(nu)){
    
    wtilde_1 = 
      sapply(1:nrow(prior_draws),
             function(i){
               exp(-0.5 * tcrossprod(prior_draws[i,],prior_draws[i,] %*% I_m_P))
             })
    
    non_nu_term = 
      -0.5 * Sigma0_logdet -
      0.5 * drop( crossprod(Mean0,Sigma0_inv %*% Mean0) )
    
    get_BF = function(nu){
      Z_nu = mean(sapply(wtilde_1,function(x)x^nu))
      
      Sigma_inv = Sigma0_inv + nu * I_m_P
      Sigma =  chol2inv(chol(Sigma_inv))
      Mean = drop(Sigma %*% Sigma0_inv %*% Mean0)
      
      non_nu_term -
        0.5 * determinant(Sigma_inv)$modulus +
        0.5 * drop( crossprod(Mean,Sigma_inv %*% Mean) ) - 
        log(Z_nu)
    }
    
    if(missing(nu_max)) nu_max = Inf
    lower_bound = 0
    nu = upper_bound = min(5,nu_max/10)
    safety = 0
    
    while( ( abs(nu - upper_bound) / upper_bound < 1e-3) & 
           (safety < 25) &
           (upper_bound < nu_max) ){
      safety = safety + 1
      upper_bound = min(nu_max,2 * upper_bound)
      opt = optimize(get_BF,
                     interval = c(lower_bound,upper_bound),
                     maximum = TRUE)
      lower_bound = upper_bound
      nu = opt$maximum
      BF_nu_0 = exp(opt$objective)
    }
    
    Sigma_inv = Sigma0_inv + nu * I_m_P
    Sigma =  chol2inv(chol(Sigma_inv))
    Mean = drop(Sigma %*% Sigma0_inv %*% Mean0)
    
  }else{
    Sigma_inv = Sigma0_inv + nu * I_m_P
    Sigma =  chol2inv(chol(Sigma_inv))
    Mean = drop(Sigma %*% Sigma0_inv %*% Mean0)
    
    # BF_nu_0 =
    #   exp(
    #     -0.5 * Sigma0_logdet -
    #       0.5 * drop( crossprod(Mean0,Sigma0_inv %*% Mean0) ) -
    #       0.5 * determinant(Sigma_inv)$modulus +
    #       0.5 * drop( crossprod(Mean,Sigma_inv %*% Mean) ) - 
    #       log(mean(sapply(wtilde_1,function(x)x^nu)))
    #   )
  }
  
  
  # Return results
  for(j in 1:n_pathogens){
    objects[[j]]$summary$`Posterior Median` = 
      Mean[(P + 2)*(j-1) + 1:(P + 2)]
    objects[[j]]$summary$`Posterior SD` = 
      sqrt(diag(Sigma[(P + 2)*(j-1) + 1:(P + 2),
                      (P + 2)*(j-1) + 1:(P + 2)]))
    objects[[j]]$summary$Lower = 
      qnorm(0.5 * (1 - objects[[j]]$CI_level),
            objects[[j]]$summary$`Posterior Median`,
            objects[[j]]$summary$`Posterior SD`)
    objects[[j]]$summary$Upper = 
      qnorm(objects[[j]]$CI_level + 0.5 * (1 - objects[[j]]$CI_level),
            objects[[j]]$summary$`Posterior Median`,
            objects[[j]]$summary$`Posterior SD`)
    objects[[j]]$summary$`Prob. Direction`[1:P] = 
      pnorm(0,
            objects[[j]]$summary$`Posterior Median`[1:P],
            objects[[j]]$summary$`Posterior SD`[1:P])
    objects[[j]]$summary$`Prob. Direction`[1:P] =
      sapply(objects[[j]]$summary$`Prob. Direction`[1:P],function(x) max(x,1-x))
    
    objects[[j]]$asymptotic_covariance = 
      Sigma[(P + 2)*(j-1) + 1:(P + 2),
            (P + 2)*(j-1) + 1:(P + 2)]
  }
  
  objects$nu = nu
  
  
  return(objects)
}
