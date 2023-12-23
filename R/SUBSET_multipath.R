#' SUBSET for multiple pathogens
#' 
#' Implement the SUBSET approach (SUBspace Shrinkage via Exponential Tilting) 
#' towards setting the regression coefficients (sans intercepts) obtained 
#' from multiple pathogens proportional to each other.
#' 
#' 
#' @import SUBSET


SUBSET_multipath = function(objects,
                            prior_phi,
                            n_posterior_draws = 5e3,
                            n_prior_draws = 500,
                            verbose = TRUE,
                            nu_max){
  # Get basic quantities
  n_pathogens = length(objects)
  
  P = 
    sapply(objects, function(X) which(X$summary$Variable == "sigma") - 1)
  if(any(diff(P) != 0)) stop("LADIE objects need to use the same covariates")
  P = P[1]
  Q = (P-1) * n_pathogens
  
  # Get posterior mean and covariance
  Mean0 = 
    lapply(objects,
           function(X) X$summary$`Posterior Median`[2:P] ) |>
    unlist()
  
  
  Sigma0_inv = Sigma0 = 
    matrix(0.0, (P-1) * length(objects), (P-1) * length(objects))
  for(i in 1:length(objects)){
    Sigma0[(P-1) * (i-1) + 1:(P-1),(P-1) * (i-1) + 1:(P-1)] = 
      objects[[i]]$asymptotic_covariance[2:P,2:P]
    Sigma0_inv[(P-1) * (i-1) + 1:(P-1),(P-1) * (i-1) + 1:(P-1)] = 
      chol2inv(chol(objects[[i]]$asymptotic_covariance[2:P,2:P]))
  }
  
  # Get prior draws of beta
  prior_draws = 
    do.call(
      cbind,
      lapply(objects,
             function(X) matrix(rnorm(n_prior_draws * (P-1),
                                      mean = X$prior_regr_coefs$mean,
                                      sd = X$prior_regr_coefs$sd),
                                n_prior_draws,
                                P - 1))
      )
  
  # Create projection matrix function
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
  
  
  # Create P(phi) function
  Proj_phi = function(phi){
    Proj(kronecker(matrix(c(1,phi),n_pathogens,1), diag(P - 1)))
  }
  
  # Get prior on phi if missing
  if(missing(prior_phi)){
    prior_phi = function() rgamma(n_pathogens - 1,
                                   shape = 3,
                                   rate = 2)
  }
  
  
  # Create array to store samples
  samples = 
    matrix(0.0,
           n_posterior_draws,P * n_pathogens - 1,
           dimnames = list(NULL,
                           c(paste(rep(paste("pathogen",1:n_pathogens,sep=""),
                                       each = P-1),
                                   rep(objects[[1]]$summary$Variable[2:P],
                                       n_pathogens),sep = "_"),
                             paste("phi",2:n_pathogens - 1,sep = "_"))))
  # if(n_pathogens > 2){
  #   samples[1,] = c(Mean0,rowMeans(replicate(250,{prior_phi()})))
  # }else{
  #   samples[1,] = c(Mean0,mean(replicate(250,{prior_phi()})))
  # }
  samples[1,] = c(Mean0,rep(1.0,n_pathogens - 1))
  
  
  # Get nu via Bayes factors
  Sigma0_logdet = determinant(Sigma0)$modulus
  {
    Proj_mat = Proj_phi(samples[1,-c(1:Q)])
    I_m_P = diag(Q) - Proj_mat
    
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
      upper_bound = min(nu_max,2 * upper_bound)
      opt = optimize(get_BF,
                     interval = c(lower_bound,upper_bound),
                     maximum = TRUE)
      lower_bound = upper_bound
      nu = opt$maximum
    }
  }
  
  
  # Get exact (actually, this is MC) method for computing Z_{\nu,\phi}
  Z_exact = function(I_m_P){
    mean(
      sapply(1:nrow(prior_draws), 
             function(i){
               exp(-0.5 * nu * drop(crossprod(prior_draws[i,],I_m_P %*% prior_draws[i,])))
             })
    )
  }
  
  phi_prior_draws = replicate(1e3,{prior_phi()})
  if(n_pathogens > 2){
    phi_range = apply(phi_prior_draws,1,range)
    
    phi_seq = 
      lapply(1:ncol(phi_range),function(i) seq(phi_range[1,i],
                                               phi_range[2,i],
                                               l = 10)) |>
      expand.grid() |>
      as.matrix()
    
    Z_vals = 
      sapply(1:nrow(phi_seq),
             function(i){
               Z_exact(diag(Q) - Proj_phi(phi_seq[i,]))
               })
    ns_fit = 
      lm(as.formula(paste0(
        "Z_vals ~ ",
        paste("splines::ns(",
              paste("Var",1:(n_pathogens-1),sep=""),
              ",df = 5)",
              collapse = " : ")
        )),
        data = as.data.frame(cbind(phi_seq,Z_vals)))
    
    Z = function(x){
      x = as.data.frame(matrix(x,nr=1))
      colnames(x) = paste("Var",1:(n_pathogens-1),sep="")
      predict(ns_fit,
              newdata = x)
    }
    
  }else{
    phi_range = range(phi_prior_draws)
    
    phi_seq = seq(Z_approximation$phi_range[1],
                  Z_approximation$phi_range[2],
                  l = 30)
    
    Z_vals = 
      sapply(phi_seq,function(phi){
        Z_exact(diag(p) - P(phi))
      })
    
    ns_fit = 
      lm(Z ~ splines::ns(phi,df = 11),
         data = data.frame(Z = Z_vals,
                           phi = phi_seq))
    
    Z = function(x){
      predict(ns_fit,
              newdata = data.frame(phi = x))
    }
    
  }
  
  
  P_old = Proj_phi(samples[1,-c(1:Q)])
  Z_old = Z(samples[1,-c(1:Q)])
  acc_rate = 0.0
  cat("\nPerforming Gibbs Sampling\n")
  
  if(verbose) pb = txtProgressBar(0,n_posterior_draws,style=3)
  
  for(it in 2:n_posterior_draws){
    
    # Draw phi from prior
    phi_proposal = prior_phi()
    P_new = Proj_phi(phi_proposal)
    # P_orth_new = diag(p) - P_new
    Z_new = Z(phi_proposal)
    acc_prob = 
      exp(-0.5 * nu * 
            drop(crossprod(samples[it - 1,1:Q],(P_old - P_new) %*% samples[it - 1,1:Q])) ) * 
      ifelse(max(Z_vals) > 1e-15,Z_old / Z_new, 1) # This ifelse is for numerical stability in case the Z_vals are too tiny
    
    if(runif(1) < acc_prob){
      samples[it,-c(1:Q)] = phi_proposal
      P_old = P_new
      Z_old = Z_new
      
      acc_rate = acc_rate + 1 / (n_posterior_draws - 1)
    }else{
      samples[it,-c(1:Q)] = samples[it - 1,-c(1:Q)]
    }
    
    # Draw theta
    Sigma_inv = Sigma0_inv + nu * (diag(Q) - P_old)
    Sigma =  chol2inv(chol(Sigma_inv))
    Mean = drop(Sigma %*% Sigma0_inv %*% Mean0)
    samples[it,1:Q] = 
      drop(
        mvtnorm::rmvnorm(1,
                         mean = Mean,
                         sigma = Sigma)
      )
    
    if(verbose) setTxtProgressBar(pb,it)
  }
  
  # asdf This is where I left off.
  
  
}
