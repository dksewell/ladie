#' ladie: LAtent Dose Incidence Estimation
#' 
#' Use longitudinal data to estimate risk ratios of accruing dose.  ladie_boot uses 
#' Bayesian bootstrapping to perform estimation.
#' 
#' @param formula an object of class "formula" (or one that can be coerced to 
#' that class): a symbolic description of the model to be fitted. 
#' @param data a data.frame containing the variables in the model.
#' @param dose_response The dose-response model to be used.  Either the 
#' beta-poisson or exponential. 
#' @param prior_regr_coefs a named list giving the hyperparameters of the 
#' regression coefficients, excluding the intercept.  Should include the following:
#' \itemize{
#'  \item mean - mean of the normal prior
#'  \item sd - standard deviation of the norma prior
#'  \item autoscale - logical used to determine whether the prior sd should be 
#'  rescaled according to the standard deviation of the columns of the design matrix
#'  }
#' @param prior_regr_intercept named list, same as prior_regr_coefs, except no autoscale
#' @param prior_sigma_df positive number giving the degrees of freedom for the half-t
#'  prior on the standard deviation of the mean dose distribution
#' @param prior_alpha_rate positive value giving the rate of the exponential 
#'  prior on alpha (for beta-poisson model only)
#' @param nonlinear_time logical.  Set to TRUE if the expected dose should be the rate 
#'  times time raised to some power to be estimated. NOTE: this loses the interpretability 
#'  of the regression coefficients somewhat.
#' @param CI_level level for the credible intervals.
#' @param n_bootstraps integer. Number of Bayesian bootstrap samples
#' @param verbose logical. Should progress be displayed?
#' @param cluster integer giving the number of threads to use, or more directly an 
#'  object of class "SOCKcluster" from parallel::makeCluster().
#'  
#' @returns An object of class "ladie", which is really a list with the following 
#'  named elements:
#' \itemize{
#'   \item summary - data.frame giving the posterior median, the posterior SD, 
#'    the credible intervals (at level CI_level), and the probability of direction.
#'   \item CI_level
#'   \item log_posterior - value of the log posterior evaluated at the MAP.  
#'   Useful for information criteria.
#'   \item asymptotic_covariance - Full posterior covariance matrix.
#' }
#'  
#' 
#' @export 
#' @import dplyr
#' @import parallel
#' @import Matrix

ladie_boot = function(formula,
                      data,
                      dose_response = c("beta-poisson","exponential")[1],
                      prior_regr_coefs = list(mean = 0, sd = 2.5, autoscale = TRUE),
                      prior_regr_intercept = list(mean = 0, sd = 10),
                      prior_sigma_df = 3,
                      prior_alpha_rate = 1,
                      nonlinear_time = FALSE,
                      CI_level = 0.95,
                      n_bootstraps = 250,
                      verbose = TRUE,
                      cluster){
  
  dose_response =
    match.arg(dose_response,
              c("beta-poisson","exponential","simple_threshold"))
  
  # Set up parallelization
  if(!missing(cluster)){
    if(is.numeric(cluster)){
      if(verbose) cat("\nSetting up parallel environment")
      cluster = makeCluster(min(detectCores(),cluster))
      on.exit({stopCluster(cluster)})
    }
  }
  
  
  # Drop NAs
  data_clean = 
    data |>
    dplyr::select(dplyr::all_of(all.vars(formula))) |>
    na.omit()
  
  
  # Extract elements from formulas
  {
    varnames = list()
    formula_string = as.character(formula)
    ## Get y variables
    varnames$y = 
      formula_string[[2]]
    ## Get ID
    varnames$id = 
      substr(formula_string[[3]],
             gregexpr("\\|",formula_string[[3]])[[1]] + 1,
             gregexpr("\\)",formula_string[[3]])[[1]] - 1) |> 
      trimws()
    ## Get time variable
    varnames$time = 
      substr(formula_string[[3]],
             gregexpr("\\(",formula_string[[3]])[[1]] + 1,
             gregexpr("\\|",formula_string[[3]])[[1]] - 1) |> 
      trimws()
    ## Get covariates
    plus_locations = 
      c(0,gregexpr("\\+",formula_string[[3]])[[1]])
    varnames$covariates = 
      sapply(1:(length(plus_locations) - 1),
             function(i){
               substr(formula_string[[3]],
                      plus_locations[i] + 1,
                      plus_locations[i + 1] - 1) |> 
                 trimws()
             })
    rm(plus_locations,formula_string)
  }
  
  
  # Get time variable shifted to start at zero
  {
    data_clean =
      data_clean |>
      dplyr::group_by(id) |> 
      dplyr::arrange(pick(varnames$time)) |> 
      dplyr::mutate(time_diff = time - min(time)) |>
      dplyr::ungroup() |>
      dplyr::arrange(pick(varnames$id,"time_diff")) |>
      dplyr::relocate(all_of(c(varnames$id,varnames$time,"time_diff",varnames$y)))
  }
  
  
  # Filter out repeated positives
  {
    helper_function = function(x){
      ret = 1
      n = length(x)
      if(n > 1){
        ret = 
          c(ret,sapply(2:n,function(i) prod(1 - x[1:(i-1)])))
      }
      ret[1] = 0
      return(ret)
    }
    data_clean = 
      data_clean |>
      rename(!!"y" := varnames$y)
    data_clean = 
      data_clean |>
      dplyr::group_by(id) |> 
      dplyr::mutate(S_it = helper_function(y)) |> 
      dplyr::ungroup() |> 
      dplyr::filter(time_diff > 0,
                    S_it == 1) |> 
      dplyr::select(-S_it)
    rm(helper_function)
  }
  
  
  # Get regression quantities
  data_clean = 
    data_clean |>
    dplyr::mutate(log_time_diff = log(time_diff))
  if(isTRUE(nonlinear_time)){
    X = 
      model.matrix(as.formula(paste0("y ~ log_time_diff + ",
                                     paste(varnames$covariates,
                                           collapse = "+"))),
                   data = data_clean)
  }else{
    X = 
      model.matrix(as.formula(paste0("y ~ ",
                                     paste(varnames$covariates,
                                           collapse = "+"))),
                   data = data_clean)
  }
  P = ncol(X)
  
  # Adjust prior hyperparameters if autoscale = T
  if(isTRUE(prior_regr_coefs$autoscale)){
    sd_x = apply(X[,-1],2,sd)
    names(sd_x) = colnames(X)[-1]
    prior_regr_coefs$sd = 
      prior_regr_coefs$sd / sd_x
    rm(sd_x)
  }
  
  
  # Split here based on the two different dose-response models
  if(dose_response == "beta-poisson"){
    
    ## Create negative log posterior function
    nlpost = function(x){
      beta = x[1:P]
      sig = exp(x[P+1])
      a = exp(x[P+2])
      
      mu_i = X %*% beta
      if(!nonlinear_time){
        mu_i = mu_i + data_clean$log_time_diff
      }
      
      p_i = 
        sapply(1:nrow(data_clean),
               function(i){
                 helper = function(dummy){
                   (1.0 - (1.0 + exp(mu_i[i] + dummy * sig))^(-a)) * dnorm(dummy)
                 }
                 integrate(helper,-4,4)$value
               })
      
      -sum(dbinom(data_clean$y,1,p_i,log = T)) -
        dnorm(beta[1],
              prior_regr_intercept$mean,
              prior_regr_intercept$sd,
              log = TRUE) -
        sum(dnorm(beta[-1],
                  prior_regr_coefs$mean,
                  prior_regr_coefs$sd,
                  log = TRUE)) - 
        dt(sig,df = prior_sigma_df, log = TRUE) - log(sig) -
        dexp(a,rate = prior_alpha_rate, log = TRUE) - log(a)
    }
    
    ## Get initial values
    if(verbose) cat("\nGetting initial values...")
    if(isTRUE(nonlinear_time)){
      initial_fit = 
        glm(as.formula(paste0("y ~ log_time_diff + ",
                              paste(varnames$covariates,
                                    collapse = "+"))),
            data = data_clean,
            family = binomial(link = "cloglog"))
    }else{
      initial_fit = 
        glm(as.formula(paste0("y ~ offset(log_time_diff) + ",
                              paste(varnames$covariates,
                                    collapse = "+"))),
            data = data_clean,
            family = binomial(link = "cloglog"))
    }
    inits = coef(initial_fit)
    partial_nlpost = function(x){
      nlpost(c(x[1],inits[-1],x[-1]))
    }
    initial_opt = 
      optim(c(inits[1],0,0),
            partial_nlpost,
            method = "Nelder-Mead",
            control = list(maxit = 200))
    inits = c(initial_opt$par[1],inits[-1],initial_opt$par[-1])
    
    ## Fit full model
    if(verbose) cat("\nPerforming optimization...")
    fit = 
      optim(inits,
            nlpost,
            method = "BFGS")
    
    ## Perform Bayesian bootstrapping
    if(verbose) cat("\nPerforming Bayesian bootstrapping...\n")
    
    IDs = unique(data_clean[[varnames$id]])
    
    ID_indices = 
      lapply(IDs, function(id) which(data_clean[[varnames$id]] == id))
    
    boot_helper = function(i){
      
      all_levels_represented = FALSE
      while(!all_levels_represented){
        boot_probabilities = 
          drop(extraDistr::rdirichlet(1,rep(1.0,length(IDs))))
        IDs_boot = sample(1:length(IDs),length(IDs),replace=T, prob = boot_probabilities)
        
        data_boot = 
          data_clean[unlist(ID_indices[IDs_boot]),]
        
        n_unique = 
          data_boot |> 
          dplyr::select(all_of(varnames$covariates)) |>
          dplyr::summarize_all(function(x){length(unique(x))})
        
        all_levels_represented = 
          !any(n_unique == 1)
      }
      
      if(isTRUE(nonlinear_time)){
        X_boot = 
          model.matrix(as.formula(paste0("y ~ log_time_diff + ",
                                         paste(varnames$covariates,
                                               collapse = "+"))),
                       data = data_boot)
      }else{
        X_boot = 
          model.matrix(as.formula(paste0("y ~ ",
                                         paste(varnames$covariates,
                                               collapse = "+"))),
                       data = data_boot)
      }
      
      nlpost_boot = function(x){
        beta = x[1:P]
        sig = exp(x[P+1])
        a = exp(x[P+2])
        
        mu_i = X_boot %*% beta
        if(!nonlinear_time){
          mu_i = mu_i + data_boot$log_time_diff
        }
        
        p_i = 
          sapply(1:nrow(data_boot),
                 function(i){
                   helper = function(dummy){
                     (1.0 - (1.0 + exp(mu_i[i] + dummy * sig))^(-a)) * dnorm(dummy)
                   }
                   integrate(helper,-4,4)$value
                 })
        
        -sum(dbinom(data_boot$y,1,p_i,log = T)) -
          dnorm(beta[1],
                prior_regr_intercept$mean,
                prior_regr_intercept$sd,
                log = TRUE) -
          sum(dnorm(beta[-1],
                    prior_regr_coefs$mean,
                    prior_regr_coefs$sd,
                    log = TRUE)) - 
          dt(sig,df = prior_sigma_df, log = TRUE) - log(sig) -
          dexp(a,rate = prior_alpha_rate, log = TRUE) - log(a)
      }
      
      fit_boot = NA
      try({
        fit_boot = 
          optim(fit$par,
                nlpost_boot,
                method = "BFGS")$par
      },silent=T)
      
      return(fit_boot)
    }
    
    if(!missing(cluster)){
      clusterExport(cluster,
                    c("P","IDs","ID_indices","data_clean","nonlinear_time",
                      "varnames","prior_regr_intercept","prior_regr_coefs",
                      "prior_sigma_df","prior_alpha_rate","fit"),
                    envir = environment())
      bootstraps = 
        parSapply(cluster,
                  1:n_bootstraps,
                  boot_helper)
    }else{
      bootstraps = 
        sapply(1:n_bootstraps,
               boot_helper)
    }
    
    ## Return results
    results = list()
    
    results$summary = 
      data.frame(Variable = c(colnames(X),"sigma","alpha"),
                 `Posterior Mean` = rowMeans(bootstraps),
                 `Posterior Mode` = fit$par,
                 `Posterior SD` = 
                   apply(bootstraps,1,sd),
                 Lower = 
                   apply(bootstraps,1,quantile,probs = 0.025),
                 Upper = 
                   apply(bootstraps,1,quantile,probs = 0.975),
                 `Prob. Direction` = 
                   c(apply(bootstraps[1:P,],1,function(x) mean(x > 0)),NA,NA),
                 check.names = FALSE,
                 row.names = NULL)
    results$summary$`Prob. Direction` =
      sapply(results$summary$`Prob. Direction`,function(x)max(x, 1.0 - x))
    results$CI_level = CI_level
    
    results$data = data_clean
    
    results$bootstrap_samples = 
      t(bootstraps)
    colnames(results$bootstrap_samples) = 
      c(colnames(X),"sigma","alpha")
    
    class(results) = "ladie"
    
    return(results)
    
  }
  if(dose_response == "exponential"){
    
    ## Create negative log posterior function
    nlpost = function(x){
      beta = x[1:P]
      sig = exp(x[P+1])
      
      mu_i = X %*% beta
      if(!nonlinear_time){
        mu_i = mu_i + data_clean$log_time_diff
      }
      
      p_i = 
        sapply(1:nrow(data_clean),
               function(i){
                 helper = function(dummy){
                   (1.0 - exp(-exp(mu_i[i] + dummy * sig))) * dnorm(dummy)
                 }
                 integrate(helper,-4,4)$value
               })
      
      -sum(dbinom(data_clean$y,1,p_i,log = T)) -
        dnorm(beta[1],
              prior_regr_intercept$mean,
              prior_regr_intercept$sd,
              log = TRUE) -
        sum(dnorm(beta[-1],
                  prior_regr_coefs$mean,
                  prior_regr_coefs$sd,
                  log = TRUE)) - 
        dt(sig,df = prior_sigma_df, log = TRUE) - log(sig)
    }
    
    ## Get initial values
    if(verbose) cat("\nGetting initial values...")
    if(isTRUE(nonlinear_time)){
      initial_fit = 
        glm(as.formula(paste0("y ~ log_time_diff + ",
                              paste(varnames$covariates,
                                    collapse = "+"))),
            data = data_clean,
            family = binomial(link = "cloglog"))
    }else{
      initial_fit = 
        glm(as.formula(paste0("y ~ offset(log_time_diff) + ",
                              paste(varnames$covariates,
                                    collapse = "+"))),
            data = data_clean,
            family = binomial(link = "cloglog"))
    }
    inits = coef(initial_fit)
    partial_nlpost = function(x){
      nlpost(c(x[1],inits[-1],x[-1]))
    }
    initial_opt = 
      optim(c(inits[1],0),
            partial_nlpost,
            method = "Nelder-Mead",
            control = list(maxit = 200))
    inits = c(initial_opt$par[1],inits[-1],initial_opt$par[-1])
    
    
    ## Fit full model
    if(verbose) cat("\nPerforming optimization...")
    fit = 
      optim(inits,
            nlpost,
            method = "BFGS")
    
    ## Perform Bayesian bootstrapping
    if(verbose) cat("\nPerforming Bayesian bootstrapping...\n")
    IDs = unique(data_clean[[varnames$id]])
    
    ID_indices = 
      lapply(IDs, function(id) which(data_clean[[varnames$id]] == id))
    
    boot_helper = function(i){
      
      all_levels_represented = FALSE
      while(!all_levels_represented){
        boot_probabilities = 
          drop(extraDistr::rdirichlet(1,rep(1.0,length(IDs))))
        IDs_boot = sample(1:length(IDs),length(IDs),replace=T, prob = boot_probabilities)
        
        data_boot = 
          data_clean[unlist(ID_indices[IDs_boot]),]
        
        n_unique = 
          data_boot |> 
          dplyr::select(all_of(varnames$covariates)) |>
          dplyr::summarize_all(function(x){length(unique(x))})
        
        all_levels_represented = 
          !any(n_unique == 1)
      }
      
      if(isTRUE(nonlinear_time)){
        X_boot = 
          model.matrix(as.formula(paste0("y ~ log_time_diff + ",
                                         paste(varnames$covariates,
                                               collapse = "+"))),
                       data = data_boot)
      }else{
        X_boot = 
          model.matrix(as.formula(paste0("y ~ ",
                                         paste(varnames$covariates,
                                               collapse = "+"))),
                       data = data_boot)
      }
      
      nlpost_boot = function(x){
        beta = x[1:P]
        sig = exp(x[P+1])
        
        mu_i = X_boot %*% beta
        if(!nonlinear_time){
          mu_i = mu_i + data_boot$log_time_diff
        }
        
        p_i = 
          sapply(1:nrow(data_boot),
                 function(i){
                   helper = function(dummy){
                     (1.0 - exp(-exp(mu_i[i] + dummy * sig))) * dnorm(dummy)
                   }
                   integrate(helper,-4,4)$value
                 })
        
        -sum(dbinom(data_boot$y,1,p_i,log = T)) -
          dnorm(beta[1],
                prior_regr_intercept$mean,
                prior_regr_intercept$sd,
                log = TRUE) -
          sum(dnorm(beta[-1],
                    prior_regr_coefs$mean,
                    prior_regr_coefs$sd,
                    log = TRUE)) - 
          dt(sig,df = prior_sigma_df, log = TRUE) - log(sig)
      }
      
      fit_boot = NA
      try({
        fit_boot = 
          optim(fit$par,
                nlpost_boot,
                method = "BFGS")$par
      },silent=T)
      
      return(fit_boot)
    }
    
    if(!missing(cluster)){
      clusterExport(cluster,
                    c("P","IDs","ID_indices","data_clean","nonlinear_time",
                      "varnames","prior_regr_intercept","prior_regr_coefs",
                      "prior_sigma_df","fit"),
                    envir = environment())
      bootstraps = 
        parSapply(cluster,
                  1:n_bootstraps,
                  boot_helper)
    }else{
      bootstraps = 
        sapply(1:n_bootstraps,
               boot_helper)
    }
    
    ## Return results
    results = list()
    
    results$summary = 
      data.frame(Variable = c(colnames(X),"sigma"),
                 `Posterior Mean` = rowMeans(bootstraps),
                 `Posterior Mode` = fit$par,
                 `Posterior SD` = 
                   apply(bootstraps,1,sd),
                 Lower = 
                   apply(bootstraps,1,quantile,probs = (1 - CI_level) / 2),
                 Upper = 
                   apply(bootstraps,1,quantile,probs = 1 - (1 - CI_level) / 2),
                 `Prob. Direction` = 
                   c(apply(bootstraps[1:P,],1,function(x) mean(x > 0)),NA),
                 check.names = FALSE,
                 row.names = NULL)
    results$summary$`Prob. Direction` =
      sapply(results$summary$`Prob. Direction`,function(x)max(x, 1.0 - x))
    
    results$CI_level = CI_level
    
    results$data = data_clean
    
    results$bootstrap_samples = 
      t(bootstraps)
    colnames(results$bootstrap_samples) = 
      c(colnames(X),"sigma")
    
    class(results) = "ladie"
    
    return(results)
    
  }
  
  
}

