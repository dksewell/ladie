#' ladie: LAtent Dose Incidence Estimation
#' 
#' 
#' @export 
#' @import dplyr
#' @import parallel
#' @import Matrix



if(FALSE){
  library(dplyr)
  library(parallel)
  library(Matrix)
  
  pathome = 
    read.csv("C:/Users/dksewell/Documents/ladie_helpers/data/micro-nairobi-top10perc.csv")
  fit_linear = 
    ladie(i_salmonella_bin ~ floor + any_flood + exclusively_breastfed + (time | id),
          data = pathome)
  fit_linear_exponential = 
    ladie(i_salmonella_bin ~ floor + any_flood + exclusively_breastfed + (time | id),
          data = pathome,
          dose_response = "exp")
  fit_nonlinear = 
    ladie(i_salmonella_bin ~ floor + any_flood + exclusively_breastfed + (time | id),
          data = pathome,
          nonlinear_time = TRUE)
  
  diag(fit$asymptotic_covariance)
  
  fit$MCMC$samples %>% matplot(type='l',lty=1)
  
  data = 
    read.csv("C:/Users/dksewell/Documents/ladie_helpers/data/micro-nairobi-top10perc.csv")
  
  formula = 
    i_salmonella_bin ~ floor + any_flood + exclusively_breastfed + (time | id)
  dose_response = c("beta","exp")[1]
  prior_regr_coefs = list(mean = 0, sd = 2.5, autoscale = TRUE)
  prior_regr_intercept = list(mean = 0, sd = 100)
  prior_sigma_df = 3
  prior_alpha_rate = 1
  nonlinear_time = FALSE
  CI_level = 0.95
  verbose = TRUE
  # cluster = 10
  
}


ladie = function(formula,
                 data,
                 dose_response = c("beta-poisson","exponential")[1],
                 prior_regr_coefs = list(mean = 0, sd = 2.5, autoscale = TRUE),
                 prior_regr_intercept = list(mean = 0, sd = 10),
                 prior_sigma_df = 3,
                 prior_alpha_rate = 1,
                 nonlinear_time = FALSE,
                 CI_level = 0.95,
                 verbose = TRUE,
                 cluster){
  
  dose_response =
    match.arg(dose_response,
              c("beta-poisson","exponential"))
  
  # Set up parallelization
  if(!missing(cluster)){
    if(is.numeric(cluster)){
      if(verbose) cat("\nSetting up parallel environment")
      cluster = makeCluster(cluster)
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
             gregexpr("\\)",formula_string[[3]])[[1]] - 1) %>% 
      trimws()
    ## Get time variable
    varnames$time = 
      substr(formula_string[[3]],
             gregexpr("\\(",formula_string[[3]])[[1]] + 1,
             gregexpr("\\|",formula_string[[3]])[[1]] - 1) %>% 
      trimws()
    ## Get covariates
    plus_locations = 
      c(0,gregexpr("\\+",formula_string[[3]])[[1]])
    varnames$covariates = 
      sapply(1:(length(plus_locations) - 1),
             function(i){
               substr(formula_string[[3]],
                      plus_locations[i] + 1,
                      plus_locations[i + 1] - 1) %>% 
                 trimws()
             })
    rm(plus_locations,formula_string)
  }
  
  
  # Get time variable shifted to start at zero
  {
    data_clean =
      data_clean |>
      dplyr::group_by(id) %>% 
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
      dplyr::group_by(id) %>% 
      dplyr::mutate(S_it = helper_function(y)) %>% 
      dplyr::ungroup() %>% 
      dplyr::filter(time_diff > 0,
                    S_it == 1) %>% 
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
    if(missing(cluster)){
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
                   1 - integrate(helper,-4,4)$value
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
    }else{
      clusterExport(cluster,
                    c("X","data_clean","nonlinear_time","P"),
                    envir=environment())
      nlpost = function(x){
        beta = x[1:P]
        sig = exp(x[P+1])
        a = exp(x[P+2])
        
        mu_i = X %*% beta
        if(!nonlinear_time){
          mu_i = mu_i + data_clean$log_time_diff
        }
        
        p_i = 
          parSapply(cluster,
                    1:nrow(data_clean),
                    function(i){
                      helper = function(dummy){
                        (1.0 - (1.0 + exp(mu_i[i] + dummy * sig))^(-a)) * dnorm(dummy)
                      }
                      1 - integrate(helper,-4,4)$value
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
      nlpost(c(inits,x))
    }
    initial_opt = 
      optim(c(0,0),
            partial_nlpost,
            method = "BFGS")
    inits = c(inits,initial_opt$par)
    
    
    ## Fit full model
    if(verbose) cat("\nPerforming optimization...")
    fit = 
      optim(inits,
            nlpost,
            method = "BFGS",
            hessian = TRUE)
    
    post_means = fit$par
    post_means[P + 1:2] = 
      exp(post_means[P + 1:2])
    names(post_means)[P + 1:2] = c("sigma","alpha")
    
    Sigma = qr.solve(fit$hessian)
    
    ## Check for numerical stability
    Sigma = 0.5 * (Sigma + t(Sigma))
    if(!matrixcalc::is.positive.definite(Sigma)){
      warning("Asymptotic covariance not positive definite.  Don't trust results.",
              immediate. = TRUE)
      Sigma = as.matrix(Matrix::nearPD(Sigma)$mat)
    }
    
    ## Get SD and CIs on original scale
    theta_draws = 
      cbind(rnorm(1e4,
                  fit$par[P+1],
                  sqrt(Sigma[P+1,P+1])),
            rnorm(1e4,
                  fit$par[P+2],
                  sqrt(Sigma[P+2,P+2]))
      )
    theta_draws = exp(theta_draws)
    
    
    ## Return results
    results = list()
    
    results$summary = 
      data.frame(Variable = names(post_means),
                 `Posterior Median` = post_means,
                 `Posterior SD` = 
                   sqrt(c(diag(Sigma)[1:P],
                          apply(theta_draws,2,var))),
                 Lower = 
                   c(qnorm((1-CI_level)/2,
                           post_means[1:P],
                           sqrt(diag(Sigma)[1:P])),
                     apply(theta_draws,2,quantile,probs = (1-CI_level)/2)),
                 Upper = 
                   c(qnorm(1 - (1-CI_level)/2,
                           post_means[1:P],
                           sqrt(diag(Sigma)[1:P])),
                     apply(theta_draws,2,quantile,probs = 1 - (1-CI_level)/2)),
                 `Prob. Direction` = 
                   c(pnorm(0,
                           post_means[1:P],
                           sqrt(diag(Sigma)[1:P])),
                     rep(NA,2)),
                 check.names = FALSE,
                 row.names = NULL)
    results$summary$`Prob. Direction` = 
      ifelse(results$summary$`Prob. Direction` < 0.5,
             1 - results$summary$`Prob. Direction`,
             results$summary$`Prob. Direction`)
    
    results$CI_level = CI_level
    
    results$log_posterior = 
      -fit$value
    
    results$asymptotic_covariance = Sigma
    
    class(results) = "ladie"
    
    return(results)
    
  }else{ # Start work for exponential DR model
    
    ## Create negative log posterior function
    if(missing(cluster)){
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
                   1 - integrate(helper,-4,4)$value
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
    }else{
      clusterExport(cluster,
                    c("X","data_clean","nonlinear_time","P"),
                    envir=environment())
      nlpost = function(x){
        beta = x[1:P]
        sig = exp(x[P+1])
        
        mu_i = X %*% beta
        if(!nonlinear_time){
          mu_i = mu_i + data_clean$log_time_diff
        }
        
        p_i = 
          parSapply(cluster,
                    1:nrow(data_clean),
                    function(i){
                      helper = function(dummy){
                        (1.0 - exp(-exp(mu_i[i] + dummy * sig))) * dnorm(dummy)
                      }
                      1 - integrate(helper,-4,4)$value
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
      nlpost(c(inits,x))
    }
    initial_opt = 
      optimize(partial_nlpost,interval = c(-100,100))
    inits = c(inits,initial_opt$minimum)
    
    
    ## Fit full model
    if(verbose) cat("\nPerforming optimization...")
    fit = 
      optim(inits,
            nlpost,
            method = "BFGS",
            hessian = TRUE)
    
    post_means = fit$par
    post_means[P + 1] = 
      exp(post_means[P + 1])
    names(post_means)[P + 1] = c("sigma")
    
    Sigma = qr.solve(fit$hessian)
    
    ## Check for numerical stability
    Sigma = 0.5 * (Sigma + t(Sigma))
    if(!matrixcalc::is.positive.definite(Sigma)){
      warning("Asymptotic covariance not positive definite.  Don't trust results.",
              immediate. = TRUE)
      Sigma = as.matrix(Matrix::nearPD(Sigma)$mat)
    }
    
    ## Get SD and CIs on original scale
    theta_draws = 
      rnorm(1e4,
            fit$par[P+1],
            sqrt(Sigma[P+1,P+1]))
    theta_draws = exp(theta_draws)
    
    
    ## Return results
    results = list()
    
    results$summary = 
      data.frame(Variable = names(post_means),
                 `Posterior Median` = post_means,
                 `Posterior SD` = 
                   sqrt(c(diag(Sigma)[1:P],
                          var(theta_draws))),
                 Lower = 
                   c(qnorm((1-CI_level)/2,
                           post_means[1:P],
                           sqrt(diag(Sigma)[1:P])),
                     quantile(theta_draws,probs = (1-CI_level)/2)),
                 Upper = 
                   c(qnorm(1 - (1-CI_level)/2,
                           post_means[1:P],
                           sqrt(diag(Sigma)[1:P])),
                     quantile(theta_draws,probs = 1 - (1-CI_level)/2)),
                 `Prob. Direction` = 
                   c(pnorm(0,
                           post_means[1:P],
                           sqrt(diag(Sigma)[1:P])),
                     NA),
                 check.names = FALSE,
                 row.names = NULL)
    results$summary$`Prob. Direction` = 
      ifelse(results$summary$`Prob. Direction` < 0.5,
             1 - results$summary$`Prob. Direction`,
             results$summary$`Prob. Direction`)
    
    results$CI_level = CI_level
    
    results$log_posterior = 
      -fit$value
    
    results$asymptotic_covariance = Sigma
    
    class(results) = "ladie"
    
    return(results)
    
  }
  
}

