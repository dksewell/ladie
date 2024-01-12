#' ladie: LAtent Dose Incidence Estimation
#' 
#' Use longitudinal data to estimate risk ratios of accruing dose
#' 
#' @param formula an object of class "formula" (or one that can be coerced to 
#' that class): a symbolic description of the model to be fitted. 
#' @param data a data.frame containing the variables in the model.
#' @param dose_response The dose-response model to be used.  Either the 
#' beta-poisson, exponential.  The simple threshold model is not finished yet.
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
#' @param verbose logical. Should progress be displayed?
#' @param optim_args named list giving arguments to optim().
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
#' @importFrom matrixcalc is.positive.definite

ladie = function(formula,
                 data,
                 dose_response = c("beta-poisson","exponential","simple_threshold")[1],
                 prior_regr_coefs = list(mean = 0, sd = 2.5, autoscale = TRUE),
                 prior_regr_intercept = list(mean = 0, sd = 10),
                 prior_sigma_df = 3,
                 prior_alpha_rate = 1,
                 nonlinear_time = FALSE,
                 CI_level = 0.95,
                 verbose = TRUE,
                 optim_args = list(method = "BFGS")){
  
  dose_response =
    match.arg(dose_response,
              c("beta-poisson","exponential","simple_threshold"))
  
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
    sd_x = apply(X[,-1,drop=F],2,sd)
    names(sd_x) = colnames(X)[-1]
    prior_regr_coefs$sd = 
      prior_regr_coefs$sd / sd_x
    rm(sd_x)
  }
  
  # Create object to store warnings/errors
  errors = NULL
  
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
    optim_args$par = inits
    optim_args$fn = nlpost
    optim_args$hessian = TRUE
    if(is.null(optim_args$method)) optim_args$method = "BFGS"
    fit = 
      do.call(optim,
              optim_args)
    
    if(fit$conv != 0) errors = "BFGS algo did not converge"
    
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
      
      errors = c(errors,"Estimated covariance matrix not PSD")
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
    
    results$log_likelihood = 
      results$log_posterior - 
      dnorm(results$summary$`Posterior Median`[1],
            prior_regr_intercept$mean,
            prior_regr_intercept$sd,
            log = TRUE) -
      sum(dnorm(results$summary$`Posterior Median`[2:P],
                prior_regr_coefs$mean,
                prior_regr_coefs$sd,
                log = TRUE)) - 
      dt(results$summary$`Posterior Median`[P+1],
         df = prior_sigma_df,
         log = TRUE) - log(results$summary$`Posterior Median`[P+1]) -
      dexp(results$summary$`Posterior Median`[P+2],
           rate = prior_alpha_rate,
           log = TRUE) - log(results$summary$`Posterior Median`[P+2])
    
    results$data = data_clean
    
    results$asymptotic_covariance = Sigma
    
    results$prior_regr_coefs = prior_regr_coefs
    results$prior_regr_intercept = prior_regr_intercept
    results$prior_sigma_df = prior_sigma_df
    results$prior_alpha_rate = prior_alpha_rate
    
    results$errors = errors
    
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
    optim_args$par = inits
    optim_args$fn = nlpost
    optim_args$hessian = TRUE
    if(is.null(optim_args$method)) optim_args$method = "BFGS"
    fit = 
      do.call(optim,
              optim_args)
    
    
    if(fit$conv != 0) errors = "BFGS algo did not converge"
    
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
      
      errors = c(errors,"Estimated covariance matrix not PSD")
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
    
    results$log_likelihood = 
      results$log_posterior - 
      dnorm(results$summary$`Posterior Median`[1],
            prior_regr_intercept$mean,
            prior_regr_intercept$sd,
            log = TRUE) -
      sum(dnorm(results$summary$`Posterior Median`[2:P],
                prior_regr_coefs$mean,
                prior_regr_coefs$sd,
                log = TRUE)) - 
      dt(results$summary$`Posterior Median`[P+1],
         df = prior_sigma_df,
         log = TRUE) - log(results$summary$`Posterior Median`[P+1])
    
    results$data = data_clean
    
    results$asymptotic_covariance = Sigma
    
    results$prior_regr_coefs = prior_regr_coefs
    results$prior_regr_intercept = prior_regr_intercept
    results$prior_sigma_df = prior_sigma_df
    
    results$errors = errors
    
    class(results) = "ladie"
    
    return(results)
    
  }
  
  
}
