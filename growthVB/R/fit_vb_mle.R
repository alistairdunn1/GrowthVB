#' Fit von Bertalanffy Growth Function using Maximum Likelihood Estimation
#'
#' This function fits a von Bertalanffy growth function to age and length data using
#' maximum likelihood estimation. The CV of the length is explicitly modelled as a function
#' of the mean length, allowing for heteroscedasticity where variance increases with fish size.
#'
#' @param age A numeric vector of ages
#' @param length A numeric vector of lengths
#' @param sex An optional factor or character vector specifying the sex for each observation
#' @param length_bins Optional numeric vector specifying length bins for analysis
#' @param sampling_prob Optional vector of sampling probabilities (defaults to 1)
#' @param ci_level Confidence interval level (default 0.95)
#' @param optim_method Optimisation method for stats::optim (default "L-BFGS-B")
#' @param maxit Maximum number of iterations for optimisation (default 1000)
#'
#' @return A list containing:
#'   \item{parameters}{Estimated parameters with confidence intervals}
#'   \item{data}{Original data with fitted values and residuals}
#'   \item{fits}{Model fit predictions for plotting}
#'   \item{model}{The fitted model object with optimisation results, including:
#'     \itemize{
#'       \item Log-likelihood of the fitted model (logLik)
#'       \item Akaike Information Criterion (AIC)
#'       \item Bayesian Information Criterion (BIC)
#'       \item Residual standard error (sigma)
#'       \item Number of observations (n_obs)
#'       \item Number of parameters (n_params)
#'       \item Model predictions (fitted_values)
#'       \item Model residuals (residuals)
#'     }
#'   }
#'
#' @examples
#' \dontrun{
#' # Simple example with simulated data
#' age <- 1:15
#' length <- 100 * (1 - exp(-0.2 * (age - (-0.5)))) + rnorm(15, 0, 5 + 0.1 * (1:15))
#' fit <- fit_vb_mle(age = age, length = length)
#' }
#'
#' @export
fit_vb_mle <- function(age, length, sex = NULL, length_bins = NULL,
                       sampling_prob = 1, ci_level = 0.95,
                       optim_method = "L-BFGS-B", maxit = 1000) {
  # Create data frame for analysis
  data <- data.frame(age = age, length = length)

  # Add sex if provided
  if (!is.null(sex)) {
    data$sex <- sex
    # Filter out any NA values
    data <- data[!is.na(data$age) & !is.na(data$length) & !is.na(data$sex), ]
  } else {
    # Filter out any NA values
    data <- data[!is.na(data$age) & !is.na(data$length), ]
  }

  # Von Bertalanffy growth function
  vb_mean <- function(age, Linf, k, t0) {
    Linf * (1 - exp(-k * (age - t0)))
  }

  # Function to fit a model to a subset of the data using maximum likelihood
  fit_model <- function(subset_data) {
    # Better initial parameter estimates based on data
    Linf_init <- max(subset_data$length, na.rm = TRUE) * 1.1 # Slightly above max observed length

    # Estimate k from growth rate - use simple linear approximation
    if (nrow(subset_data) > 1) {
      # Simple linear model to get rough k estimate
      lm_rough <- lm(log(Linf_init - subset_data$length + 1) ~ subset_data$age)
      k_init <- max(0.05, min(2.0, -coef(lm_rough)[2])) # Bound k between 0.05 and 2.0
    } else {
      k_init <- 0.2
    }

    # Estimate t0 from intercept or use reasonable default
    t0_init <- min(subset_data$age) - 1 # Set t0 to be before the youngest age

    cv_init <- 0.15 # Initial coefficient of variation

    # Negative log-likelihood function with CV as a function of predicted length
    neg_log_likelihood <- function(params) {
      # Extract parameters
      Linf <- params[1]
      k <- params[2]
      t0 <- params[3]
      cv <- params[4]

      # Check for invalid parameter values
      if (Linf <= 0 || k <= 0 || cv <= 0) {
        return(1e10) # Return large value for invalid parameters
      }

      # Calculate predicted lengths
      pred_lengths <- vb_mean(subset_data$age, Linf, k, t0)

      # Calculate standard deviations as CV * predicted length
      # Add a small constant to avoid zeros that can cause NaNs
      sds <- cv * pred_lengths

      # Ensure standard deviations are valid (minimum value of 0.001)
      sds <- pmax(sds, 0.001)

      # Calculate log-likelihoods for each observation
      log_likes <- rep(NA_real_, length(pred_lengths))
      for (i in seq_along(pred_lengths)) {
        # Safe calculation of log-likelihood to avoid NaNs
        x <- subset_data$length[i]
        mu <- pred_lengths[i]
        sigma <- sds[i]

        # If sigma is valid, calculate log-likelihood properly
        if (is.finite(sigma) && sigma > 0) {
          log_likes[i] <- stats::dnorm(x = x, mean = mu, sd = sigma, log = TRUE)
        } else {
          # If sigma is invalid, use a penalized value
          log_likes[i] <- -1e5
        }
      }

      # Compute negative log-likelihood (sum of individual log-likelihoods)
      neg_ll <- -sum(log_likes, na.rm = TRUE)

      # If result is NA or NaN, return a large value
      if (!is.finite(neg_ll)) {
        neg_ll <- 1e10
      }

      # Return negative log-likelihood (for minimisation)
      return(neg_ll)
    }

    # Initial parameter vector
    init_params <- c(Linf = Linf_init, k = k_init, t0 = t0_init, cv = cv_init)

    # Try optimisation with bounds to ensure reasonable parameter values
    if (optim_method %in% c("L-BFGS-B")) {
      # Use bounded optimisation
      lower_bounds <- c(Linf = max(subset_data$length) * 0.5, k = 0.001, t0 = min(subset_data$age) - 5, cv = 0.01)
      upper_bounds <- c(Linf = max(subset_data$length) * 3, k = 3.0, t0 = max(subset_data$age) + 2, cv = 1.0)

      optim_result <- try(
        stats::optim(
          par = init_params,
          fn = neg_log_likelihood,
          method = optim_method,
          lower = lower_bounds,
          upper = upper_bounds,
          hessian = TRUE,
          control = list(maxit = maxit, trace = FALSE)
        ),
        silent = TRUE
      )
    } else {
      # Use unbounded optimisation
      optim_result <- try(
        stats::optim(
          par = init_params,
          fn = neg_log_likelihood,
          method = optim_method,
          hessian = TRUE,
          control = list(maxit = maxit, trace = FALSE)
        ),
        silent = TRUE
      )
    }

    # If first attempt fails, try with Nelder-Mead method
    if (inherits(optim_result, "try-error") || optim_result$convergence != 0) {
      warning("First optimisation attempt failed. Trying with Nelder-Mead method.")
      optim_result <- try(
        stats::optim(
          par = init_params,
          fn = neg_log_likelihood,
          method = "Nelder-Mead",
          hessian = TRUE,
          control = list(maxit = 2000, trace = FALSE)
        ),
        silent = TRUE
      )

      # If still failing, try with different starting values
      if (inherits(optim_result, "try-error") || optim_result$convergence != 0) {
        warning("Second optimisation attempt failed. Trying with different starting values.")
        init_params <- c(
          Linf = max(subset_data$length, na.rm = TRUE) * 1.5,
          k = 0.15,
          t0 = 0,
          cv = 0.2
        )

        optim_result <- try(
          stats::optim(
            par = init_params,
            fn = neg_log_likelihood,
            method = "BFGS",
            hessian = TRUE,
            control = list(maxit = 2000, trace = FALSE)
          ),
          silent = TRUE
        )

        if (inherits(optim_result, "try-error") || optim_result$convergence != 0) {
          stop("All optimisation attempts failed. Try different starting values or check data quality.")
        }
      }
    }

    # Create a model object with the optimisation results and other useful information
    model <- list(
      parameters = optim_result$par,
      convergence = optim_result$convergence,
      value = optim_result$value,
      counts = optim_result$counts,
      message = optim_result$message,
      hessian = optim_result$hessian,
      data = subset_data,
      vb_mean = vb_mean
    )
    
    # Calculate model diagnostics
    n_obs <- nrow(subset_data)
    n_params <- length(optim_result$par)
    
    # Calculate log-likelihood (negative of the minimised objective)
    logLik <- -optim_result$value
    
    # Calculate information criteria
    AIC <- -2 * logLik + 2 * n_params
    BIC <- -2 * logLik + log(n_obs) * n_params
    
    # Calculate fitted values and residuals for additional diagnostics
    fitted_values <- vb_mean(subset_data$age, 
                            optim_result$par["Linf"], 
                            optim_result$par["k"], 
                            optim_result$par["t0"])
    residuals <- subset_data$length - fitted_values
    
    # Calculate residual standard error
    sigma <- sqrt(sum(residuals^2) / (n_obs - n_params))
    
    # Add diagnostics to model object
    model$logLik <- logLik
    model$AIC <- AIC
    model$BIC <- BIC
    model$sigma <- sigma
    model$n_obs <- n_obs
    model$n_params <- n_params
    model$fitted_values <- fitted_values
    model$residuals <- residuals

    # Calculate parameter standard errors if hessian is available
    if (!is.null(model$hessian)) {
      tryCatch(
        {
          # Check if Hessian is positive definite
          eigenvalues <- eigen(model$hessian, only.values = TRUE)$values
          if (any(eigenvalues <= 0)) {
            warning("Hessian matrix is not positive definite. Standard errors may be unreliable.")
          }

          vcov_matrix <- solve(model$hessian)
          model$se <- sqrt(diag(vcov_matrix))
          model$vcov <- vcov_matrix
        },
        error = function(e) {
          warning("Could not compute standard errors from Hessian matrix: ", e$message)
          warning("This often indicates convergence issues or poorly conditioned optimisation.")
          model$se <- rep(NA, length(model$parameters))
          model$vcov <- matrix(NA, nrow = length(model$parameters), ncol = length(model$parameters))
        }
      )
    } else {
      warning("Hessian matrix not available. Cannot compute confidence intervals.")
      model$se <- rep(NA, length(model$parameters))
      model$vcov <- matrix(NA, nrow = length(model$parameters), ncol = length(model$parameters))
    }

    class(model) <- c("vb_optim", "list")
    return(model)
  }

  # Results list
  results <- list()

  # Split data by sex if provided, otherwise fit single model
  if (!is.null(sex) && length(unique(sex[!is.na(sex)])) > 1) {
    # Process each sex separately
    sex_levels <- unique(sex[!is.na(sex)])
    params_list <- list()
    data_list <- list()
    fits_list <- list()
    models_list <- list()

    for (s in sex_levels) {
      subset_data <- data[data$sex == s, ]
      model <- fit_model(subset_data)
      models_list[[s]] <- model

      # Extract parameter estimates and create a row for each parameter with sex prefix
      # The tests expect parameters in sex-specific models to be named like "M_Linf", "F_k", etc.
      param_names <- c("Linf", "k", "t0", "cv")
      for (p in param_names) {
        param_name <- paste0(s, "_", p)
        estimate <- model$parameters[p]
        std_error <- if (!is.null(model$se)) model$se[p] else NA

        params_row <- c(estimate = estimate, std.error = std_error)
        params_list[[param_name]] <- params_row
      }

      # Also add the original parameter structure to the model for internal use
      model$orig_parameters <- model$parameters

      # Calculate fitted values and residuals
      fitted_vals <- model$vb_mean(
        subset_data$age,
        model$parameters["Linf"],
        model$parameters["k"],
        model$parameters["t0"]
      )
      residuals <- subset_data$length - fitted_vals

      # Calculate standardized residuals using model CV
      sds <- model$parameters["cv"] * fitted_vals
      # Ensure no zero or NA standard deviations
      sds <- pmax(sds, 0.001)
      std_residuals <- residuals / sds

      # Ensure no NaN or Inf values in residuals
      std_residuals[!is.finite(std_residuals)] <- NA

      subset_data$fitted <- fitted_vals
      subset_data$residual <- residuals
      subset_data$student <- std_residuals

      data_list[[s]] <- subset_data

      # Create data for plotting smooth curves
      age_seq <- seq(0, max(subset_data$age, na.rm = TRUE), length.out = 100)
      new_data <- data.frame(age = age_seq)

      # Predict from model
      pred_mean <- model$vb_mean(
        new_data$age,
        model$parameters["Linf"],
        model$parameters["k"],
        model$parameters["t0"]
      )

      # Calculate confidence intervals if standard errors are available
      if (!is.null(model$se) && !is.null(model$vcov)) {
        # Simplified calculation of confidence intervals based on delta method
        # Note: This is an approximation
        t_value <- stats::qt(1 - (1 - ci_level) / 2, df = nrow(subset_data) - length(model$parameters))

        # Get CI width for each age point (simplified approach)
        se_pred <- rep(NA_real_, length(age_seq))
        for (i in seq_along(age_seq)) {
          # Gradient of the VB function with respect to parameters
          grad <- c(
            1 - exp(-model$parameters["k"] * (age_seq[i] - model$parameters["t0"])), # dL/dLinf
            model$parameters["Linf"] * exp(-model$parameters["k"] * (age_seq[i] - model$parameters["t0"])) * (age_seq[i] - model$parameters["t0"]), # dL/dk
            model$parameters["Linf"] * exp(-model$parameters["k"] * (age_seq[i] - model$parameters["t0"])) * model$parameters["k"] # dL/dt0
          )

          # Subset vcov matrix for growth parameters (exclude cv)
          vcov_growth <- model$vcov[1:3, 1:3]

          # Calculate standard error of prediction
          se_pred[i] <- sqrt(t(grad) %*% vcov_growth %*% grad)
        }

        lowerCI <- pred_mean - t_value * se_pred
        upperCI <- pred_mean + t_value * se_pred
      } else {
        se_pred <- rep(NA_real_, length(age_seq))
        lowerCI <- rep(NA_real_, length(age_seq))
        upperCI <- rep(NA_real_, length(age_seq))
      }

      fit_data <- data.frame(
        age = age_seq,
        mean = pred_mean,
        se = se_pred,
        lowerCI = lowerCI,
        upperCI = upperCI,
        sex = s
      )

      fits_list[[s]] <- fit_data
    }

    # Combine results
    # Create parameter matrix from the list, with appropriate structure for tests
    param_matrix <- matrix(
      unlist(params_list),
      nrow = length(params_list),
      ncol = 2,
      byrow = TRUE,
      dimnames = list(names(params_list), c("estimate", "std.error"))
    )

    results$parameters <- param_matrix
    results$data <- do.call(rbind, data_list)
    results$fits <- do.call(rbind, fits_list)
    results$model <- models_list
    results$models <- models_list # Add models alias for backward compatibility
  } else {
    # Fit a single model to all data
    model <- fit_model(data)

    # Create a parameter matrix for compatibility with tests
    params <- matrix(nrow = 4, ncol = 2)
    rownames(params) <- c("Linf", "k", "t0", "cv")
    colnames(params) <- c("estimate", "std.error")

    # Fill with parameter estimates
    params["Linf", "estimate"] <- model$parameters["Linf"]
    params["k", "estimate"] <- model$parameters["k"]
    params["t0", "estimate"] <- model$parameters["t0"]
    params["cv", "estimate"] <- model$parameters["cv"]

    # Add standard errors if available
    if (!is.null(model$se)) {
      params["Linf", "std.error"] <- model$se["Linf"]
      params["k", "std.error"] <- model$se["k"]
      params["t0", "std.error"] <- model$se["t0"]
      params["cv", "std.error"] <- model$se["cv"]
    } else {
      params[, "std.error"] <- NA
    }

    # Calculate fitted values and residuals
    fitted_vals <- model$vb_mean(
      data$age,
      model$parameters["Linf"],
      model$parameters["k"],
      model$parameters["t0"]
    )
    residuals <- data$length - fitted_vals

    # Calculate standardized residuals using model CV
    sds <- model$parameters["cv"] * fitted_vals
    # Ensure no zero or NA standard deviations
    sds <- pmax(sds, 0.001)
    std_residuals <- residuals / sds

    # Ensure no NaN or Inf values in residuals
    std_residuals[!is.finite(std_residuals)] <- NA

    data$fitted <- fitted_vals
    data$residual <- residuals
    data$student <- std_residuals

    # Create data for plotting smooth curves
    age_seq <- seq(0, max(data$age, na.rm = TRUE), length.out = 100)
    new_data <- data.frame(age = age_seq)

    # Predict from model
    pred_mean <- model$vb_mean(
      new_data$age,
      model$parameters["Linf"],
      model$parameters["k"],
      model$parameters["t0"]
    )

    # Calculate confidence intervals if standard errors are available
    if (!is.null(model$se) && !is.null(model$vcov)) {
      # Simplified calculation of confidence intervals based on delta method
      t_value <- stats::qt(1 - (1 - ci_level) / 2, df = nrow(data) - length(model$parameters))

      # Get CI width for each age point (simplified approach)
      se_pred <- rep(NA_real_, length(age_seq))
      for (i in seq_along(age_seq)) {
        # Gradient of the VB function with respect to parameters
        grad <- c(
          1 - exp(-model$parameters["k"] * (age_seq[i] - model$parameters["t0"])), # dL/dLinf
          model$parameters["Linf"] * exp(-model$parameters["k"] * (age_seq[i] - model$parameters["t0"])) * (age_seq[i] - model$parameters["t0"]), # dL/dk
          model$parameters["Linf"] * exp(-model$parameters["k"] * (age_seq[i] - model$parameters["t0"])) * model$parameters["k"] # dL/dt0
        )

        # Subset vcov matrix for growth parameters (exclude cv)
        vcov_growth <- model$vcov[1:3, 1:3]

        # Calculate standard error of prediction
        se_pred[i] <- sqrt(t(grad) %*% vcov_growth %*% grad)
      }

      lowerCI <- pred_mean - t_value * se_pred
      upperCI <- pred_mean + t_value * se_pred
    } else {
      se_pred <- rep(NA_real_, length(age_seq))
      lowerCI <- rep(NA_real_, length(age_seq))
      upperCI <- rep(NA_real_, length(age_seq))
    }

    fit_data <- data.frame(
      age = age_seq,
      mean = pred_mean,
      se = se_pred,
      lowerCI = lowerCI,
      upperCI = upperCI
    )

    # Store results
    results$parameters <- params
    results$data <- data
    results$fits <- fit_data
    results$model <- model
  }

  class(results) <- c("vb_mle", "list")
  return(results)
}
