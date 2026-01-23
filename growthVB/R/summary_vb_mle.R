#' Summary method for vb_mle objects
#'
#' @param object A fitted model object of class vb_mle
#' @param ... Additional arguments (not used)
#'
#' @return A summary object
#'
#' @export
summary.vb_mle <- function(object, ...) {
  cat("von Bertalanffy Growth Model (Maximum Likelihood Estimation)\n")
  cat("============================================================\n\n")

  # Check if we have multiple models (sex-specific)
  if (is.list(object$model) && !inherits(object$model, "vb_optim")) {
    cat("Sex-specific models fitted\n\n")
    
    # Summary table for all sexes
    all_params <- data.frame()
    all_diagnostics <- data.frame()
    
    for (s in names(object$model)) {
      cat("Model for", s, ":\n")
      cat(rep("-", nchar(paste("Model for", s, ":"))), "\n", sep = "")
      summary_single_model(object$model[[s]])
      cat("\n")
      
      # Collect parameters for comparison table
      params <- object$model[[s]]$parameters
      se <- object$model[[s]]$se
      
      param_row <- data.frame(
        Sex = s,
        Linf = params["Linf"],
        k = params["k"],
        t0 = params["t0"],
        CV = params["cv"],
        row.names = NULL
      )
      
      if (!is.null(se)) {
        param_row$Linf_SE <- se["Linf"]
        param_row$k_SE <- se["k"]
        param_row$t0_SE <- se["t0"]
        param_row$CV_SE <- se["cv"]
      }
      
      all_params <- rbind(all_params, param_row)
      
      # Collect diagnostics
      diag_row <- data.frame(
        Sex = s,
        n_obs = if (!is.null(object$model[[s]]$n_obs)) object$model[[s]]$n_obs else NA,
        logLik = if (!is.null(object$model[[s]]$logLik)) object$model[[s]]$logLik else (-object$model[[s]]$value),
        AIC = if (!is.null(object$model[[s]]$AIC)) object$model[[s]]$AIC else NA,
        BIC = if (!is.null(object$model[[s]]$BIC)) object$model[[s]]$BIC else NA,
        sigma = if (!is.null(object$model[[s]]$sigma)) object$model[[s]]$sigma else NA,
        row.names = NULL
      )
      
      all_diagnostics <- rbind(all_diagnostics, diag_row)
    }
    
    # Print comparison tables
    cat("Parameter Comparison:\n")
    cat("====================\n")
    print(all_params, digits = 4, row.names = FALSE)
    cat("\n")
    
    cat("Model Fit Comparison:\n")
    cat("====================\n")
    print(all_diagnostics, digits = 4, row.names = FALSE)
    
  } else {
    # Single model
    summary_single_model(object$model)
  }

  invisible(object)
}

#' Helper function for summarising a single model
#'
#' @param model A single vb_optim model
#'
#' @keywords internal
summary_single_model <- function(model) {
  # Model fit statistics
  cat("Model Fit:\n")
  if (!is.null(model$n_obs)) {
    cat(sprintf("  Observations: %d\n", model$n_obs))
  }
  if (!is.null(model$n_params)) {
    cat(sprintf("  Parameters: %d\n", model$n_params))
  }
  if (!is.null(model$sigma)) {
    cat(sprintf("  Residual std error: %.4f\n", model$sigma))
  }
  if (!is.null(model$logLik)) {
    cat(sprintf("  Log-likelihood: %.2f\n", model$logLik))
  } else {
    cat(sprintf("  Log-likelihood: %.2f\n", -model$value))
  }
  if (!is.null(model$AIC)) {
    cat(sprintf("  AIC: %.2f\n", model$AIC))
  }
  if (!is.null(model$BIC)) {
    cat(sprintf("  BIC: %.2f\n", model$BIC))
  }
  cat("\n")
  
  # Parameter estimates with detailed statistics
  params <- model$parameters
  se <- model$se
  
  cat("Parameter Estimates:\n")
  if (!is.null(se) && all(!is.na(se))) {
    # Create detailed parameter table
    param_table <- data.frame(
      Parameter = names(params),
      Estimate = params,
      Std.Error = se,
      t.value = params / se,
      p.value = 2 * stats::pnorm(-abs(params / se)),
      row.names = NULL
    )
    
    # Add confidence intervals
    alpha <- 0.05  # 95% CI
    t_crit <- stats::qnorm(1 - alpha/2)
    param_table$CI.lower <- param_table$Estimate - t_crit * param_table$Std.Error
    param_table$CI.upper <- param_table$Estimate + t_crit * param_table$Std.Error
    
    # Add significance stars
    param_table$Signif <- symnum(param_table$p.value,
      corr = FALSE,
      na = FALSE,
      cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
      symbols = c("***", "**", "*", ".", " ")
    )
    
    print(param_table, digits = 4, row.names = FALSE)
    cat("---\n")
    cat("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
    
  } else {
    # Simple parameter table without standard errors
    param_table <- data.frame(
      Parameter = names(params),
      Estimate = params,
      row.names = NULL
    )
    print(param_table, digits = 4, row.names = FALSE)
    cat("(Standard errors not available)\n")
  }
  cat("\n")
  
  # Convergence diagnostics
  cat("Optimisation:\n")
  cat(sprintf("  Method: Maximum Likelihood\n"))
  cat(sprintf("  Convergence: %s\n", 
              if (model$convergence == 0) "Successful" else "Failed"))
  if (model$convergence != 0) {
    cat(sprintf("  Convergence code: %d\n", model$convergence))
  }
  cat(sprintf("  Function evaluations: %d\n", model$counts[1]))
  if (!is.null(model$message)) {
    cat(sprintf("  Message: %s\n", model$message))
  }
  
  # Residual diagnostics if available
  if (!is.null(model$residuals)) {
    cat("\nResidual Summary:\n")
    residual_summary <- summary(model$residuals)
    print(residual_summary)
  }
}
