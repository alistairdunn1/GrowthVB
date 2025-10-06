#' Summarize von Bertalanffy Growth Model Results
#'
#' This function summarizes the results of a von Bertalanffy growth model fit,
#' providing parameter estimates and their confidence/credible intervals.
#'
#' @param model A model object returned by fit_vb_nls() or fit_vb_brms()
#' @param digits Number of decimal places to round to (default 3)
#'
#' @return A data frame with parameter estimates and their intervals
#'
#' @examples
#' \dontrun{
#' # Simple example with simulated data
#' age <- 1:15
#' length <- 100 * (1 - exp(-0.2 * (age - (-0.5)))) + rnorm(15, 0, 5)
#' fit <- fit_vb_nls(age = age, length = length)
#' summarize_vb(fit)
#' }
#'
#' @export
summarize_vb <- function(model, digits = 3) {
  if (inherits(model, "vb_nls")) {
    # Handle nls model summary
    # We compute confidence intervals from the stored nls model(s) when possible
    make_row <- function(nls_model, sex = NULL) {
      cf <- try(stats::coef(nls_model), silent = TRUE)
      # Default NA CIs
      lower <- rep(NA_real_, 3)
      upper <- rep(NA_real_, 3)
      nm <- c("Linf", "k", "t0")
      if (!inherits(cf, "try-error")) {
        ci <- try(stats::confint(nls_model), silent = TRUE)
        if (!inherits(ci, "try-error") && all(nm %in% rownames(ci))) {
          lower <- ci[nm, 1]
          upper <- ci[nm, 2]
        }
        est <- unname(cf[nm])
      } else {
        est <- rep(NA_real_, 3)
      }
      out <- data.frame(
        Parameter = nm,
        Estimate = est,
        Lower_CI = lower,
        Upper_CI = upper,
        stringsAsFactors = FALSE
      )
      if (!is.null(sex)) out$Sex <- sex
      out
    }

    if (is.list(model$model) && !inherits(model$model, "nls")) {
      # Multiple models by sex
      res_list <- list()
      for (s in names(model$model)) {
        res_list[[s]] <- make_row(model$model[[s]], sex = s)
      }
      result <- do.call(rbind, res_list)
      # Ensure column order
      result <- result[, c("Sex", "Parameter", "Estimate", "Lower_CI", "Upper_CI")]
    } else {
      # Single model
      result <- make_row(model$model)
    }

  # Round only numeric columns
  num_cols <- vapply(result, is.numeric, logical(1))
  result[num_cols] <- lapply(result[num_cols], round, digits = digits)
  } else if (inherits(model, "vb_brms")) {
    # Handle brms model summary
    if (is.list(model$models) && !inherits(model$models, "brmsfit")) {
      # Multiple models by sex
      result <- data.frame()
      for (s in names(model$models)) {
        params <- model$parameters[[s]]
        temp <- data.frame(
          Sex = s,
          Parameter = c("Linf", "k", "t0", "tau"),
          Estimate = params[c("b_Linf_Intercept", "b_k_Intercept", "b_t0_Intercept", "b_tau_Intercept"), "Estimate"],
          Lower_CI = params[c("b_Linf_Intercept", "b_k_Intercept", "b_t0_Intercept", "b_tau_Intercept"), "Q2.5"],
          Upper_CI = params[c("b_Linf_Intercept", "b_k_Intercept", "b_t0_Intercept", "b_tau_Intercept"), "Q97.5"]
        )
        result <- rbind(result, temp)
      }
    } else {
      # Single model
      params <- model$parameters
      result <- data.frame(
        Parameter = c("Linf", "k", "t0", "tau"),
        Estimate = params[c("b_Linf_Intercept", "b_k_Intercept", "b_t0_Intercept", "b_tau_Intercept"), "Estimate"],
        Lower_CI = params[c("b_Linf_Intercept", "b_k_Intercept", "b_t0_Intercept", "b_tau_Intercept"), "Q2.5"],
        Upper_CI = params[c("b_Linf_Intercept", "b_k_Intercept", "b_t0_Intercept", "b_tau_Intercept"), "Q97.5"]
      )
    }

  num_cols <- vapply(result, is.numeric, logical(1))
  result[num_cols] <- lapply(result[num_cols], round, digits = digits)
  } else {
    stop("Input must be a model object from fit_vb_nls() or fit_vb_brms()")
  }

  return(result)
}
