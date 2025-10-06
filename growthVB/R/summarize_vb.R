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
    params <- model$parameters

    if (is.null(params$sex)) {
      # Single model
      result <- data.frame(
        Parameter = c("Linf", "k", "t0"),
        Estimate = c(params$Linf, params$k, params$t0),
        Lower_CI = c(params$lowerCI[1], params$lowerCI[2], params$lowerCI[3]),
        Upper_CI = c(params$upperCI[1], params$upperCI[2], params$upperCI[3])
      )
    } else {
      # Multiple models by sex
      result <- data.frame()
      for (s in unique(params$sex)) {
        subset_params <- params[params$sex == s, ]
        temp <- data.frame(
          Sex = s,
          Parameter = c("Linf", "k", "t0"),
          Estimate = c(subset_params$Linf, subset_params$k, subset_params$t0),
          Lower_CI = c(subset_params$lowerCI[1], subset_params$lowerCI[2], subset_params$lowerCI[3]),
          Upper_CI = c(subset_params$upperCI[1], subset_params$upperCI[2], subset_params$upperCI[3])
        )
        result <- rbind(result, temp)
      }
    }

    result <- round(result, digits)
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

    result <- round(result, digits)
  } else {
    stop("Input must be a model object from fit_vb_nls() or fit_vb_brms()")
  }

  return(result)
}
