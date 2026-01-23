#' Summarise von Bertalanffy Growth Model Results
#'
#' This function summarises the results of a von Bertalanffy growth model fit,
#' providing parameter estimates and their confidence/credible intervals.
#'
#' @param model A model object returned by fit_vb_mle() or fit_vb_brms()
#' @param digits Number of decimal places to round to (default 3)
#'
#' @return A data frame with parameter estimates and their intervals
#'
#' @examples
#' \dontrun{
#' # Simple example with simulated data
#' age <- 1:15
#' length <- 100 * (1 - exp(-0.2 * (age - (-0.5)))) + rnorm(15, 0, 5)
#' fit <- fit_vb_mle(age = age, length = length)
#' summarise_vb(fit)
#' }
#'
#' @export
summarise_vb <- function(model, digits = 3) {
  if (inherits(model, "vb_mle")) {
    # Handle MLE model summary
    if ((!is.null(model$models) && is.list(model$models)) ||
      (!is.null(model$model) && is.list(model$model) && !inherits(model$model, "vb_optim"))) {
      # Sex-specific models
      sex_models <- if (!is.null(model$models)) model$models else model$model
      result <- data.frame(
        Sex = character(),
        Parameter = character(),
        Estimate = numeric(),
        Lower = numeric(),
        Upper = numeric(),
        stringsAsFactors = FALSE
      )

      for (s in names(sex_models)) {
        for (p in c("Linf", "k", "t0", "cv")) {
          param_name <- paste0(s, "_", p)
          if (param_name %in% rownames(model$parameters)) {
            est <- model$parameters[param_name, "estimate"]
            se <- model$parameters[param_name, "std.error"]
            ci_delta <- ifelse(!is.na(se), 1.96 * se, NA_real_)

            result <- rbind(result, data.frame(
              Sex = s,
              Parameter = p,
              Estimate = round(est, digits),
              Lower = round(est - ci_delta, digits),
              Upper = round(est + ci_delta, digits),
              stringsAsFactors = FALSE
            ))
          }
        }
      }

      return(result)
    } else {
      params <- model$parameters

      if (is.matrix(params)) {
        est <- params[, "estimate"]
        se <- params[, "std.error"]
        ci_delta <- ifelse(!is.na(se), 1.96 * se, NA_real_)

        result <- data.frame(
          Parameter = rownames(params),
          Estimate = round(est, digits),
          Lower = round(est - ci_delta, digits),
          Upper = round(est + ci_delta, digits),
          stringsAsFactors = FALSE
        )
        return(result)
      }

      if (is.data.frame(params)) {
        result <- data.frame(
          Parameter = colnames(params),
          Estimate = as.numeric(params[1, ]),
          Lower = NA,
          Upper = NA,
          stringsAsFactors = FALSE
        )
        num_cols <- vapply(result, is.numeric, logical(1))
        result[num_cols] <- lapply(result[num_cols], round, digits = digits)
        return(result)
      }

      if (!is.null(model$model) && is.numeric(model$model$parameters)) {
        params_vec <- model$model$parameters
        result <- data.frame(
          Parameter = names(params_vec),
          Estimate = round(as.numeric(params_vec), digits),
          Lower = NA,
          Upper = NA,
          stringsAsFactors = FALSE
        )
        return(result)
      }

      stop("Unsupported vb_mle parameter structure")
    }
  } else if (inherits(model, "vb_brms")) {
    # Handle brms model summary
    if (is.list(model$models) && !inherits(model$models, "brmsfit")) {
      # Multiple models by sex
      result <- data.frame()
      for (s in names(model$models)) {
        params <- model$parameters[[s]]
        temp <- data.frame(
          Sex = s,
          Parameter = c("Linf", "k", "t0", "CV"),
          Estimate = params[c("Linf_Intercept", "k_Intercept", "t0_Intercept", "CV_Intercept"), "Estimate"],
          Lower = params[c("Linf_Intercept", "k_Intercept", "t0_Intercept", "CV_Intercept"), "Q2.5"],
          Upper = params[c("Linf_Intercept", "k_Intercept", "t0_Intercept", "CV_Intercept"), "Q97.5"]
        )
        result <- rbind(result, temp)
      }
    } else {
      # Single model
      params <- model$parameters
      result <- data.frame(
        Parameter = c("Linf", "k", "t0", "CV"),
        Estimate = params[c("Linf_Intercept", "k_Intercept", "t0_Intercept", "CV_Intercept"), "Estimate"],
        Lower = params[c("Linf_Intercept", "k_Intercept", "t0_Intercept", "CV_Intercept"), "Q2.5"],
        Upper = params[c("Linf_Intercept", "k_Intercept", "t0_Intercept", "CV_Intercept"), "Q97.5"]
      )
    }

    num_cols <- vapply(result, is.numeric, logical(1))
    result[num_cols] <- lapply(result[num_cols], round, digits = digits)

    return(result)
  } else {
    stop("Input must be a model object from fit_vb_mle() or fit_vb_brms()")
  }
}
