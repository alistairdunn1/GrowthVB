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
      # Use models if available, otherwise use model
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
        # Extract parameters for this sex
        for (p in c("Linf", "k", "t0", "cv")) {
          # Get parameter name with sex prefix
          param_name <- paste0(s, "_", p)
          if (param_name %in% rownames(model$parameters)) {
            est <- model$parameters[param_name, "estimate"]
            se <- model$parameters[param_name, "std.error"]
            lower <- est - 1.96 * se
            upper <- est + 1.96 * se

            result <- rbind(result, data.frame(
              Sex = s,
              Parameter = p,
              Estimate = round(est, digits),
              Lower = round(lower, digits),
              Upper = round(upper, digits),
              stringsAsFactors = FALSE
            ))
          }
        }
      }

      return(result)
    } else {
      # Single model - new MLE implementation
      if (is.data.frame(model$parameters)) {
        # Handle data.frame parameters format
        result <- data.frame(
          Parameter = colnames(model$parameters),
          Estimate = as.numeric(model$parameters[1, ]),
          Lower = NA,
          Upper = NA,
          stringsAsFactors = FALSE
        )

        # Round numeric values
        num_cols <- vapply(result, is.numeric, logical(1))
        result[num_cols] <- lapply(result[num_cols], round, digits = digits)

        return(result)
      } else if (is.vector(model$model$parameters) && !is.null(names(model$model$parameters))) {
        # Handle case where parameters are in model$model$parameters as a named vector
        params <- model$model$parameters
        result <- data.frame(
          Parameter = names(params),
          Estimate = as.numeric(params),
          Lower = NA,
          Upper = NA,
          stringsAsFactors = FALSE
        )

        # Round numeric values
        num_cols <- vapply(result, is.numeric, logical(1))
        result[num_cols] <- lapply(result[num_cols], round, digits = digits)

        return(result)
      }
    }
  } else if (inherits(model, "vb_mle")) {
    # MLE model handling
    make_row <- function(model, sex = NULL) {
      cf <- try(stats::coef(model), silent = TRUE)
      # Default NA CIs
      lower <- rep(NA_real_, 3)
      upper <- rep(NA_real_, 3)
      nm <- c("Linf", "k", "t0")

      if (!inherits(cf, "try-error")) {
        ci <- try(stats::confint(model), silent = TRUE)
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
        Lower = lower,
        Upper = upper,
        stringsAsFactors = FALSE
      )

      if (!is.null(sex)) out$Sex <- sex
      return(out)
    }

    if (is.list(model$model) && !inherits(model$model, "vb_optim")) {
      # Multiple models by sex
      res_list <- list()
      for (s in names(model$model)) {
        res_list[[s]] <- make_row(model$model[[s]], sex = s)
      }
      result <- do.call(rbind, res_list)
      # Ensure column order
      result <- result[, c("Sex", "Parameter", "Estimate", "Lower", "Upper")]
    } else {
      # Single model
      result <- make_row(model$model)
    }

    # Round only numeric columns
    num_cols <- vapply(result, is.numeric, logical(1))
    result[num_cols] <- lapply(result[num_cols], round, digits = digits)

    return(result)
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
