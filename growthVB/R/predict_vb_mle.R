#' Predict method for vb_mle objects
#'
#' Generates predicted lengths for a given set of ages using a fitted von Bertalanffy growth model.
#'
#' @param object A fitted model object of class vb_mle
#' @param newdata A data frame containing required variables for prediction:
#'   \itemize{
#'     \item For single models: 'age' (numeric)
#'     \item For sex-specific models: 'age' (numeric) and 'sex' (factor/character)
#'   }
#' @param interval Type of interval calculation: "none", "confidence", or "prediction"
#' @param level Confidence/prediction interval level (default 0.95)
#' @param ... Additional arguments (not used)
#'
#' @return A data frame with predicted values and optional intervals
#'
#' @export
predict.vb_mle <- function(object, newdata = NULL, interval = c("confidence", "none", "prediction"),
                           level = 0.95, ...) {
  interval <- match.arg(interval)

  # If no newdata provided, use the original data
  if (is.null(newdata)) {
    if (!is.null(object$data)) {
      newdata <- data.frame(age = object$data$age)
    } else {
      stop("No newdata provided and no data stored in model object")
    }
  }

  # Validate newdata structure
  if (!is.data.frame(newdata)) {
    stop("newdata must be a data frame")
  }

  if (nrow(newdata) == 0) {
    stop("newdata must contain at least one row")
  }

  # Check if we have multiple models (sex-specific)
  is_sex_specific <- is.list(object$model) && !inherits(object$model, "vb_optim")

  # Validate required variables based on model type
  if (is_sex_specific) {
    # Sex-specific models require both 'age' and 'sex'
    if (!"age" %in% names(newdata)) {
      stop("Sex-specific model requires 'age' variable in newdata for prediction")
    }
    if (!"sex" %in% names(newdata)) {
      stop("Sex-specific model requires 'sex' variable in newdata for prediction")
    }

    # Check for missing values in required variables
    if (any(is.na(newdata$age))) {
      stop("newdata contains missing values in 'age' variable")
    }
    if (any(is.na(newdata$sex))) {
      stop("newdata contains missing values in 'sex' variable")
    }

    # Check if any sex levels in newdata exist in the fitted model
    available_sex <- names(object$model)
    requested_sex <- unique(newdata$sex)
    missing_sex <- setdiff(requested_sex, available_sex)

    if (length(missing_sex) == length(requested_sex)) {
      stop(
        "No models available for any of the requested sex level(s): ",
        paste(requested_sex, collapse = ", "),
        ". Available sex levels: ", paste(available_sex, collapse = ", ")
      )
    }

    if (length(missing_sex) > 0) {
      warning(
        "Some sex levels not available in fitted model and will be ignored: ",
        paste(missing_sex, collapse = ", "),
        ". Available sex levels: ", paste(available_sex, collapse = ", ")
      )
    }
  } else {
    # Single models require only 'age'
    if (!"age" %in% names(newdata)) {
      stop("Single model requires 'age' variable in newdata for prediction")
    }

    # Check for missing values in age
    if (any(is.na(newdata$age))) {
      stop("newdata contains missing values in 'age' variable")
    }
  }

  if (is_sex_specific) {
    # Predict for each sex
    pred_list <- list()
    for (s in names(object$model)) {
      subset_data <- newdata[newdata$sex == s, , drop = FALSE]
      if (nrow(subset_data) > 0) {
        subset_pred <- predict_single_model(object$model[[s]], subset_data, interval, level)
        subset_pred$sex <- s
        pred_list[[s]] <- subset_pred
      }
    }

    # Combine predictions
    if (length(pred_list) > 0) {
      return(do.call(rbind, pred_list))
    } else {
      # Check what sex levels were requested vs available
      requested_sex <- unique(newdata$sex)
      available_sex <- names(object$model)
      missing_sex <- setdiff(requested_sex, available_sex)

      if (length(missing_sex) > 0) {
        stop(
          "No models available for sex level(s): ", paste(missing_sex, collapse = ", "),
          ". Available sex levels: ", paste(available_sex, collapse = ", ")
        )
      } else {
        return(NULL)
      }
    }
  } else {
    # Single model
    return(predict_single_model(object$model, newdata, interval, level))
  }
}

#' Helper function for prediction
#'
#' @param model A single vb_optim model
#' @param newdata Data frame with ages
#' @param interval Type of interval
#' @param level Confidence level
#'
#' @return Data frame with predictions
#'
#' @keywords internal
predict_single_model <- function(model, newdata, interval, level) {
  # Check inputs
  if (is.null(model) || is.null(newdata)) {
    stop("model and newdata cannot be NULL")
  }

  if (is.null(model$parameters)) {
    stop("model must contain 'parameters' component")
  }

  if (!"age" %in% names(newdata)) {
    stop("newdata must contain 'age' column")
  }

  # Extract parameters
  Linf <- unname(model$parameters["Linf"])
  k <- unname(model$parameters["k"])
  t0 <- unname(model$parameters["t0"])
  cv <- unname(model$parameters["cv"])

  # Check for missing parameters
  if (any(is.na(c(Linf, k, t0, cv)))) {
    stop("Missing required parameters: Linf, k, t0, cv")
  }

  # Calculate mean predictions using von Bertalanffy equation
  pred_mean <- Linf * (1 - exp(-k * (newdata$age - t0)))

  # Prepare output
  result <- data.frame(
    age = newdata$age,
    mean = pred_mean
  )

  # Calculate intervals if requested
  if (interval != "none") {
    # Degrees of freedom for t-distribution
    df <- nrow(model$data) - length(model$parameters)

    # Critical value for the desired confidence level
    t_value <- stats::qt(1 - (1 - level) / 2, df = df)

    # Calculate confidence intervals using delta method
    se_pred <- rep(NA_real_, length(newdata$age))

    # Check if we have a Hessian matrix for proper standard errors
    if (!is.null(model$hessian) && !is.null(model$vcov)) {
      # Use proper delta method with variance-covariance matrix
      vcov_mat <- model$vcov[1:3, 1:3] # Only for Linf, k, t0

      for (i in seq_along(newdata$age)) {
        # Gradient of the VB function with respect to parameters
        grad <- c(
          1 - exp(-k * (newdata$age[i] - t0)), # dL/dLinf
          Linf * exp(-k * (newdata$age[i] - t0)) * (newdata$age[i] - t0), # dL/dk
          Linf * exp(-k * (newdata$age[i] - t0)) * k # dL/dt0
        )

        # Delta method: var(g(theta)) = grad^T * Sigma * grad
        se_pred[i] <- sqrt(as.numeric(t(grad) %*% vcov_mat %*% grad))
      }
    } else {
      # Fallback to simple approximation
      for (i in seq_along(newdata$age)) {
        se_pred[i] <- 0.05 * pred_mean[i]
      }
    }

    # Add to result
    result$se <- se_pred

    if (interval == "confidence") {
      # Confidence interval for the mean
      result$lowerCI <- pred_mean - t_value * se_pred
      result$upperCI <- pred_mean + t_value * se_pred
    } else if (interval == "prediction") {
      # Prediction interval includes observation error
      # For prediction intervals, we need to add the variance of individual observations
      # SD for individual observations = CV * mean
      sd_obs <- cv * pred_mean

      # Total prediction variance = variance of mean estimate + variance of individual observation
      pred_var <- se_pred^2 + sd_obs^2

      result$lowerCI <- pred_mean - t_value * sqrt(pred_var)
      result$upperCI <- pred_mean + t_value * sqrt(pred_var)
    }
  }

  return(result)
}
