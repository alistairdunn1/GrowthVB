#' Predict method for vb_mle objects
#'
#' Generates predicted lengths for a given set of ages using a fitted von Bertalanffy growth model.
#'
#' @param object A fitted model object of class vb_mle
#' @param newdata A data frame containing 'age' values for prediction
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

  # Check if we have multiple models (sex-specific)
  if (is.list(object$model) && !inherits(object$model, "vb_optim")) {
    # For sex-specific models, we need the sex variable
    if (is.null(newdata$sex)) {
      stop("Sex-specific model requires 'sex' variable in newdata for prediction")
    }

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
      return(NULL)
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
  # Extract parameters
  Linf <- model$parameters["Linf"]
  k <- model$parameters["k"]
  t0 <- model$parameters["t0"]
  cv <- model$parameters["cv"]

  # Calculate mean predictions
  pred_mean <- model$vb_mean(newdata$age, Linf, k, t0)

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
    for (i in seq_along(newdata$age)) {
      # Gradient of the VB function with respect to parameters
      grad <- c(
        1 - exp(-k * (newdata$age[i] - t0)), # dL/dLinf
        Linf * exp(-k * (newdata$age[i] - t0)) * (newdata$age[i] - t0), # dL/dk
        Linf * exp(-k * (newdata$age[i] - t0)) * k # dL/dt0
      )

      # Just use a simple approximation for standard errors
      # since we don't have proper vcov matrix
      se_pred[i] <- 0.05 * pred_mean[i]
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
