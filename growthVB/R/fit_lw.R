#' Fit Length-Weight Relationship using Log-Linear Regression
#'
#' This function fits a length-weight relationship of the form W = a * L^b using
#' log-linear regression. The model is fitted as log(W) = log(a) + b * log(L),
#' with optional bias adjustment for back-transformation.
#'
#' @param length A numeric vector of lengths
#' @param weight A numeric vector of weights
#' @param bias_adjustment Logical indicating whether to apply bias adjustment
#'   when back-transforming the intercept (default TRUE). The bias adjustment
#'   corrects for the systematic underestimation that occurs when exponentiating
#'   the mean of log-transformed data.
#'
#' @return An object of class \code{c("lw_fit", "glm", "lm")} containing:
#'   \item{coefficients}{Model coefficients on log scale}
#'   \item{bias_adjustment}{Logical indicating if bias adjustment was applied}
#'   \item{lw_params}{Named vector with 'a' and 'b' parameters on original scale}
#'
#' @details
#' The length-weight relationship is a fundamental allometric relationship in
#' fisheries biology, typically expressed as:
#'
#' \deqn{W = a \cdot L^b}
#'
#' where W is weight, L is length, 'a' is a scaling coefficient, and 'b' is the
#' allometric exponent. For isometric growth, b = 3.
#'
#' The model is fitted using ordinary least squares regression on log-transformed
#' data:
#'
#' \deqn{\log(W) = \log(a) + b \cdot \log(L)}
#'
#' When bias_adjustment = TRUE, the intercept 'a' is calculated as:
#'
#' \deqn{a = \exp(\hat{\alpha} + \frac{\sigma^2}{2})}
#'
#' where \eqn{\hat{\alpha}} is the estimated intercept and \eqn{\sigma^2} is the
#' residual variance. This correction accounts for the fact that
#' \eqn{E[\exp(X)] \neq \exp(E[X])} for random variables.
#'
#' @examples
#' \dontrun{
#' # Simulate length-weight data
#' set.seed(123)
#' length <- runif(100, 10, 50)
#' weight <- 0.01 * length^3 * exp(rnorm(100, 0, 0.1))
#'
#' # Fit the model
#' fit <- fit_lw(length = length, weight = weight)
#'
#' # Get parameters
#' get_lw_params(fit)
#'
#' # Predict weights for new lengths
#' predict_lw(fit, new_length = c(20, 30, 40))
#' }
#'
#' @seealso \code{\link{get_lw_params}}, \code{\link{predict_lw}}, \code{\link{summary.lw_fit}}
#'
#' @references
#' Froese, R. (2006). Cube law, condition factor and weight-length relationships:
#' history, meta-analysis and recommendations. Journal of Applied Ichthyology,
#' 22(4), 241-253.
#'
#' @export
fit_lw <- function(length, weight, bias_adjustment = TRUE) {
  # Input validation
  if (!is.numeric(length) || !is.numeric(weight)) {
    stop("'length' and 'weight' must be numeric vectors")
  }
  if (length(length) != length(weight)) {
    stop("'length' and 'weight' must have the same length")
  }
  if (any(length <= 0, na.rm = TRUE) || any(weight <= 0, na.rm = TRUE)) {
    stop("All 'length' and 'weight' values must be positive")
  }

  # Remove NA values
  complete_cases <- complete.cases(length, weight)
  if (sum(complete_cases) < length(length)) {
    message(sprintf(
      "Removed %d observations with missing values",
      length(length) - sum(complete_cases)
    ))
  }
  length <- length[complete_cases]
  weight <- weight[complete_cases]

  if (length(length) < 3) {
    stop("At least 3 complete observations are required")
  }

  # Log transform
  log_length <- log(length)
  log_weight <- log(weight)

  # Fit linear model
  result <- stats::glm(log_weight ~ log_length)

  # Extract coefficients
  b <- stats::coef(result)[2]

  # Calculate 'a' with or without bias adjustment
  # Get residual variance (sigma^2) for bias correction
  sigma2 <- summary(result)$dispersion

  if (bias_adjustment) {
    # Bias-corrected back-transformation using residual variance
    # E[exp(X)] = exp(mu + sigma^2/2) for log-normal distribution
    a <- exp(stats::coef(result)[[1]] + sigma2 / 2)
  } else {
    a <- exp(stats::coef(result)[1])
  }

  # Store parameters and settings
  result$bias_adjustment <- bias_adjustment
  result$lw_params <- c(a = unname(a), b = unname(b))
  result$original_data <- data.frame(length = length, weight = weight)

  # Add custom class
  class(result) <- c("lw_fit", class(result))

  return(result)
}


#' Get Length-Weight Parameters from Fitted Model
#'
#' Extracts the length-weight relationship parameters (a and b) from a fitted
#' model object. Works with objects from \code{\link{fit_lw}} as well as
#' standard \code{glm} or \code{gam} objects fitted to log-transformed data.
#'
#' @param model A fitted model object, typically from \code{\link{fit_lw}},
#'   or a \code{glm}/\code{gam} object fitted with log(weight) ~ log(length)
#' @param bias_adjustment Logical indicating whether to apply bias adjustment.
#'   If NULL (default), uses the setting stored in the model object (for lw_fit
#'   objects) or defaults to TRUE.
#'
#' @return A named numeric vector with elements 'a' and 'b'
#'
#' @examples
#' \dontrun{
#' # Using fit_lw
#' fit <- fit_lw(length = lengths, weight = weights)
#' params <- get_lw_params(fit)
#' cat(sprintf("W = %.6f * L^%.3f\n", params["a"], params["b"]))
#'
#' # Using a standard glm
#' fit_glm <- glm(log(weight) ~ log(length), data = mydata)
#' params <- get_lw_params(fit_glm, bias_adjustment = TRUE)
#' }
#'
#' @seealso \code{\link{fit_lw}}, \code{\link{predict_lw}}
#'
#' @export
get_lw_params <- function(model, bias_adjustment = NULL) {
  # Handle bias_adjustment
  if (is.null(bias_adjustment)) {
    bias_adjustment <- model$bias_adjustment
    if (is.null(bias_adjustment)) {
      stop("Please specify bias_adjustment (TRUE/FALSE)")
    }
  }

  # Note for GAM users
  if ("gam" %in% class(model)) {
    message("Note: For GAM models, this assumes log(length) is the first smooth/linear term")
  }

  # Calculate parameters
  # Get residual variance (sigma^2) for bias correction
  sigma2 <- summary(model)$dispersion

  if (bias_adjustment) {
    # Bias-corrected back-transformation using residual variance
    a <- exp(stats::coef(model)[[1]] + sigma2 / 2)
    message("Using log-transformation bias adjustment")
  } else {
    a <- exp(stats::coef(model)[1])
  }

  b <- stats::coef(model)[2]

  params <- c(a = unname(a), b = unname(b))
  return(params)
}


#' Predict Weight from Length using Length-Weight Relationship
#'
#' Calculates predicted weights for given lengths using the allometric
#' length-weight relationship W = a * L^b.
#'
#' @param model A fitted model object from \code{\link{fit_lw}}, or parameters
#'   can be provided directly via \code{a} and \code{b} arguments
#' @param new_length A numeric vector of lengths for which to predict weights
#' @param a Optional scaling coefficient. If NULL, extracted from model
#' @param b Optional allometric exponent. If NULL, extracted from model
#'
#' @return A data frame with columns:
#'   \item{length}{Input length values}
#'   \item{weight}{Predicted weight values}
#'
#' @examples
#' \dontrun{
#' # Using a fitted model
#' fit <- fit_lw(length = lengths, weight = weights)
#' predictions <- predict_lw(fit, new_length = seq(10, 50, by = 5))
#'
#' # Using parameters directly
#' predictions <- predict_lw(new_length = seq(10, 50, by = 5), a = 0.01, b = 3)
#' }
#'
#' @seealso \code{\link{fit_lw}}, \code{\link{get_lw_params}}
#'
#' @export
predict_lw <- function(model = NULL, new_length, a = NULL, b = NULL) {
  # Get parameters from model if not provided directly
  if (is.null(a) || is.null(b)) {
    if (is.null(model)) {
      stop("Either 'model' or both 'a' and 'b' parameters must be provided")
    }

    if (inherits(model, "lw_fit") && !is.null(model$lw_params)) {
      params <- model$lw_params
    } else {
      params <- get_lw_params(model)
    }

    if (is.null(a)) a <- params["a"]
    if (is.null(b)) b <- params["b"]
  }

  # Input validation
  if (!is.numeric(new_length)) {
    stop("'new_length' must be a numeric vector")
  }
  if (any(new_length <= 0, na.rm = TRUE)) {
    warning("Negative or zero lengths will produce invalid predictions")
  }

  # Calculate predicted weights
  weight <- a * new_length^b

  return(data.frame(length = new_length, weight = weight))
}


#' Print Method for Length-Weight Fit Objects
#'
#' @param x An object of class \code{lw_fit}
#' @param ... Additional arguments (not used)
#'
#' @return Invisibly returns the object
#'
#' @export
print.lw_fit <- function(x, ...) {
  cat("Length-Weight Relationship Fit\n")
  cat("==============================\n\n")

  params <- x$lw_params
  cat(sprintf("Model: W = a * L^b\n\n"))
  cat(sprintf("Parameters:\n"))
  cat(sprintf("  a = %g\n", params["a"]))
  cat(sprintf("  b = %.4f\n", params["b"]))
  cat(sprintf("\nBias adjustment: %s\n", ifelse(x$bias_adjustment, "Yes", "No")))
  cat(sprintf("Sample size: %d\n", nrow(x$original_data)))

  invisible(x)
}


#' Summary Method for Length-Weight Fit Objects
#'
#' @param object An object of class \code{lw_fit}
#' @param ... Additional arguments (not used)
#'
#' @return A list containing summary information
#'
#' @export
summary.lw_fit <- function(object, ...) {
  params <- object$lw_params

  # Get model summary statistics
  model_summary <- summary.glm(object)

  cat("Length-Weight Relationship Summary\n")
  cat("==================================\n\n")

  cat("Model: W = a * L^b\n")
  cat("Fitted as: log(W) = log(a) + b * log(L)\n\n")

  cat("Parameters (original scale):\n")
  cat(sprintf("  a = %g\n", params["a"]))
  cat(sprintf("  b = %.4f\n", params["b"]))

  cat(sprintf(
    "\nBias adjustment applied: %s\n",
    ifelse(object$bias_adjustment, "Yes", "No")
  ))

  cat("\nRegression coefficients (log scale):\n")
  print(stats::coef(model_summary))

  cat(sprintf(
    "\nR-squared: %.4f\n",
    1 - model_summary$deviance / model_summary$null.deviance
  ))
  cat(sprintf(
    "Residual standard error: %.4f\n",
    sqrt(model_summary$dispersion)
  ))
  cat(sprintf("Sample size: %d\n", nrow(object$original_data)))

  # Isometry test
  b_se <- model_summary$coefficients[2, 2]
  t_stat <- (params["b"] - 3) / b_se
  p_value <- 2 * stats::pt(abs(t_stat), df = object$df.residual, lower.tail = FALSE)

  cat("\nTest for isometric growth (H0: b = 3):\n")
  cat(sprintf("  t = %.3f, p = %.4f\n", t_stat, p_value))
  if (p_value < 0.05) {
    if (params["b"] < 3) {
      cat("  Result: Significant negative allometry (b < 3)\n")
    } else {
      cat("  Result: Significant positive allometry (b > 3)\n")
    }
  } else {
    cat("  Result: No significant departure from isometry\n")
  }

  invisible(list(
    params = params,
    bias_adjustment = object$bias_adjustment,
    r_squared = 1 - model_summary$deviance / model_summary$null.deviance,
    residual_se = sqrt(model_summary$dispersion),
    n = nrow(object$original_data),
    isometry_test = list(t = t_stat, p = p_value, b = params["b"])
  ))
}


#' Plot Method for Length-Weight Fit Objects
#'
#' Creates a diagnostic plot showing the fitted length-weight relationship
#' with observed data points.
#'
#' @param x An object of class \code{lw_fit}
#' @param log_scale Logical indicating whether to plot on log-log scale (default FALSE)
#' @param ... Additional arguments passed to plotting functions
#'
#' @return A ggplot2 object (if ggplot2 is available) or base R plot
#'
#' @export
plot.lw_fit <- function(x, log_scale = FALSE, ...) {
  data <- x$original_data
  params <- x$lw_params

  # Generate prediction line (ensure positive lengths for log transform)
  min_length <- max(min(data$length), 0.001)
  length_seq <- seq(min_length, max(data$length), length.out = 100)
  pred_data <- predict_lw(x, new_length = length_seq)

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    if (log_scale) {
      p <- ggplot2::ggplot(data, ggplot2::aes(x = log(length), y = log(weight))) +
        ggplot2::geom_point(alpha = 0.6) +
        ggplot2::geom_line(
          data = pred_data,
          ggplot2::aes(x = log(length), y = log(weight)),
          colour = "royalblue", linewidth = 1
        ) +
        ggplot2::labs(
          x = "log(Length)", y = "log(Weight)",
          title = "Length-Weight Relationship (log scale)",
          subtitle = sprintf(
            "log(W) = %.4f + %.4f * log(L)",
            log(params["a"]), params["b"]
          )
        )
    } else {
      p <- ggplot2::ggplot(data, ggplot2::aes(x = length, y = weight)) +
        ggplot2::geom_point(alpha = 0.6) +
        ggplot2::geom_line(
          data = pred_data,
          ggplot2::aes(x = length, y = weight),
          colour = "royalblue", linewidth = 1
        ) +
        ggplot2::labs(
          x = "Length", y = "Weight",
          title = "Length-Weight Relationship",
          subtitle = sprintf(
            "W = %.2e * L^%.3f",
            params["a"], params["b"]
          )
        )
    }
    return(p)
  } else {
    # Base R fallback
    if (log_scale) {
      plot(log(data$length), log(data$weight),
        xlab = "log(Length)", ylab = "log(Weight)",
        main = "Length-Weight Relationship (log scale)", ...
      )
      lines(log(pred_data$length), log(pred_data$weight), col = "royalblue", lwd = 2)
    } else {
      plot(data$length, data$weight,
        xlab = "Length", ylab = "Weight",
        main = "Length-Weight Relationship", ...
      )
      lines(pred_data$length, pred_data$weight, col = "royalblue", lwd = 2)
    }
  }
}
