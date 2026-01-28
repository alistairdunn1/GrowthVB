#' Plot Length-Weight Model Diagnostics
#'
#' This function creates diagnostic plots for length-weight relationship fits,
#' including residual plots, Q-Q plots, and other diagnostics to assess model
#' fit quality.
#'
#' @param model A model object returned by fit_lw()
#'
#' @return A list of ggplot2 objects:
#'   \item{residuals_vs_fitted}{Residuals vs fitted values plot}
#'   \item{qq_plot}{Normal Q-Q plot of residuals}
#'   \item{scale_location}{Scale-Location plot (sqrt of standardized residuals vs fitted)}
#'   \item{residuals_vs_length}{Residuals vs length plot}
#'
#' @details
#' The diagnostic plots help assess:
#' \itemize{
#'   \item \strong{Residuals vs Fitted}: Should show random scatter around zero.
#'     Patterns indicate potential model misspecification.
#'   \item \strong{Q-Q Plot}: Points should follow the diagonal line if residuals
#'     are normally distributed.
#'   \item \strong{Scale-Location}: Should show roughly constant spread.
#'     A funnel shape indicates heteroscedasticity.
#'   \item \strong{Residuals vs Length}: Should show random scatter.
#'     Patterns may indicate that the power relationship is inadequate.
#' }
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
#' # Generate diagnostics
#' diagnostics <- plot_lw_diagnostics(fit)
#'
#' # View individual plots
#' diagnostics$residuals_vs_fitted
#' diagnostics$qq_plot
#' diagnostics$scale_location
#' diagnostics$residuals_vs_length
#' }
#'
#' @seealso \code{\link{fit_lw}}, \code{\link{plot.lw_fit}}, \code{\link{summary.lw_fit}}
#'
#' @export
plot_lw_diagnostics <- function(model) {
  if (!inherits(model, "lw_fit")) {
    stop("This function only works with models from fit_lw()")
  }

  # Extract data and calculate residuals
  # Use [[ ]] to avoid partial matching issues with glm class
  data <- model[["original_data"]]
  params <- model[["lw_params"]]

  # Calculate fitted values (on original scale) and residuals
  fitted_log <- stats::fitted(model)
  residuals_log <- stats::residuals(model)

  # Calculate standardized residuals manually to avoid summary.lw_fit dispatch
  # rstandard uses: residuals / (sigma * sqrt(1 - h)) where h is leverage
  model_summary <- stats::summary.glm(model)
  sigma <- sqrt(model_summary$dispersion)
  leverage <- stats::hatvalues(model)
  student <- residuals_log / (sigma * sqrt(1 - leverage))

  # Create data frame for plotting
  # Use explicit [[ ]] extraction for data columns
  plot_data <- data.frame(
    length = data[["length"]],
    weight = data[["weight"]],
    log_length = log(data[["length"]]),
    log_weight = log(data[["weight"]]),
    fitted_log = fitted_log,
    fitted = exp(fitted_log),
    residual = residuals_log,
    student = student
  )

  # Initialize list for plots
  plots <- list()

  # 1. Residuals vs Fitted (on log scale)
  plots$residuals_vs_fitted <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = fitted_log, y = residual)
  ) +
    ggplot2::geom_point(colour = "royalblue", alpha = 0.8, size = 1.2) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
    ggplot2::geom_smooth(se = FALSE, colour = "black", method = "loess") +
    ggplot2::labs(
      x = "Fitted values (log scale)",
      y = "Residuals",
      title = "Residuals vs Fitted"
    ) +
    ggplot2::theme_minimal()

  # 2. Q-Q plot of residuals
  plots$qq_plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(sample = residual)
  ) +
    ggplot2::stat_qq(colour = "royalblue", alpha = 0.8, size = 1.2) +
    ggplot2::stat_qq_line(colour = "black") +
    ggplot2::labs(
      x = "Theoretical Quantiles",
      y = "Sample Quantiles",
      title = "Normal Q-Q Plot"
    ) +
    ggplot2::theme_minimal()

  # 3. Scale-Location Plot (sqrt of abs standardized residuals vs. fitted)
  plots$scale_location <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = fitted_log, y = sqrt(abs(student)))
  ) +
    ggplot2::geom_point(colour = "royalblue", alpha = 0.8, size = 1.2) +
    ggplot2::geom_smooth(se = FALSE, colour = "black", method = "loess") +
    ggplot2::labs(
      x = "Fitted values (log scale)",
      y = "sqrt(|Standardized residuals|)",
      title = "Scale-Location Plot"
    ) +
    ggplot2::theme_minimal()

  # 4. Residuals vs Length (to check for length-dependent patterns)
  plots$residuals_vs_length <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = length, y = residual)
  ) +
    ggplot2::geom_point(colour = "royalblue", alpha = 0.8, size = 1.2) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
    ggplot2::geom_smooth(se = FALSE, colour = "black", method = "loess") +
    ggplot2::labs(
      x = "Length",
      y = "Residuals",
      title = "Residuals vs Length"
    ) +
    ggplot2::theme_minimal()

  message("Available diagnostic plots: residuals_vs_fitted, qq_plot, scale_location, residuals_vs_length")

  return(plots)
}
