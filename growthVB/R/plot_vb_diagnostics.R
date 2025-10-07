#' Plot von Bertalanffy Growth Model Diagnostics
#'
#' This function creates diagnostic plots for von Bertalanffy growth model fits,
#' including residual plots, Q-Q plots, and other diagnostics.
#'
#' @param model A model object returned by fit_vb_mle()
#' @param theme_fn Optional ggplot2 theme function to apply (default theme_minimal)
#'
#' @return A list of ggplot2 objects
#'
#' @examples
#' \dontrun{
#' # Simple example with simulated data
#' age <- 1:15
#' length <- 100 * (1 - exp(-0.2 * (age - (-0.5)))) + rnorm(15, 0, 5)
#' fit <- fit_vb_mle(age = age, length = length)
#' diagnostics <- plot_vb_diagnostics(fit)
#' diagnostics$residuals_vs_fitted
#' }
#'
#' @export
plot_vb_diagnostics <- function(model, theme_fn = ggplot2::theme_minimal()) {
  if (!inherits(model, "vb_mle")) {
    stop("This function only works with models from fit_vb_mle()")
  }

  # Initialize list for plots
  plots <- list()

  # Function to create diagnostic plots for a specific model/sex
  create_diagnostic_plots <- function(data, model_obj, sex_label = NULL) {
    # Prepare plot title suffix
    title_suffix <- if (is.null(sex_label)) "" else paste(" -", sex_label)

    # Residuals vs Fitted
    p1 <- ggplot2::ggplot(data, ggplot2::aes(x = fitted, y = residual)) +
      ggplot2::geom_point() +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      ggplot2::geom_smooth(se = FALSE, color = "blue", method = "loess") +
      ggplot2::labs(
        x = "Fitted values",
        y = "Residuals",
        title = paste("Residuals vs Fitted", title_suffix)
      ) +
      theme_fn

    # Q-Q plot of residuals
    p2 <- ggplot2::ggplot(data, ggplot2::aes(sample = residual)) +
      ggplot2::stat_qq() +
      ggplot2::stat_qq_line() +
      ggplot2::labs(
        x = "Theoretical Quantiles",
        y = "Sample Quantiles",
        title = paste("Normal Q-Q Plot", title_suffix)
      ) +
      theme_fn

    # Scale-Location Plot (sqrt of abs residuals vs. fitted)
    p3 <- ggplot2::ggplot(data, ggplot2::aes(x = fitted, y = sqrt(abs(student)))) +
      ggplot2::geom_point() +
      ggplot2::geom_smooth(se = FALSE, color = "blue", method = "loess") +
      ggplot2::labs(
        x = "Fitted values",
        y = "sqrt(|Standardized residuals|)",
        title = paste("Scale-Location Plot", title_suffix)
      ) +
      theme_fn

    # Residuals vs Age
    p4 <- ggplot2::ggplot(data, ggplot2::aes(x = age, y = residual)) +
      ggplot2::geom_point() +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      ggplot2::geom_smooth(se = FALSE, color = "blue", method = "loess") +
      ggplot2::labs(
        x = "Age",
        y = "Residuals",
        title = paste("Residuals vs Age", title_suffix)
      ) +
      theme_fn

    return(list(
      residuals_vs_fitted = p1,
      qq_plot = p2,
      scale_location = p3,
      residuals_vs_age = p4
    ))
  }

  # Check if we have multiple models by sex
  if (is.list(model$model) && !inherits(model$model, "vb_optim")) {
    # Create diagnostic plots for each sex
    for (s in names(model$model)) {
      subset_data <- model$data[model$data$sex == s, ]
      subset_model <- model$model[[s]]

      plots[[s]] <- create_diagnostic_plots(subset_data, subset_model, s)
    }
  } else {
    # Single model
    plots <- create_diagnostic_plots(model$data, model$model)
  }

  return(plots)
}
