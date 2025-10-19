#' Plot von Bertalanffy Growth Curve
#'
#' This function creates a plot of the von Bertalanffy growth curve fitted to age and length data.
#'
#' @param model A model object returned by fit_vb_mle() or fit_vb_brms()
#' @param point_alpha Transparency of data points (default 0.5)
#' @param show_ci Whether to display confidence intervals (default TRUE)
#'
#' @return A ggplot2 object
#'
#' @examples
#' \dontrun{
#' # Simple example with simulated data
#' age <- 1:15
#' length <- 100 * (1 - exp(-0.2 * (age - (-0.5)))) + rnorm(15, 0, 5)
#' fit <- fit_vb_mle(age = age, length = length)
#' plot_vb(fit)
#' }
#'
#' @export
plot_vb <- function(model, point_alpha = 0.5, show_ci = TRUE) {
  if (inherits(model, "vb_mle")) {
    # Plot for MLE model
    if ("sex" %in% names(model$data)) {
      # Model with sex
      p <- ggplot2::ggplot(model$data, ggplot2::aes(x = age, y = length, colour = sex)) +
        ggplot2::geom_point(alpha = 0.8, size = 1.2) +
        ggplot2::geom_line(data = model$fits, ggplot2::aes(x = age, y = mean), linewidth = 1)

      if (show_ci) {
        # Only add confidence intervals if they're not all NA
        if (sum(!is.na(model$fits$lowerCI)) > 0 && sum(!is.na(model$fits$upperCI)) > 0) {
          p <- p + ggplot2::geom_line(
            data = model$fits[!is.na(model$fits$lowerCI), ],
            ggplot2::aes(x = age, y = lowerCI),
            linetype = "dashed"
          ) +
            ggplot2::geom_line(
              data = model$fits[!is.na(model$fits$upperCI), ],
              ggplot2::aes(x = age, y = upperCI),
              linetype = "dashed"
            )
        }
      }

      p <- p + ggplot2::facet_wrap(~sex) +
        ggplot2::labs(x = "Age", y = "Length") +
        ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
        ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
    } else {
      # Single model
      p <- ggplot2::ggplot(model$data, ggplot2::aes(x = age, y = length)) +
        ggplot2::geom_point(alpha = 0.8, size = 1.2, colour = "royalblue") +
        ggplot2::geom_line(
          data = model$fits, ggplot2::aes(x = age, y = mean),
          colour = "royalblue", linewidth = 1
        )

      if (show_ci) {
        # Only add confidence intervals if they're not all NA
        if (sum(!is.na(model$fits$lowerCI)) > 0 && sum(!is.na(model$fits$upperCI)) > 0) {
          p <- p + ggplot2::geom_line(
            data = model$fits[!is.na(model$fits$lowerCI), ],
            ggplot2::aes(x = age, y = lowerCI),
            colour = "royalblue", linetype = "dashed"
          ) +
            ggplot2::geom_line(
              data = model$fits[!is.na(model$fits$upperCI), ],
              ggplot2::aes(x = age, y = upperCI),
              colour = "royalblue", linetype = "dashed"
            )
        }
      }

      p <- p + ggplot2::labs(x = "Age", y = "Length") +
        ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
        ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
    }
  } else if (inherits(model, "vb_brms")) {
    # Plot for brms model
    if (is.list(model$predictions) && !is.data.frame(model$predictions)) {
      # Multiple models by sex
      # Combine predictions
      preds_combined <- do.call(rbind, model$predictions)

      # Create plot
      p <- ggplot2::ggplot(data = preds_combined) +
        ggplot2::geom_ribbon(ggplot2::aes(x = age, ymin = Q2.5, ymax = Q97.5, fill = Sex), alpha = 0.3) +
        ggplot2::geom_line(ggplot2::aes(x = age, y = Estimate, colour = Sex), linewidth = 1)

      # Add original data points if available
      if (!is.null(model$models[[1]]$data)) {
        data_list <- list()
        for (s in names(model$models)) {
          df <- data.frame(model$models[[s]]$data, Sex = s)
          data_list[[s]] <- df
        }
        data_combined <- do.call(rbind, data_list)

        p <- p + ggplot2::geom_point(
          data = data_combined,
          ggplot2::aes(x = age, y = length, colour = Sex),
          alpha = 0.8, size = 1.2
        )
      }

      p <- p + ggplot2::labs(x = "Age", y = "Length") +
        ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
        ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
    } else {
      # Single model
      p <- ggplot2::ggplot(data = model$predictions) +
        ggplot2::geom_ribbon(ggplot2::aes(x = age, ymin = Q2.5, ymax = Q97.5),
          alpha = 0.3, fill = "royalblue"
        ) +
        ggplot2::geom_line(ggplot2::aes(x = age, y = Estimate),
          colour = "royalblue", linewidth = 1
        )

      # Add original data points if available
      if (!is.null(model$models$data)) {
        p <- p + ggplot2::geom_point(
          data = model$models$data,
          ggplot2::aes(x = age, y = length),
          alpha = 0.8, size = 1.2, colour = "royalblue"
        )
      }

      p <- p + ggplot2::labs(x = "Age", y = "Length") +
        ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
        ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
    }
  } else {
    stop("Input must be a model object from fit_vb_mle() or fit_vb_brms()")
  }

  return(p)
}
