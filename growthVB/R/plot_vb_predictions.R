#' Plot von Bertalanffy Growth Model Predictions
#'
#' Creates a comprehensive plot showing the fitted von Bertalanffy growth curve(s),
#' prediction intervals, and original data points from predict.vb_mle() output.
#'
#' @param predictions A data frame output from predict.vb_mle() containing predictions
#' @param original_data Optional data frame with original data points to overlay.
#'   Should contain 'age', 'length', and optionally 'sex' columns.
#' @param show_points Logical. Whether to show original data points (default TRUE)
#' @param show_intervals Logical. Whether to show prediction/confidence intervals (default TRUE)
#' @param facet_by_sex Logical. Whether to create separate panels for each sex (default TRUE when sex present)
#' @param alpha_ribbon Numeric. Transparency for confidence/prediction intervals (default 0.3)
#' @param alpha_points Numeric. Transparency for data points (default 0.5)
#' @param point_size Numeric. Size of data points (default 0.9)
#' @param line_size Numeric. Size of fitted curve lines (default 1)
#' @param xlab Character. X-axis label (default "Age")
#' @param ylab Character. Y-axis label (default "Length (cm)")
#' @param title Character. Plot title (default NULL)
#'
#' @return A ggplot2 object
#'
#' @examples
#' \dontrun{
#' # Fit model and make predictions
#' fit <- fit_vb_mle(age = age, length = length, sex = sex)
#' new_ages <- data.frame(
#'   age = rep(seq(0, 50, by = 1), each = 2),
#'   sex = rep(c("M", "F"), times = 51)
#' )
#' predictions <- predict(fit, newdata = new_ages, interval = "prediction")
#'
#' # Create plot
#' plot_vb_predictions(predictions, original_data = fit$data)
#' }
#'
#' @export
plot_vb_predictions <- function(predictions,
                                original_data = NULL,
                                show_points = TRUE,
                                show_intervals = TRUE,
                                facet_by_sex = NULL,
                                alpha_ribbon = 0.3,
                                alpha_points = 0.5,
                                point_size = 0.9,
                                line_size = 1,
                                xlab = "Age",
                                ylab = "Length (cm)",
                                title = NULL) {
  # Check if ggplot2 is available
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for this function. Please install it.",
      call. = FALSE
    )
  }

  # Validate input
  if (!is.data.frame(predictions)) {
    stop("predictions must be a data frame")
  }

  required_cols <- c("age", "mean")
  missing_cols <- setdiff(required_cols, names(predictions))
  if (length(missing_cols) > 0) {
    stop("predictions must contain columns: ", paste(missing_cols, collapse = ", "))
  }

  # Check if we have sex-specific predictions
  has_sex <- "sex" %in% names(predictions)

  # Check if we have confidence/prediction intervals
  has_intervals <- all(c("lowerCI", "upperCI") %in% names(predictions))

  if (show_intervals && !has_intervals) {
    warning("Confidence/prediction intervals not available in predictions. Setting show_intervals = FALSE.")
    show_intervals <- FALSE
  }

  # Determine faceting
  if (is.null(facet_by_sex)) {
    facet_by_sex <- has_sex
  }

  # Validate original_data if provided
  if (!is.null(original_data)) {
    if (!is.data.frame(original_data)) {
      stop("original_data must be a data frame")
    }

    if (!"age" %in% names(original_data) || !"length" %in% names(original_data)) {
      stop("original_data must contain 'age' and 'length' columns")
    }

    # Check sex consistency
    if (has_sex && !"sex" %in% names(original_data)) {
      warning("predictions has sex variable but original_data does not. Data points may not align properly.")
    }

    if (!has_sex && "sex" %in% names(original_data)) {
      warning("original_data has sex variable but predictions does not. Consider providing sex-specific predictions.")
    }
  }

  # Create base plot
  if (has_sex) {
    p <- ggplot2::ggplot(predictions, ggplot2::aes(x = age, y = mean, colour = sex, fill = sex))
  } else {
    p <- ggplot2::ggplot(predictions, ggplot2::aes(x = age, y = mean))
  }

  # Add confidence/prediction intervals
  if (show_intervals && has_intervals) {
    if (has_sex) {
      p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = lowerCI, ymax = upperCI),
        alpha = alpha_ribbon, colour = NA
      )
    } else {
      p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = lowerCI, ymax = upperCI),
        alpha = alpha_ribbon, colour = NA, fill = "grey50"
      )
    }
  }

  # Add fitted curve
  p <- p + ggplot2::geom_line(linewidth = line_size)

  # Add original data points
  if (show_points && !is.null(original_data)) {
    if (has_sex && "sex" %in% names(original_data)) {
      p <- p + ggplot2::geom_point(
        data = original_data,
        ggplot2::aes(x = age, y = length, colour = sex),
        size = point_size, alpha = alpha_points
      )
    } else {
      p <- p + ggplot2::geom_point(
        data = original_data,
        ggplot2::aes(x = age, y = length),
        size = point_size, alpha = alpha_points,
        colour = "black", inherit.aes = FALSE
      )
    }
  }

  # Add faceting for sex if requested
  if (facet_by_sex && has_sex) {
    p <- p + ggplot2::facet_wrap(~sex)
  }

  # Add labels and formatting
  p <- p +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab) +
    ggplot2::ylim(0, NA) +
    ggplot2::xlim(0, NA)

  # Add title if provided
  if (!is.null(title)) {
    p <- p + ggplot2::ggtitle(title)
  }

  return(p)
}
