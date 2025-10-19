#' Plot Age-Length Heatmap
#'
#' This function creates a heatmap of the age-length observations to visualize
#' the distribution of samples.
#'
#' @param age A numeric vector of ages
#' @param length A numeric vector of lengths
#' @param sex An optional factor or character vector specifying the sex for each observation
#' @param bin_width_age Width of age bins (default 1)
#' @param bin_width_length Width of length bins (default 5)
#' @param add_smoother Logical, whether to add a smoother line showing the age-length relationship (default TRUE)
#'
#' @return A ggplot2 object
#'
#' @examples
#' \dontrun{
#' # Simple example with simulated data
#' age <- 1:15
#' length <- 100 * (1 - exp(-0.2 * (age - (-0.5)))) + rnorm(15, 0, 5)
#' plot_age_length_heatmap(age = age, length = length)
#' }
#'
#' @export
plot_age_length_heatmap <- function(age, length, sex = NULL,
                                    bin_width_age = 1, bin_width_length = 5,
                                    add_smoother = TRUE) {
  # Create data frame
  df <- data.frame(age = age, length = length)

  # Add sex if provided
  if (!is.null(sex)) {
    df$sex <- sex
    # Filter out NA values
    df <- df[!is.na(df$age) & !is.na(df$length) & !is.na(df$sex), ]
  } else {
    # Filter out NA values
    df <- df[!is.na(df$age) & !is.na(df$length), ]
  }

  # Define bins
  age_breaks <- seq(floor(min(df$age)), ceiling(max(df$age)) + bin_width_age, by = bin_width_age)
  length_breaks <- seq(floor(min(df$length)), ceiling(max(df$length)) + bin_width_length, by = bin_width_length)

  # Create binned data
  df$age_bin <- cut(df$age, breaks = age_breaks, include.lowest = TRUE, right = FALSE)
  df$length_bin <- cut(df$length, breaks = length_breaks, include.lowest = TRUE, right = FALSE)

  # Calculate midpoints for bins
  age_midpoints <- (age_breaks[-length(age_breaks)] + age_breaks[-1]) / 2
  length_midpoints <- (length_breaks[-length(length_breaks)] + length_breaks[-1]) / 2

  # Convert bins to midpoint values
  df$age_midpoint <- age_midpoints[as.numeric(df$age_bin)]
  df$length_midpoint <- length_midpoints[as.numeric(df$length_bin)]

  # Count observations in each bin
  if (!is.null(sex)) {
    counts <- as.data.frame(table(df$age_bin, df$length_bin, df$sex))
    names(counts) <- c("age_bin", "length_bin", "sex", "count")

    # Add midpoint values to counts data
    counts$age_midpoint <- age_midpoints[as.numeric(counts$age_bin)]
    counts$length_midpoint <- length_midpoints[as.numeric(counts$length_bin)]

    # Create plot using midpoint coordinates
    p <- ggplot2::ggplot(counts, ggplot2::aes(x = age_midpoint, y = length_midpoint, fill = count)) +
      ggplot2::geom_tile(width = bin_width_age * 0.9, height = bin_width_length * 0.9) +
      ggplot2::scale_fill_gradient(name = "Count", low = "white", high = "steelblue") +
      ggplot2::labs(x = "Age", y = "Length") +
      ggplot2::facet_wrap(~sex)

    # Add smoother if requested - now using same coordinate system
    if (add_smoother) {
      p <- p + ggplot2::geom_smooth(
        data = df,
        ggplot2::aes(x = age, y = length),
        method = "loess", se = TRUE, color = "red", linewidth = 1, alpha = 0.7,
        inherit.aes = FALSE
      )
    }
  } else {
    counts <- as.data.frame(table(df$age_bin, df$length_bin))
    names(counts) <- c("age_bin", "length_bin", "count")

    # Add midpoint values to counts data
    counts$age_midpoint <- age_midpoints[as.numeric(counts$age_bin)]
    counts$length_midpoint <- length_midpoints[as.numeric(counts$length_bin)]

    # Create plot using midpoint coordinates
    p <- ggplot2::ggplot(counts, ggplot2::aes(x = age_midpoint, y = length_midpoint, fill = count)) +
      ggplot2::geom_tile(width = bin_width_age * 0.9, height = bin_width_length * 0.9) +
      ggplot2::scale_fill_gradient(name = "Count", low = "white", high = "steelblue") +
      ggplot2::labs(x = "Age", y = "Length")

    # Add smoother if requested - now using same coordinate system
    if (add_smoother) {
      p <- p + ggplot2::geom_smooth(
        data = df,
        ggplot2::aes(x = age, y = length),
        method = "loess", se = TRUE, color = "red", linewidth = 1, alpha = 0.7,
        inherit.aes = FALSE
      )
    }
  }

  return(p)
}

#' Plot Empirical Coefficient of Variation by Age
#'
#' This function creates a plot showing the empirical coefficient of variation (CV)
#' of length measurements by age. This is useful for understanding how variability
#' in length changes with age and for validating CV models in von Bertalanffy fits.
#'
#' @param age A numeric vector of ages
#' @param length A numeric vector of lengths
#' @param sex An optional factor or character vector specifying the sex for each observation
#' @param min_n Minimum number of observations required per age group to calculate CV (default 3)
#' @param add_smoother Logical, whether to add a smoother line showing the CV trend (default TRUE)
#'
#' @return A ggplot2 object
#'
#' @examples
#' \dontrun{
#' # Simple example with simulated data
#' age <- rep(1:15, each = 10)
#' length <- 100 * (1 - exp(-0.2 * (age - (-0.5)))) +
#'   rnorm(length(age), 0, 5 + 0.1 * age) # CV increases with age
#' plot_empirical_cv(age = age, length = length)
#' }
#'
#' @export
plot_empirical_cv <- function(age, length, sex = NULL, min_n = 3, add_smoother = TRUE) {
  # Create data frame
  df <- data.frame(age = age, length = length)

  # Add sex if provided
  if (!is.null(sex)) {
    df$sex <- sex
    # Filter out NA values
    df <- df[!is.na(df$age) & !is.na(df$length) & !is.na(df$sex), ]
  } else {
    # Filter out NA values
    df <- df[!is.na(df$age) & !is.na(df$length), ]
  }

  # Function to calculate CV statistics by age
  calculate_cv_by_age <- function(data, group_var = NULL) {
    if (is.null(group_var)) {
      # Group by age only
      age_groups <- split(data, data$age)
    } else {
      # Group by both sex and age
      data$group_key <- paste(data[[group_var]], data$age, sep = "_")
      age_groups <- split(data, data$group_key)
    }

    # Calculate CV for each age group
    cv_list <- list()
    for (i in seq_along(age_groups)) {
      group_data <- age_groups[[i]]
      if (nrow(group_data) >= min_n) {
        n_obs <- nrow(group_data)
        mean_len <- mean(group_data$length, na.rm = TRUE)
        sd_len <- stats::sd(group_data$length, na.rm = TRUE)
        cv <- sd_len / mean_len

        if (!is.na(cv) && is.finite(cv)) {
          cv_row <- data.frame(
            age = unique(group_data$age),
            n = n_obs,
            mean_length = mean_len,
            sd_length = sd_len,
            cv = cv
          )

          # Add sex if grouping by sex
          if (!is.null(group_var)) {
            cv_row[[group_var]] <- unique(group_data[[group_var]])
          }

          cv_list[[i]] <- cv_row
        }
      }
    }

    if (length(cv_list) > 0) {
      cv_data <- do.call(rbind, cv_list)
      rownames(cv_data) <- NULL
      return(cv_data)
    } else {
      return(data.frame())
    }
  }

  # Calculate CV by age, with or without sex
  if (!is.null(sex) && length(unique(sex[!is.na(sex)])) > 1) {
    # Calculate CV separately for each sex
    cv_data <- calculate_cv_by_age(df, "sex")

    # Create plot with sex facets
    p <- ggplot2::ggplot(cv_data, ggplot2::aes(x = age, y = cv)) +
      ggplot2::geom_point(size = 1.2, alpha = 0.8, colour = "royalblue") +
      ggplot2::geom_line(alpha = 0.6, colour = "royalblue") +
      ggplot2::facet_wrap(~sex) +
      ggplot2::labs(
        x = "Age",
        y = "Coefficient of Variation (CV)",
        title = "Empirical CV of Length by Age",
        subtitle = paste("Minimum", min_n, "observations per age group")
      )

    # Add smoother if requested
    if (add_smoother) {
      p <- p + ggplot2::geom_smooth(
        method = "loess", se = TRUE, color = "red", linewidth = 1, alpha = 0.7
      )
    }
  } else {
    # Calculate CV for combined data
    cv_data <- calculate_cv_by_age(df, NULL)

    # Create plot without sex facets
    p <- ggplot2::ggplot(cv_data, ggplot2::aes(x = age, y = cv)) +
      ggplot2::geom_point(size = 1.2, alpha = 0.8, colour = "royalblue") +
      ggplot2::geom_line(alpha = 0.6, colour = "royalblue") +
      ggplot2::labs(
        x = "Age",
        y = "Coefficient of Variation (CV)",
        title = "Empirical CV of Length by Age",
        subtitle = paste("Minimum", min_n, "observations per age group")
      )

    # Add smoother if requested
    if (add_smoother) {
      p <- p + ggplot2::geom_smooth(
        method = "loess", se = TRUE, color = "red", linewidth = 1, alpha = 0.7
      )
    }
  }

  return(p)
}
