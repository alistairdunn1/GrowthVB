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
#' @param theme_fn Optional ggplot2 theme function to apply (default theme_minimal)
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
                                    theme_fn = ggplot2::theme_minimal()) {
  # Create data frame
  data <- data.frame(age = age, length = length)

  # Add sex if provided
  if (!is.null(sex)) {
    data$sex <- sex
    # Filter out NA values
    data <- data[!is.na(data$age) & !is.na(data$length) & !is.na(data$sex), ]
  } else {
    # Filter out NA values
    data <- data[!is.na(data$age) & !is.na(data$length), ]
  }

  # Define bins
  age_breaks <- seq(floor(min(data$age)), ceiling(max(data$age)) + bin_width_age, by = bin_width_age)
  length_breaks <- seq(floor(min(data$length)), ceiling(max(data$length)) + bin_width_length, by = bin_width_length)

  # Create binned data
  data$age_bin <- cut(data$age, breaks = age_breaks, include.lowest = TRUE, right = FALSE)
  data$length_bin <- cut(data$length, breaks = length_breaks, include.lowest = TRUE, right = FALSE)

  # Count observations in each bin
  if (!is.null(sex)) {
    counts <- as.data.frame(table(data$age_bin, data$length_bin, data$sex))
    names(counts) <- c("age_bin", "length_bin", "sex", "count")

    # Create plot
    p <- ggplot2::ggplot(counts, ggplot2::aes(x = age_bin, y = length_bin, fill = count)) +
      ggplot2::geom_tile() +
      viridis::scale_fill_viridis_c(name = "Count") +
      ggplot2::labs(x = "Age", y = "Length") +
      ggplot2::facet_wrap(~sex) +
      theme_fn +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
  } else {
    counts <- as.data.frame(table(data$age_bin, data$length_bin))
    names(counts) <- c("age_bin", "length_bin", "count")

    # Create plot
    p <- ggplot2::ggplot(counts, ggplot2::aes(x = age_bin, y = length_bin, fill = count)) +
      ggplot2::geom_tile() +
      viridis::scale_fill_viridis_c(name = "Count") +
      ggplot2::labs(x = "Age", y = "Length") +
      theme_fn +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
  }

  return(p)
}

#' Plot Age Count Summary
#'
#' This function creates a plot summarizing the number of observations by age (and optionally by year).
#'
#' @param age A numeric vector of ages
#' @param year An optional vector of years for each observation
#' @param sex An optional factor or character vector specifying the sex for each observation
#' @param theme_fn Optional ggplot2 theme function to apply (default theme_minimal)
#'
#' @return A ggplot2 object
#'
#' @examples
#' \dontrun{
#' # Simple example with simulated data
#' age <- 1:15
#' years <- sample(2020:2023, 15, replace = TRUE)
#' plot_vb_age_counts(age = age, year = years)
#' }
#'
#' @export
plot_vb_age_counts <- function(age, year = NULL, sex = NULL,
                               theme_fn = ggplot2::theme_minimal()) {
  # Create data frame
  data <- data.frame(age = age)

  # Add year and/or sex if provided
  if (!is.null(year)) {
    data$year <- year
  }

  if (!is.null(sex)) {
    data$sex <- sex
  }

  # Filter out NA values for age
  data <- data[!is.na(data$age), ]

  # Create plot based on available data
  if (!is.null(year) && !is.null(sex)) {
    # Group by year and sex
    counts <- as.data.frame(table(data$age, data$year, data$sex))
    names(counts) <- c("age", "year", "sex", "count")

    p <- ggplot2::ggplot(counts, ggplot2::aes(x = age, y = count, fill = year)) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::facet_wrap(~sex) +
      ggplot2::labs(x = "Age", y = "Count", fill = "Year") +
      theme_fn
  } else if (!is.null(year)) {
    # Group by year only
    counts <- as.data.frame(table(data$age, data$year))
    names(counts) <- c("age", "year", "count")

    p <- ggplot2::ggplot(counts, ggplot2::aes(x = age, y = count, fill = year)) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::labs(x = "Age", y = "Count", fill = "Year") +
      theme_fn
  } else if (!is.null(sex)) {
    # Group by sex only
    counts <- as.data.frame(table(data$age, data$sex))
    names(counts) <- c("age", "sex", "count")

    p <- ggplot2::ggplot(counts, ggplot2::aes(x = age, y = count, fill = sex)) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::labs(x = "Age", y = "Count", fill = "Sex") +
      theme_fn
  } else {
    # Simple age counts
    counts <- as.data.frame(table(data$age))
    names(counts) <- c("age", "count")

    p <- ggplot2::ggplot(counts, ggplot2::aes(x = age, y = count)) +
      ggplot2::geom_col(fill = "steelblue") +
      ggplot2::labs(x = "Age", y = "Count") +
      theme_fn
  }

  return(p)
}
