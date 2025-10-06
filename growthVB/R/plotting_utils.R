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

  # Count observations in each bin
  if (!is.null(sex)) {
  counts <- as.data.frame(table(df$age_bin, df$length_bin, df$sex))
    names(counts) <- c("age_bin", "length_bin", "sex", "count")

    # Create plot
    p <- ggplot2::ggplot(counts, ggplot2::aes(x = age_bin, y = length_bin, fill = count)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient(name = "Count", low = "white", high = "steelblue") +
      ggplot2::labs(x = "Age", y = "Length") +
      ggplot2::facet_wrap(~sex) +
      theme_fn +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
  } else {
  counts <- as.data.frame(table(df$age_bin, df$length_bin))
    names(counts) <- c("age_bin", "length_bin", "count")

    # Create plot
    p <- ggplot2::ggplot(counts, ggplot2::aes(x = age_bin, y = length_bin, fill = count)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient(name = "Count", low = "white", high = "steelblue") +
      ggplot2::labs(x = "Age", y = "Length") +
      theme_fn +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
  }

  return(p)
}


