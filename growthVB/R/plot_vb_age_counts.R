#' Plot frequency distribution of age samples
#'
#' This function creates bar charts of age frequency distributions, optionally separated by year and/or sex.
#'
#' @param age A numeric vector of age values.
#' @param year An optional vector of years (same length as age).
#' @param sex An optional vector of sex values (same length as age).
#' @param theme_fn A ggplot2 theme function (default: theme_minimal).
#' @return A ggplot2 object.
#' @export
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
  # Basic input validation - make sure we have something to plot
  if(is.null(age) || length(age) == 0) {
    stop("'age' must be provided and have length > 0")
  }
  
  # Basic plot without considering year or sex
  if(is.null(year) && is.null(sex)) {
    # Count occurrences of each age value
    age_values <- sort(unique(age))
    age_counts <- numeric(length(age_values))
    
    for(i in seq_along(age_values)) {
      age_counts[i] <- sum(age == age_values[i], na.rm = TRUE)
    }
    
    # Create data frame for plotting
    plot_data <- data.frame(
      age = age_values,
      count = age_counts
    )
    
    # Create the plot
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = age, y = count)) +
      ggplot2::geom_col(fill = "steelblue") +
      ggplot2::labs(x = "Age", y = "Count") +
      theme_fn
    
    return(p)
  }
  
  # With year but no sex
  if(!is.null(year) && is.null(sex)) {
    # Ensure year has same length as age
    if(length(year) != length(age)) {
      stop("'year' must have the same length as 'age'")
    }
    
    # Get unique combinations of age and year
    age_values <- sort(unique(age))
    year_values <- sort(unique(year))
    
    # Pre-allocate results matrix
    counts_matrix <- matrix(0, nrow = length(age_values), ncol = length(year_values))
    rownames(counts_matrix) <- age_values
    colnames(counts_matrix) <- year_values
    
    # Count occurrences
    for(i in seq_along(age)) {
      if(!is.na(age[i]) && !is.na(year[i])) {
        row_idx <- match(age[i], age_values)
        col_idx <- match(year[i], year_values)
        counts_matrix[row_idx, col_idx] <- counts_matrix[row_idx, col_idx] + 1
      }
    }
    
    # Convert to long format for plotting
    plot_data <- data.frame()
    for(i in seq_along(age_values)) {
      for(j in seq_along(year_values)) {
        plot_data <- rbind(plot_data, data.frame(
          age = age_values[i],
          year = factor(year_values[j]),
          count = counts_matrix[i, j]
        ))
      }
    }
    
    # Create plot
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = age, y = count, fill = year)) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::labs(x = "Age", y = "Count", fill = "Year") +
      theme_fn
    
    return(p)
  }
  
  # With sex but no year
  if(is.null(year) && !is.null(sex)) {
    # Ensure sex has same length as age
    if(length(sex) != length(age)) {
      stop("'sex' must have the same length as 'age'")
    }
    
    # Get unique combinations of age and sex
    age_values <- sort(unique(age))
    sex_values <- sort(unique(sex))
    
    # Pre-allocate results matrix
    counts_matrix <- matrix(0, nrow = length(age_values), ncol = length(sex_values))
    rownames(counts_matrix) <- age_values
    colnames(counts_matrix) <- sex_values
    
    # Count occurrences
    for(i in seq_along(age)) {
      if(!is.na(age[i]) && !is.na(sex[i])) {
        row_idx <- match(age[i], age_values)
        col_idx <- match(sex[i], sex_values)
        counts_matrix[row_idx, col_idx] <- counts_matrix[row_idx, col_idx] + 1
      }
    }
    
    # Convert to long format for plotting
    plot_data <- data.frame()
    for(i in seq_along(age_values)) {
      for(j in seq_along(sex_values)) {
        plot_data <- rbind(plot_data, data.frame(
          age = age_values[i],
          sex = factor(sex_values[j]),
          count = counts_matrix[i, j]
        ))
      }
    }
    
    # Create plot
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = age, y = count, fill = sex)) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::labs(x = "Age", y = "Count", fill = "Sex") +
      theme_fn
    
    return(p)
  }
  
  # With both year and sex
  if(!is.null(year) && !is.null(sex)) {
    # Ensure lengths match
    if(length(year) != length(age) || length(sex) != length(age)) {
      stop("'year' and 'sex' must have the same length as 'age'")
    }
    
    # Get unique values
    age_values <- sort(unique(age))
    year_values <- sort(unique(year))
    sex_values <- sort(unique(sex))
    
    # For simplicity, we'll create a data frame with counts
    plot_data <- data.frame()
    
    # For each combination, count occurrences
    for(a in age_values) {
      for(y in year_values) {
        for(s in sex_values) {
          count <- sum(age == a & year == y & sex == s, na.rm = TRUE)
          plot_data <- rbind(plot_data, data.frame(
            age = a,
            year = factor(y),
            sex = factor(s),
            count = count
          ))
        }
      }
    }
    
    # Create plot
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = age, y = count, fill = year)) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::facet_wrap(~sex) +
      ggplot2::labs(x = "Age", y = "Count", fill = "Year") +
      theme_fn
    
    return(p)
  }
}
