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
  # Convert age to numeric if it isn't already
  age <- as.numeric(age)
  
  # Basic input validation
  if(is.null(age) || length(age) == 0) {
    stop("'age' must be provided and have length > 0")
  }
  
  # Create a base data frame with the age data
  df <- data.frame(age = age)
  
  # Add year and sex if provided
  if(!is.null(year)) {
    if(length(year) != length(age)) {
      stop("'year' must have the same length as 'age'")
    }
    df$year <- as.factor(year)
  }
  
  if(!is.null(sex)) {
    if(length(sex) != length(age)) {
      stop("'sex' must have the same length as 'age'")
    }
    df$sex <- as.factor(sex)
  }
  
  # Remove NA values for age
  df <- df[!is.na(df$age), ]
  
  # Basic plot (age only)
  if(is.null(year) && is.null(sex)) {
    # Calculate counts
    counts <- as.data.frame(table(df$age))
    names(counts) <- c("age", "count")
    counts$age <- as.numeric(as.character(counts$age))
    
    # Create plot
    p <- ggplot2::ggplot(counts, ggplot2::aes(x = age, y = count)) +
      ggplot2::geom_col(fill = "steelblue") +
      ggplot2::labs(x = "Age", y = "Count") +
      theme_fn
  }
  # Plot by year
  else if(!is.null(year) && is.null(sex)) {
    # Calculate counts
    counts <- as.data.frame(table(df$age, df$year))
    names(counts) <- c("age", "year", "count")
    counts$age <- as.numeric(as.character(counts$age))
    
    # Create plot
    p <- ggplot2::ggplot(counts, ggplot2::aes(x = age, y = count, fill = year)) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::labs(x = "Age", y = "Count", fill = "Year") +
      theme_fn
  }
  # Plot by sex
  else if(is.null(year) && !is.null(sex)) {
    # Calculate counts
    counts <- as.data.frame(table(df$age, df$sex))
    names(counts) <- c("age", "sex", "count")
    counts$age <- as.numeric(as.character(counts$age))
    
    # Create plot
    p <- ggplot2::ggplot(counts, ggplot2::aes(x = age, y = count, fill = sex)) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::labs(x = "Age", y = "Count", fill = "Sex") +
      theme_fn
  }
  # Plot by year and sex
  else {
    # Calculate counts
    counts <- as.data.frame(table(df$age, df$year, df$sex))
    names(counts) <- c("age", "year", "sex", "count")
    counts$age <- as.numeric(as.character(counts$age))
    
    # Create plot
    p <- ggplot2::ggplot(counts, ggplot2::aes(x = age, y = count, fill = year)) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::facet_wrap(~sex) +
      ggplot2::labs(x = "Age", y = "Count", fill = "Year") +
      theme_fn
  }
  
  return(p)
}
