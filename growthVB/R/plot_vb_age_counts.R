#' Plot frequency distribution of age samples
#'
#' This function creates bar charts of age frequency distributions, optionally separated by group and/or sex.
#'
#' @param age A numeric vector of age values.
#' @param group An optional vector of grouping variables (same length as age).
#' @param sex An optional vector of sex values (same length as age).
#' @return A ggplot2 object.
#' @export
#' @examples
#' \dontrun{
#' # Simple example with simulated data
#' age <- 1:15
#' groups <- sample(c("A", "B", "C"), 15, replace = TRUE)
#' plot_vb_age_counts(age = age, group = groups)
#' }
#'
#' @export
plot_vb_age_counts <- function(age, group = NULL, sex = NULL) {
  require(ggplot2)

  # Basic input validation
  if (is.null(age) || length(age) == 0) {
    stop("'age' must be provided and have length > 0")
  }

  # Convert age to numeric if it isn't already
  age <- as.numeric(age)

  # Basic input validation for group and sex
  if (!is.null(group) && length(group) != length(age)) {
    stop("'group' must have the same length as 'age'")
  }

  if (!is.null(sex) && length(sex) != length(age)) {
    stop("'sex' must have the same length as 'age'")
  }

  # Remove NA values from age
  valid_indices <- !is.na(age)
  age <- age[valid_indices]

  if (!is.null(group)) {
    group <- as.factor(group[valid_indices])
  }

  if (!is.null(sex)) {
    sex <- as.factor(sex[valid_indices])
  }

  # Basic plot (age only)
  if (is.null(group) && is.null(sex)) {
    # Create table directly from age vector
    age_table <- table(age)

    # Extract values and create data frame
    counts <- data.frame(
      age = as.numeric(names(age_table)),
      count = as.numeric(age_table)
    )

    # Create plot
    p <- ggplot2::ggplot(counts, ggplot2::aes(x = age, y = count)) +
      ggplot2::geom_col(fill = "steelblue") +
      ggplot2::labs(x = "Age", y = "Count")
  }
  # Plot by group
  else if (!is.null(group) && is.null(sex)) {
    # Calculate counts
    group_table <- table(age, group)
    counts <- as.data.frame(group_table)
    names(counts) <- c("age", "group", "count")
    counts$age <- as.numeric(as.character(counts$age))
    # Create plot
    p <- ggplot2::ggplot(counts, ggplot2::aes(x = age, y = count, fill = group)) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::labs(x = "Age", y = "Count", fill = "Group")
  }
  # Plot by sex
  else if (is.null(group) && !is.null(sex)) {
    # Calculate counts
    sex_table <- table(age, sex)
    counts <- as.data.frame(sex_table)
    names(counts) <- c("age", "sex", "count")
    counts$age <- as.numeric(as.character(counts$age))

    # Create plot
    p <- ggplot2::ggplot(counts, ggplot2::aes(x = age, y = count, fill = sex)) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::labs(x = "Age", y = "Count", fill = "Sex")
  }

  # With both group and sex
  else if (!is.null(group) && !is.null(sex)) {
    # Calculate counts
    group_sex_table <- table(age, group, sex)
    counts <- as.data.frame(group_sex_table)
    names(counts) <- c("age", "group", "sex", "count")
    counts$age <- as.numeric(as.character(counts$age))

    # Create plot
    p <- ggplot2::ggplot(counts, ggplot2::aes(x = age, y = count, fill = group)) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::facet_wrap(~sex) +
      ggplot2::labs(x = "Age", y = "Count", fill = "Group")
  }

  return(p)
}
