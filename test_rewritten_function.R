#!/usr/bin/env Rscript

# Complete rewrite of plot_vb_age_counts
plot_vb_age_counts_test <- function(age, year = NULL, sex = NULL, theme_fn = ggplot2::theme_minimal()) {
  # Load required package
  library(ggplot2)
  
  # Basic input validation
  if(is.null(age) || length(age) == 0) {
    stop("'age' must be provided and have length > 0")
  }
  
  # Convert age to numeric
  age <- as.numeric(age)
  cat("Age converted to numeric\n")
  
  # Basic plot (age only)
  if(is.null(year) && is.null(sex)) {
    cat("Creating basic plot (age only)\n")
    
    # Create table with counts
    cat("Creating table...\n")
    counts <- table(age)
    cat("Table created\n")
    
    # Convert to data frame
    cat("Converting to data frame...\n")
    df <- data.frame(
      age = as.numeric(names(counts)),
      count = as.numeric(counts)
    )
    cat("Data frame created\n")
    cat("Data frame structure:\n")
    print(str(df))
    
    # Create plot
    cat("Creating ggplot...\n")
    p <- ggplot(df, aes(x = age, y = count)) +
      geom_col(fill = "steelblue") +
      labs(x = "Age", y = "Count") +
      theme_fn
    cat("Plot created\n")
    
    return(p)
  }
  
  # Other cases omitted for brevity
  else {
    cat("Not testing other cases in this function\n")
    return(NULL)
  }
}

# Test the function
print("Loading AGE dataset...")
load("Example/age.rdata")
print("Dataset loaded")

print("AGE dataset structure:")
str(AGE)

print("\nTesting with new function:")
result <- tryCatch({
  p <- plot_vb_age_counts_test(AGE$age)
  print(p)
  print("Test successful!")
  TRUE
}, error = function(e) {
  print(paste("Error:", e$message))
  FALSE
})

if(result) {
  print("Basic plot works with new implementation!")
} else {
  print("Basic plot failed with new implementation.")
}
