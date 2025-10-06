#!/usr/bin/env Rscript

# Debug script focused on basic plot
print("Loading AGE dataset...")
load("Example/age.rdata")
print("Dataset loaded")

# Function to test basic plot directly
test_basic_plot <- function() {
  print("Original age data type:")
  print(class(AGE$age))
  print("Sample values:")
  print(head(AGE$age))
  
  # Try the direct approach from our function
  print("Creating table and converting to data frame...")
  age <- as.numeric(AGE$age)
  
  # Debug detailed steps
  print("Age after conversion to numeric:")
  print(head(age))
  print(class(age))
  
  # Create table directly
  print("Creating table...")
  age_table <- table(age)
  print("Table created")
  print("Table structure:")
  print(str(age_table))
  
  # Try extracting names and counts
  print("Extracting names and counts...")
  age_values <- as.numeric(names(age_table))
  count_values <- as.numeric(age_table)
  print("Values extracted")
  
  # Create data frame
  print("Creating data frame...")
  counts_df <- data.frame(
    age = age_values,
    count = count_values
  )
  print("Data frame created")
  print("Structure:")
  print(str(counts_df))
  
  # Try plotting
  print("Creating plot...")
  library(ggplot2)
  p <- ggplot2::ggplot(counts_df, ggplot2::aes(x = age, y = count)) +
    ggplot2::geom_col(fill = "steelblue") +
    ggplot2::labs(x = "Age", y = "Count") +
    ggplot2::theme_minimal()
  
  print("Plot created successfully!")
  return(p)
}

# Run the test
tryCatch({
  test_basic_plot()
  print("Success!")
}, error = function(e) {
  print(paste("Error:", e$message))
})
