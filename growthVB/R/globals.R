# Declare global variables used in non-standard evaluation (e.g., ggplot2 aesthetics)
utils::globalVariables(c(
  # Common data columns
  "age", "length", "sex", "year", "Sex", "weight",
  # Fitted outputs
  "fitted", "residual", "student", "mean", "se",
  # Confidence intervals and predictions
  "lowerCI", "upperCI", "Estimate", "Q2.5", "Q97.5", "Model",
  # Binned summaries
  "age_bin", "length_bin", "count", "value"
))
