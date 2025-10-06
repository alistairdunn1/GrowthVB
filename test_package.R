# Load the growthVB package
library(growthVB)

# Set working directory to the project root
# setwd("C:/Users/alist/OneDrive/Projects/Software/GrowthVB")

# Load the test data
load("Example/age.rdata")

# Print a summary of the loaded data
cat("Data summary:\n")
str(data)

# 1. Plot the age count distribution
cat("\n1. Testing plot_vb_age_counts:\n")
p_counts <- plot_vb_age_counts(age = data$age)
print(p_counts)

# 2. Plot with year
cat("\n2. Testing plot_vb_age_counts with year:\n")
p_counts_year <- plot_vb_age_counts(age = data$age, year = data$year)
print(p_counts_year)

# 3. Plot with sex
cat("\n3. Testing plot_vb_age_counts with sex:\n")
p_counts_sex <- plot_vb_age_counts(age = data$age, sex = data$sex)
print(p_counts_sex)

# 4. Plot with both year and sex
cat("\n4. Testing plot_vb_age_counts with year and sex:\n")
p_counts_both <- plot_vb_age_counts(age = data$age, year = data$year, sex = data$sex)
print(p_counts_both)

# 5. Age-length heatmap
cat("\n5. Testing plot_age_length_heatmap:\n")
p_heatmap <- plot_age_length_heatmap(age = data$age, length = data$length)
print(p_heatmap)

# 6. Fit VB with nls
cat("\n6. Testing fit_vb_nls:\n")
vb_fit <- try(fit_vb_nls(age = data$age, length = data$length))

if(!inherits(vb_fit, "try-error")) {
  # 7. Summarize the model
  cat("\n7. Testing summarize_vb:\n")
  summary_nls <- summarize_vb(vb_fit)
  print(summary_nls)
  
  # 8. Plot the model fit
  cat("\n8. Testing plot_vb:\n")
  p_fit <- plot_vb(vb_fit)
  print(p_fit)
  
  # 9. Plot diagnostics
  cat("\n9. Testing plot_vb_diagnostics:\n")
  p_diag <- plot_vb_diagnostics(vb_fit)
  print(p_diag)
}

# 10. Save all plots to a PDF
cat("\n10. Saving all plots to test_results.pdf\n")
pdf("test_results.pdf")
print(p_counts)
print(p_counts_year)
print(p_counts_sex)
print(p_counts_both)
print(p_heatmap)
if(!inherits(vb_fit, "try-error")) {
  print(p_fit)
  print(p_diag)
}
dev.off()

cat("\nTest complete! Results saved to test_results.pdf\n")