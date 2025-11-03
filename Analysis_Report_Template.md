# Von Bertalanffy Growth Analysis Report Template

## 1. Executive Summary

*[Brief overview of the analysis, key findings, and conclusions in 3-5 sentences]*

**Key Results:**
- Asymptotic length (L∞): [value] ± [error] [units]
- Growth coefficient (k): [value] ± [error] [year⁻¹]
- Age at zero length (t₀): [value] ± [error] [years]
- *[Any other key findings]*

## 2. Introduction

### 2.1 Study Background
*[Provide context for why this analysis was conducted]*

### 2.2 Objectives
*[List specific objectives of the growth analysis]*

- To estimate growth parameters for [species name]
- To compare growth between [different groups, if applicable]
- To examine [any specific hypotheses]
- *[Other objectives]*

### 2.3 Species Biology
*[Brief description of the species and relevant biological information]*

## 3. Materials and Methods

### 3.1 Data Collection

**Sampling period:** [dates]

**Geographic area:** [location]

**Sampling methods:** *[Brief description of how samples were collected]*

**Sample size:**
- Total: [n] individuals
- Males: [n] individuals
- Females: [n] individuals
- Unknown sex: [n] individuals
- Age range: [min] to [max] years
- Length range: [min] to [max] [units]

### 3.2 Age and Length Determination

**Length measurements:**
- Type: [total length/fork length/standard length]
- Precision: [value] [units]
- Method: [measuring board/calipers/etc.]

**Age determination:**
- Structure used: [otoliths/scales/spines/vertebrae]
- Preparation method: *[Brief description]*
- Reading method: *[Brief description]*
- Quality control: *[Describe any age validation procedures]*

### 3.3 Statistical Analysis

**Software used:**
- R version [x.x.x]
- growthVB package version [x.x.x]
- *[Other relevant software]*

**Growth models applied:**
- Von Bertalanffy growth function (VBGF)
- *[Any modified models or alternative models]*

**Estimation methods:**
- Maximum likelihood estimation (MLE) with heteroscedasticity
- Bayesian estimation using BRMS
- *[Other methods if applicable]*

**Data traceability:**
- ID indexing preserves original input order for data matching
- Fitted values and residuals calculated for all observations
- Full result traceability maintained throughout analysis

**Model assumptions:**
- *[List any assumptions made in the analysis]*

**Model selection criteria:**
- AIC/BIC/WAIC
- Likelihood ratio tests
- *[Other criteria used]*

## 4. Results

### 4.1 Data Overview

```
[Insert table summarising sample characteristics]
```

#### 4.1.1 Age Distribution

*[Insert age frequency plot using plot_vb_age_counts(age)]*

```r
# Basic age distribution plot
plot_vb_age_counts(age = age_data)
```

#### 4.1.2 Age Distribution by Sex

*[Insert age frequency by sex plot using plot_vb_age_counts(age, sex=sex_data)]*

```r
# Age distribution by sex
plot_vb_age_counts(age = age_data, sex = sex_data)
```

#### 4.1.3 Age Distribution by Group

*[Insert age frequency by group plot using plot_vb_age_counts(age, group=group_data)]*

```r
# Age distribution by group (e.g., year, study, site)
plot_vb_age_counts(age = age_data, group = group_data)
```

*[Additional visualisations as needed]*

#### 4.1.4 Age-Length Relationship Overview

*[Insert age-length heatmap showing sample distribution]*

```r
# Age-length distribution heatmap with smoother
plot_age_length_heatmap(age = age_data, length = length_data, add_smoother = TRUE)
```

#### 4.1.5 Empirical Coefficient of Variation

*[Insert empirical CV plot to examine variance patterns]*

```r
# Empirical CV by age
plot_empirical_cv(age = age_data, length = length_data)

# CV by age and sex
plot_empirical_cv(age = age_data, length = length_data, sex = sex_data)
```

### 4.2 Growth Parameter Estimates

**Data structure note:** Both `fit_vb_mle()` and `fit_vb_brms()` return enhanced data with:
- `ID`: Unique identifier preserving original input order
- `fitted`: Model-predicted length values  
- `residual`: Residuals (observed - fitted values)

```r
# Access model data with fitted values and residuals
head(fitted_model$data)

# Check original order preservation
all(fitted_model$data$ID == seq_along(original_data))

# Summary of model residuals
summary(fitted_model$data$residual)
```

#### 4.2.1 Combined Sex Analysis

**Maximum likelihood estimates:**

| Parameter | Estimate | Standard Error | 95% CI Lower | 95% CI Upper |
| --------- | -------- | -------------- | ------------ | ------------ |
| L∞        | [value]  | [value]        | [value]      | [value]      |
| k         | [value]  | [value]        | [value]      | [value]      |
| t₀        | [value]  | [value]        | [value]      | [value]      |
| CV        | [value]  | [value]        | [value]      | [value]      |

**Bayesian estimates:**

| Parameter | Posterior Mean | Posterior SD | 95% CI Lower | 95% CI Upper |
| --------- | -------------- | ------------ | ------------ | ------------ |
| L∞        | [value]        | [value]      | [value]      | [value]      |
| k         | [value]        | [value]      | [value]      | [value]      |
| t₀        | [value]        | [value]      | [value]      | [value]      |
| CV        | [value]        | [value]      | [value]      | [value]      |

*[Insert von Bertalanffy growth curve plot showing fitted curve and observed data points]*

#### 4.2.2 Sex-Specific Analysis

**Male growth parameters:**

| Parameter | Estimate | Standard Error | 95% CI Lower | 95% CI Upper |
| --------- | -------- | -------------- | ------------ | ------------ |
| L∞        | [value]  | [value]        | [value]      | [value]      |
| k         | [value]  | [value]        | [value]      | [value]      |
| t₀        | [value]  | [value]        | [value]      | [value]      |
| CV        | [value]  | [value]        | [value]      | [value]      |

**Female growth parameters:**

| Parameter | Estimate | Standard Error | 95% CI Lower | 95% CI Upper |
| --------- | -------- | -------------- | ------------ | ------------ |
| L∞        | [value]  | [value]        | [value]      | [value]      |
| k         | [value]  | [value]        | [value]      | [value]      |
| t₀        | [value]  | [value]        | [value]      | [value]      |
| CV        | [value]  | [value]        | [value]      | [value]      |

*[Insert sex-specific growth curve plot with different colours/symbols for each sex]*

### 4.3 Model Diagnostics

#### 4.3.1 Residual Analysis

*[Insert comprehensive diagnostic plots]*

```r
# Standard diagnostic plots
diagnostics <- plot_vb_mle_diagnostics(fitted_model)

# View individual diagnostic plots
diagnostics$residuals_vs_fitted  # Residuals vs fitted values
diagnostics$qq_plot              # Q-Q plot of residuals  
diagnostics$scale_location       # Scale-location plot
diagnostics$residuals_vs_age     # Residuals vs age
```

*[Describe any patterns or concerns in the residuals]*

#### 4.3.2 Model Comparison

| Model     | Parameters | Log-likelihood | AIC     | BIC     | ΔAIC    |
| --------- | ---------- | -------------- | ------- | ------- | ------- |
| [Model 1] | [n]        | [value]        | [value] | [value] | [value] |
| [Model 2] | [n]        | [value]        | [value] | [value] | [value] |
| [Model 3] | [n]        | [value]        | [value] | [value] | [value] |

#### 4.3.3 Group Comparisons

*[Include this section if comparing growth between multiple groups/populations]*

**Bootstrap permutation test approach:**
- Statistical method: Bootstrap permutation testing
- Permutations: [n] iterations  
- Age-stratified: [Yes/No] with [n]-year bins
- Parameters tested: [list parameters]
- Growth curve testing: [Yes/No]

```r
# Compare growth parameters between groups
group_comparison <- compare_vb_mle(
  age = age_data,
  length = length_data, 
  group = group_data,
  n_bootstrap = 1000,
  min_obs = 50,          # Minimum observations per group
  age_stratified = TRUE, # Preserve age distributions
  age_bin_width = 2,     # 2-year age bins
  test_curves = TRUE,    # Test overall growth curves
  seed = 123
)

# Professional formatted summary
summary(group_comparison)

# Visualise permutation test results
plot(group_comparison)

# Plot specific parameters
plot(group_comparison, parameters = c("Linf", "k"))
```

**Overall test results:**
- Total observations: [n]
- Groups compared: [list groups]
- Group sizes: [list sizes]
- Significant parameters: [n] of [total]
- Growth curve difference: [Significant/Not significant]

**Parameter comparison table:**

| Parameter | Group A | Group B | Observed Difference | P-value | Significant |
| --------- | ------- | ------- | ------------------ | ------- | ----------- |
| L∞        | [value] | [value] | [value]            | [value] | [Yes/No]    |
| k         | [value] | [value] | [value]            | [value] | [Yes/No]    |
| t₀        | [value] | [value] | [value]            | [value] | [Yes/No]    |
| CV        | [value] | [value] | [value]            | [value] | [Yes/No]    |

**Growth curve comparison:**
- Maximum curve deviance: [value]
- Curve p-value: [value]  
- Interpretation: [Growth trajectories are/are not significantly different]

*[Discuss the biological significance of any detected differences]*

**Sex-specific group comparisons:**

*[If using sex-specific models, include results for sex-specific parameters]*

| Parameter | Group A | Group B | Difference | P-value | Significant |
| --------- | ------- | ------- | ---------- | ------- | ----------- |
| L∞_M      | [value] | [value] | [value]    | [value] | [Yes/No]    |
| L∞_F      | [value] | [value] | [value]    | [value] | [Yes/No]    |
| k_M       | [value] | [value] | [value]    | [value] | [Yes/No]    |
| k_F       | [value] | [value] | [value]    | [value] | [Yes/No]    |

#### 4.3.4 Bayesian Diagnostic Plots

*[For Bayesian analyses only - Insert convergence diagnostics and posterior distribution plots]*

*[Describe any issues with convergence or posterior distributions]*

## 5. Discussion

### 5.1 Parameter Interpretation

*[Discuss the biological meaning of the estimated parameters]*

### 5.2 Comparison with Previous Studies

| Study           | Region     | Year   | L∞      | k       | t₀      | Method   |
| --------------- | ---------- | ------ | ------- | ------- | ------- | -------- |
| Current study   | [location] | [year] | [value] | [value] | [value] | [method] |
| [Author et al.] | [location] | [year] | [value] | [value] | [value] | [method] |
| [Author et al.] | [location] | [year] | [value] | [value] | [value] | [method] |

*[Discuss similarities and differences with previous studies]*

### 5.3 Sex-Specific Growth Differences

*[If applicable, discuss differences in growth parameters between sexes and their biological significance]*

### 5.4 Limitations

*[Discuss any limitations of the data or analytical approach]*

## 6. Conclusions

*[Summarise the main findings and their implications for fisheries management, stock assessment, or other applications]*

## 7. References

*[List all references in appropriate citation format]*

## Appendix A: Supplementary Figures and Tables

*[Include any additional figures or tables that support the analysis but aren't essential to the main narrative]*

## Appendix B: R Code

```r
# Example code used for the analysis
library(growthVB)

# Data preparation
# [code snippet]

# Fitting the model
# [code snippet]

# Plotting results
# [code snippet]
```

## Appendix C: Raw Data

*[Include a link to or description of where the raw data can be accessed, if applicable]*
