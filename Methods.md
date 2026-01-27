# Methods description for growthVB

This document provides both the mathematical framework underlying the `growthVB` package and ready-to-use methods descriptions for scientific reports and publications.

---

# Part A: Mathematical Framework

## 1. Introduction

The `growthVB` package implements methods for two fundamental relationships in fisheries biology:

1. **Length-weight relationships**: Allometric models relating fish length to weight
2. **Von Bertalanffy growth function (VBGF)**: Models describing fish growth over time

This section outlines the mathematical framework for both frequentist and Bayesian implementations.

---

## 2. Length-Weight Relationships

### 2.1 The Allometric Model

The relationship between fish length ($L$) and weight ($W$) follows the standard allometric equation:

$$
W = aL^b
$$

where:

- $a$ is the scaling coefficient (condition factor)
- $b$ is the allometric exponent

For isometric growth (fish maintain the same shape as they grow), $b = 3$. Departures from $b = 3$ indicate allometric growth:

- $b < 3$: Negative allometry (fish become relatively lighter/more elongated with size)
- $b > 3$: Positive allometry (fish become relatively heavier/more rotund with size)

### 2.2 Log-Linear Regression

The model is linearised by log-transformation:

$$
\log(W) = \log(a) + b \cdot \log(L) + \varepsilon
$$

where $\varepsilon \sim N(0, \sigma^2)$ represents random error. This is fitted using ordinary least squares regression.

### 2.3 Bias Correction for Back-Transformation

Back-transformation of the intercept from log scale introduces systematic negative bias because for a random variable $X$:

$$
E[\exp(X)] = \exp\left(\mu + \frac{\sigma^2}{2}\right) \neq \exp(\mu)
$$

Following Sprugel (1983) and Miller (1984), the bias-corrected estimate of $a$ is:

$$
\hat{a} = \exp\left(\hat{\alpha} + \frac{\hat{\sigma}^2}{2}\right)
$$

where:

- $\hat{\alpha}$ is the estimated intercept on the log scale
- $\hat{\sigma}^2$ is the residual variance (mean squared error) from the regression

### 2.4 Isometry Testing

Isometric growth is tested using a t-test:

$$
t = \frac{\hat{b} - 3}{SE(\hat{b})}
$$

where $\hat{b}$ is the estimated allometric exponent and $SE(\hat{b})$ is its standard error. The test statistic follows a t-distribution with $n-2$ degrees of freedom under the null hypothesis $H_0: b = 3$.

---

## 3. The Von Bertalanffy Growth Function

### 3.1 Basic Equation

The standard von Bertalanffy growth function is defined as:

$$
L_t = L_\infty \cdot (1 - e^{-k(t-t_0)})
$$

where:

- $L_t$ is the length at age $t$
- $L_\infty$ is the asymptotic length (theoretical maximum length)
- $k$ is the growth coefficient (rate at which $L_\infty$ is approached, year⁻¹)
- $t_0$ is the theoretical age at which length would be zero (years)

### 3.2 Biological Interpretation

- **$L_\infty$**: Represents the asymptotic average length of old fish in the population. It is not the maximum observed length, but the mean length that fish would approach if they lived indefinitely.
- **$k$**: Reflects how quickly fish approach their maximum size. Higher values indicate faster growth towards $L_\infty$.
- **$t_0$**: Theoretical age at zero length. Typically negative and has limited biological meaning; it is a mathematical artefact of the model.

### 3.3 Heteroscedastic Error Structure

In fish growth data, variance often increases with fish size. We model this heteroscedasticity by specifying that the standard deviation is proportional to the expected length:

$$
\sigma_i = CV \cdot \hat{L}_i
$$

where:

- $\sigma_i$ is the standard deviation for observation $i$
- $\hat{L}_i$ is the predicted length for observation $i$
- $CV$ is the coefficient of variation (constant)

The observation model becomes:

$$
L_i \sim N\left(L_\infty \cdot (1 - e^{-k(t_i-t_0)}), (CV \cdot \hat{L}_i)^2\right)
$$

---

## 4. Frequentist Implementation (Maximum Likelihood)

### 4.1 Likelihood Function

The negative log-likelihood function is:

$$
-\ln(\mathcal{L}) = \sum_{i=1}^{n} \left[ \ln(\sigma_i) + \frac{1}{2} \ln(2\pi) + \frac{(L_i - \hat{L}_i)^2}{2\sigma_i^2} \right]
$$

Substituting $\sigma_i = CV \cdot \hat{L}_i$:

$$
-\ln(\mathcal{L}) = \sum_{i=1}^{n} \left[ \ln(CV \cdot \hat{L}_i) + \frac{1}{2} \ln(2\pi) + \frac{(L_i - \hat{L}_i)^2}{2(CV \cdot \hat{L}_i)^2} \right]
$$

### 4.2 Optimisation

Parameters ($L_\infty$, $k$, $t_0$, $CV$) are estimated by minimising the negative log-likelihood using the L-BFGS-B algorithm (Byrd et al. 1995) with bounded constraints to ensure biologically realistic values.

### 4.3 Parameter Uncertainty

The variance-covariance matrix is estimated from the inverse Hessian:

$$
\text{Var-Cov} \approx H^{-1}
$$

Standard errors are the square roots of the diagonal elements:

$$
SE(\hat{\theta}_j) = \sqrt{(H^{-1})_{jj}}
$$

### 4.4 Confidence Intervals

Asymptotic 95% confidence intervals:

$$
\hat{\theta}_j \pm 1.96 \times SE(\hat{\theta}_j)
$$

For derived quantities (e.g., predicted lengths), uncertainty is propagated using the delta method.

---

## 5. Bayesian Implementation (brms/Stan)

### 5.1 Model Formulation

$$
L_i \sim N(\mu_i, \sigma_i^2)
$$

$$
\mu_i = L_\infty \cdot (1 - e^{-k(t_i-t_0)})
$$

$$
\sigma_i = CV \cdot \mu_i
$$

### 5.2 Prior Distributions

$$
L_\infty \sim N(\mu_{L_\infty}, \sigma_{L_\infty}^2)
$$

$$
k \sim N^+(\mu_k, \sigma_k^2) \quad \text{(truncated at 0)}
$$

$$
t_0 \sim N(\mu_{t_0}, \sigma_{t_0}^2)
$$

$$
CV \sim N^+(\mu_{CV}, \sigma_{CV}^2) \quad \text{(truncated at 0)}
$$

Hyperparameters are chosen based on biological knowledge or set to be weakly informative.

### 5.3 Posterior Estimation

Posterior distributions are estimated using Hamiltonian Monte Carlo via the No-U-Turn Sampler (NUTS; Hoffman & Gelman 2014) implemented in Stan (Carpenter et al. 2017):

$$
p(L_\infty, k, t_0, CV | \text{data}) \propto p(\text{data}|L_\infty, k, t_0, CV) \cdot p(L_\infty) \cdot p(k) \cdot p(t_0) \cdot p(CV)
$$

### 5.4 Convergence Diagnostics

Convergence is assessed using:

- $\hat{R}$ statistic (values < 1.01 indicate convergence; Vehtari et al. 2021)
- Effective sample size (ESS)
- Visual inspection of trace plots

### 5.5 Credible Intervals

Bayesian 95% credible intervals are the 2.5th and 97.5th percentiles of the posterior distributions.

---

## 6. Sex-Specific Growth Models

### 6.1 Separate Models

For sex-specific growth, separate VBGF parameters are estimated:

**Males**: $L_t^M = L_\infty^M \cdot (1 - e^{-k^M(t-t_0^M)})$

**Females**: $L_t^F = L_\infty^F \cdot (1 - e^{-k^F(t-t_0^F)})$

### 6.2 Parameter Set Expansion

For sex models, the parameter set expands from:

- Single model: $\{L_\infty, k, t_0, CV\}$
- Sex model: $\{L_\infty^M, L_\infty^F, k^M, k^F, t_0^M, t_0^F, CV^M, CV^F\}$

---

## 7. Bootstrap Permutation Testing for Group Comparisons

### 7.1 Test Framework

The `compare_vb_mle()` function implements bootstrap permutation testing to compare growth parameters between groups under the null hypothesis of no difference.

**Procedure:**

1. Fit VBGF models to each group separately
2. Calculate observed parameter differences (maximum pairwise difference)
3. Permute group labels $B$ times (typically 1000)
4. Refit models to permuted data
5. Calculate p-value as proportion of permuted differences ≥ observed

### 7.2 Age-Stratified Permutation

To ensure valid testing when age distributions differ between groups, permutations are performed within age strata:

$$
\text{Age bins} = \{\text{bin}_1, \text{bin}_2, \ldots, \text{bin}_m\}
$$

Each bin spans a fixed width (default 2 years). Group labels are permuted only within each bin, preserving age distributions.

### 7.3 Growth Curve Comparison

Overall trajectory differences are assessed using maximum absolute deviance:

$$
D_{max} = \max_{t \in [t_{min}, t_{max}]} |L_A(t) - L_B(t)|
$$

This captures the largest difference between curves regardless of which parameter drives it.

### 7.4 Multiple Group Comparisons

For $G$ groups, the test statistic for parameter $\theta$ is:

$$
T_\theta = \max_{i,j \in \{1,...,G\}} |\hat{\theta}_i - \hat{\theta}_j|
$$

---

## 8. Model Comparison and Selection

### 8.1 Information Criteria

- **AIC**: $AIC = -2\ln(L) + 2p$
- **BIC**: $BIC = -2\ln(L) + p\ln(n)$

where $L$ is the likelihood, $p$ is the number of parameters, and $n$ is the sample size.

### 8.2 Bayesian Model Comparison

- **WAIC**: Widely Applicable Information Criterion
- **LOO-CV**: Leave-One-Out Cross-Validation

---

## 9. Residuals and Model Diagnostics

### 9.1 MLE Residuals

$$
r_i = L_i - \hat{L}_i
$$

### 9.2 Bayesian Residuals

Using posterior mean fitted values:

$$
r_i = L_i - E[\hat{L}_i | \text{data}]
$$

where:

$$
E[\hat{L}_i | \text{data}] \approx \frac{1}{S} \sum_{s=1}^{S} L_\infty^{(s)} \cdot (1 - e^{-k^{(s)}(t_i-t_0^{(s)})})
$$

### 9.3 Posterior Predictive Checks

For Bayesian models:

- Density overlays comparing observed vs simulated data
- Q-Q plots comparing quantiles
- Residual plots examining systematic deviations

---

# Part B: Report Templates

This section provides ready-to-use methods descriptions for scientific reports and publications.

---

## 10. Length-Weight Relationship Templates

### 10.1 Brief Version

> Length-weight relationships were estimated using the allometric model $W = aL^b$, fitted via linear regression on log-transformed data. A bias correction factor was applied when back-transforming the intercept to the original scale (Sprugel 1983). Isometric growth was tested using a t-test (H₀: b = 3). All analyses were conducted using the `growthVB` package (Dunn 2025) in R (R Core Team 2024).

### 10.2 Detailed Version

> The relationship between fish length ($L$) and weight ($W$) was modelled using the standard allometric equation $W = aL^b$, where $a$ is the scaling coefficient and $b$ is the allometric exponent. This model was linearised by log-transformation and fitted using ordinary least squares regression.
>
> Back-transformation of the intercept from log scale introduces a systematic negative bias because $E[\exp(X)] \neq \exp(E[X])$ for random variables. Following Sprugel (1983) and Miller (1984), a bias correction was applied: $\hat{a} = \exp(\hat{\alpha} + \hat{\sigma}^2/2)$, where $\hat{\alpha}$ is the estimated intercept on the log scale and $\hat{\sigma}^2$ is the residual variance from the regression.
>
> Isometric growth (b = 3) was tested using a t-test. Significant departure from b = 3 indicates allometric growth: negative allometry (b < 3) implies fish become relatively lighter as they grow, while positive allometry (b > 3) implies they become relatively heavier.
>
> All length-weight analyses were conducted using the `fit_lw()` function in the `growthVB` package version X.X.X (Dunn 2025) in R version X.X.X (R Core Team 2024).

---

## 11. Von Bertalanffy Growth Curve Templates

### 11.1 Brief Version (MLE)

> Growth was modelled using the von Bertalanffy growth function (VBGF), $L_t = L_\infty(1 - e^{-k(t-t_0)})$, fitted by maximum likelihood estimation assuming heteroscedastic errors where the coefficient of variation (CV) is constant across sizes. Parameter estimates and 95% confidence intervals were obtained using the delta method. All analyses were conducted using the `growthVB` package (Dunn 2025) in R (R Core Team 2024).

### 11.2 Brief Version (Bayesian)

> Growth was modelled using a Bayesian implementation of the von Bertalanffy growth function (VBGF). The model incorporated heteroscedastic errors with CV modelled as a function of predicted length. Posterior distributions were estimated using Hamiltonian Monte Carlo via Stan (Carpenter et al. 2017), with X chains of Y iterations (Z warmup). Convergence was assessed using $\hat{R}$ statistics and visual inspection of trace plots. All analyses were conducted using the `growthVB` package (Dunn 2025), which interfaces with `brms` (Bürkner 2017).

### 11.3 Detailed Version

> Individual growth was modelled using the von Bertalanffy growth function (VBGF; von Bertalanffy 1938):
>
> $$
> L_t = L_\infty (1 - e^{-k(t - t_0)})
> $$
>
> where $L_t$ is the expected length at age $t$, $L_\infty$ is the asymptotic length, $k$ is the Brody growth coefficient (year⁻¹), and $t_0$ is the theoretical age at which length equals zero.
>
> The model assumed heteroscedastic normally distributed errors, where $\sigma_t = CV \cdot L_t$ and $CV$ is the coefficient of variation. This error structure accounts for the common observation that variability in fish length increases with size (Kimura 2008).
>
> **For MLE:** Parameters were estimated by maximising the log-likelihood using the L-BFGS-B algorithm with bounded parameter constraints. Standard errors were estimated from the Hessian matrix, and 95% confidence intervals were calculated using the delta method.
>
> **For Bayesian:** Posterior distributions were estimated using the No-U-Turn Sampler (NUTS; Hoffman & Gelman 2014) implemented via Stan. Prior distributions were: [specify priors]. X Markov chains were run for Y iterations each, with Z warmup iterations discarded. Convergence was assessed using $\hat{R}$ < 1.01 (Vehtari et al. 2021).
>
> Growth analyses were conducted using the `fit_vb_mle()` [or `fit_vb_brms()`] function in the `growthVB` package version X.X.X (Dunn 2025).

---

## 12. Group Comparison Templates

### 12.1 Brief Version

> Differences in growth parameters between [groups/populations/sexes] were tested using bootstrap permutation tests (n = X permutations). Age-stratified permutation was used to maintain the age structure within groups. Significance was assessed at α = 0.05. Additionally, overall growth curve differences were evaluated using the maximum absolute deviance between fitted curves.

### 12.2 Detailed Version

> To test for differences in VBGF parameters between groups, a bootstrap permutation approach was employed. Under the null hypothesis of no difference between groups, group labels were randomly permuted and VBGF models were refitted to generate null distributions of parameter differences.
>
> For each of X permutations: (1) group labels were randomly reassigned maintaining age-stratified structure, (2) VBGF parameters were estimated for each permuted group, and (3) the difference in each parameter between groups was calculated. Two-sided p-values were calculated as the proportion of permuted differences with absolute values exceeding the observed absolute difference.
>
> To maintain valid age distributions during permutation, data were stratified by age bins of Y years. Within each stratum, group labels were permuted independently, ensuring the age structure of each group was preserved.
>
> Overall differences in growth trajectories were assessed by calculating the maximum absolute deviance: $D_{max} = \max_t |L_A(t) - L_B(t)|$, which captures the largest difference between growth curves regardless of the specific parameter driving the difference.
>
> Group comparisons were performed using the `compare_vb_mle()` function in the `growthVB` package version X.X.X (Dunn 2025).

---

## 13. Results Table Templates

### 13.1 Length-Weight Results

| Parameter | Estimate | SE    | 95% CI         |
| --------- | -------- | ----- | -------------- |
| a         | X.XXe-X  | -     | -              |
| b         | X.XXX    | X.XXX | X.XXX – X.XXX |
| n         | XXX      | -     | -              |
| R²       | X.XXX    | -     | -              |

Isometry test: t = X.XX, df = XX, p = X.XXX [significant/not significant departure from isometry]

### 13.2 Growth Parameter Results

| Sex      | Parameter         | Estimate | SE    | 95% CI         |
| -------- | ----------------- | -------- | ----- | -------------- |
| Combined | $L_\infty$ (cm) | XXX.X    | X.X   | XXX.X – XXX.X |
| Combined | $k$ (year⁻¹)  | X.XXX    | X.XXX | X.XXX – X.XXX |
| Combined | $t_0$ (years)   | -X.XX    | X.XX  | -X.XX – -X.XX |
| Combined | CV                | X.XXX    | X.XXX | X.XXX – X.XXX |

### 13.3 Group Comparison Results

| Parameter        | Group A | Group B | Difference | p-value |
| ---------------- | ------- | ------- | ---------- | ------- |
| $L_\infty$     | XXX.X   | XXX.X   | XX.X       | X.XXX   |
| $k$            | X.XXX   | X.XXX   | X.XXX      | X.XXX   |
| $t_0$          | -X.XX   | -X.XX   | X.XX       | X.XXX   |
| Curve$D_{max}$ | -       | -       | XX.X cm    | X.XXX   |

---

## 14. References

Bürkner, P. C. (2017). brms: An R package for Bayesian multilevel models using Stan. *Journal of Statistical Software*, 80(1), 1-28.

Byrd, R. H., Lu, P., Nocedal, J., & Zhu, C. (1995). A limited memory algorithm for bound constrained optimization. *SIAM Journal on Scientific Computing*, 16(5), 1190-1208.

Carpenter, B., Gelman, A., Hoffman, M. D., Lee, D., Goodrich, B., Betancourt, M., ... & Riddell, A. (2017). Stan: A probabilistic programming language. *Journal of Statistical Software*, 76(1).

Dunn, A. (2025). growthVB: Von Bertalanffy Growth Curve Estimation with Heteroscedastic Errors. R package. https://github.com/alistairdunn1/growthVB

Froese, R. (2006). Cube law, condition factor and weight-length relationships: history, meta-analysis and recommendations. *Journal of Applied Ichthyology*, 22(4), 241-253.

Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., & Rubin, D. B. (2013). *Bayesian Data Analysis* (3rd ed.). CRC Press.

Hoffman, M. D., & Gelman, A. (2014). The No-U-Turn sampler: adaptively setting path lengths in Hamiltonian Monte Carlo. *Journal of Machine Learning Research*, 15(1), 1593-1623.

Katsanevakis, S. (2006). Modelling fish growth: Model selection, multi-model inference and model selection uncertainty. *Fisheries Research*, 81(2-3), 229-235.

Kimura, D. K. (2008). Extending the von Bertalanffy growth model using explanatory variables. *Canadian Journal of Fisheries and Aquatic Sciences*, 65(9), 1879-1891.

Miller, D. M. (1984). Reducing transformation bias in curve fitting. *The American Statistician*, 38(2), 124-126.

R Core Team (2024). R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing, Vienna, Austria. https://www.R-project.org/

Sprugel, D. G. (1983). Correcting for bias in log-transformed allometric equations. *Ecology*, 64(1), 209-210.

Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Bürkner, P. C. (2021). Rank-normalization, folding, and localization: An improved $\hat{R}$ for assessing convergence of MCMC. *Bayesian Analysis*, 16(2), 667-718.

von Bertalanffy, L. (1938). A quantitative theory of organic growth. *Human Biology*, 10(2), 181-213.

---

## 15. Usage Notes

1. Replace X.X.X with actual version numbers and parameter values
2. Adjust prior specifications based on your specific analysis
3. Include/exclude sex-specific sections as appropriate
4. Modify confidence interval levels if different from 95%
5. Add species-specific context and biological interpretation
6. Cite the package and relevant methodological references
