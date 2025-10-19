# Von Bertalanffy Growth Model: Mathematical Methods

## 1. Introduction

The von Bertalanffy growth function (VBGF) is one of the most widely used models in fisheries science to describe the growth of individual fish. This document outlines the mathematical framework for both frequentist and Bayesian implementations of the VBGF as implemented in the `growthVB` package.

## 2. The Von Bertalanffy Growth Function

### 2.1 Basic Equation

The standard von Bertalanffy growth function is defined as:

$$L_t = L_\infty \cdot (1 - e^{-k(t-t_0)})$$

Where:
- $L_t$ is the length at age $t$
- $L_\infty$ is the asymptotic length (theoretical maximum length)
- $k$ is the growth coefficient (rate at which $L_\infty$ is approached)
- $t_0$ is the theoretical age at which length would be zero

### 2.2 Biological Interpretation

- $L_\infty$: Represents the asymptotic average length of very old fish in the population
- $k$: Reflects how quickly the fish approaches its maximum size (higher values indicate faster growth)
- $t_0$: Theoretical age at zero length (typically negative and has limited biological meaning)

## 3. Frequentist Implementation (Maximum Likelihood Method)

### 3.1 Model Formulation

In the frequentist approach, we estimate the VBGF parameters using maximum likelihood estimation. We assume that observed lengths follow the model:

$$L_i = L_\infty \cdot (1 - e^{-k(t_i-t_0)}) + \varepsilon_i$$

Where $\varepsilon_i \sim \mathcal{N}(0, \sigma_i^2)$ represents random error with heteroscedastic variance.

### 3.2 Heteroscedasticity Modeling

In fish growth data, variance often increases with fish size. We explicitly model this heteroscedasticity by specifying that the standard deviation is proportional to the expected length:

$$\sigma_i = \text{CV} \cdot \hat{L}_i$$

Where:
- $\sigma_i$ is the standard deviation for observation $i$
- $\hat{L}_i$ is the predicted length for observation $i$
- $\text{CV}$ is the coefficient of variation

This means the observation model becomes:

$$L_i \sim \mathcal{N}(L_\infty \cdot (1 - e^{-k(t_i-t_0)}), (\text{CV} \cdot \hat{L}_i)^2)$$

### 3.3 Maximum Likelihood Estimation

We use maximum likelihood estimation to find the parameters ($L_\infty$, $k$, $t_0$, and $\text{CV}$) that best fit the data. The negative log-likelihood function is:

$$-\ln(\mathcal{L}) = \sum_{i=1}^{n} \left[ \ln(\sigma_i) + \frac{1}{2} \ln(2\pi) + \frac{(L_i - \hat{L}_i)^2}{2\sigma_i^2} \right]$$

Substituting $\sigma_i = \text{CV} \cdot \hat{L}_i$:

$$-\ln(\mathcal{L}) = \sum_{i=1}^{n} \left[ \ln(\text{CV} \cdot \hat{L}_i) + \frac{1}{2} \ln(2\pi) + \frac{(L_i - \hat{L}_i)^2}{2(\text{CV} \cdot \hat{L}_i)^2} \right]$$

We minimize this function using the `optim()` function in R, with the BFGS algorithm as the default optimization method.

### 3.4 Parameter Uncertainty

We estimate the variance-covariance matrix of the parameters using the inverse of the Hessian matrix from the optimization procedure:

$$\text{Var-Cov} \approx H^{-1}$$

The standard errors of the parameters are the square roots of the diagonal elements of this matrix:

$$\text{SE}(\hat{\theta}_j) = \sqrt{(H^{-1})_{jj}}$$

### 3.5 Confidence Intervals

We compute 95% confidence intervals for the parameters using the asymptotic normality of maximum likelihood estimators:

$$\hat{\theta}_j \pm 1.96 \times \text{SE}(\hat{\theta}_j)$$

For derived quantities like predicted lengths at given ages, we use the delta method to propagate uncertainty.

### 3.3 Parameter Estimation via Maximum Likelihood

The likelihood function for this model is:

$$L(L_\infty, k, t_0, \tau | \text{data}) = \prod_{i=1}^{n} \frac{1}{\tau \cdot \hat{L}_i \cdot \sqrt{2\pi}} \exp\left(-\frac{(L_i - \hat{L}_i)^2}{2 \cdot (\tau \cdot \hat{L}_i)^2}\right)$$

Where $\hat{L}_i = L_\infty \cdot (1 - e^{-k(t_i-t_0)})$

We maximize the log-likelihood:

$$\ell(L_\infty, k, t_0, \tau | \text{data}) = -\sum_{i=1}^{n} \left[\log(\tau \cdot \hat{L}_i) + \frac{(L_i - \hat{L}_i)^2}{2 \cdot (\tau \cdot \hat{L}_i)^2} + \frac{1}{2}\log(2\pi)\right]$$

This is accomplished using numerical optimization methods such as BFGS or Nelder-Mead.

### 3.4 Confidence Intervals

For the maximum likelihood approach, confidence intervals are derived using:

1. **Hessian-based estimates**: The inverse of the negative Hessian matrix at the MLE provides an estimate of the variance-covariance matrix of parameters
2. **Delta method**: For derived quantities (like predicted lengths), we propagate uncertainty using the delta method
3. **Profile likelihood methods**: For more accurate intervals, profile likelihood can be used (computing likelihood-based confidence regions)

## 4. Bayesian Implementation (BRMS Method)

### 4.1 Model Formulation

In the Bayesian approach, we specify a full probability model:

$$L_i \sim \mathcal{N}(\mu_i, \sigma_i^2)$$
$$\mu_i = L_\infty \cdot (1 - e^{-k(t_i-t_0)})$$
$$\sigma_i = \tau \cdot \mu_i$$

### 4.2 Prior Distributions

We specify prior distributions for all parameters:

$$L_\infty \sim \mathcal{N}(\mu_{L_\infty}, \sigma_{L_\infty}^2)$$
$$k \sim \mathcal{N}^+(\mu_k, \sigma_k^2) \quad \text{(Half-normal, k > 0)}$$
$$t_0 \sim \mathcal{N}(\mu_{t_0}, \sigma_{t_0}^2)$$
$$\tau \sim \mathcal{N}^+(\mu_{\tau}, \sigma_{\tau}^2) \quad \text{(Half-normal, \tau > 0)}$$

Where the hyperparameters ($\mu_{L_\infty}$, $\sigma_{L_\infty}$, etc.) are chosen based on biological knowledge or set to be weakly informative.

### 4.3 Posterior Estimation

Using Markov Chain Monte Carlo (MCMC) sampling, we estimate the posterior distribution:

$$p(L_\infty, k, t_0, \tau | \text{data}) \propto p(\text{data}|L_\infty, k, t_0, \tau) \cdot p(L_\infty) \cdot p(k) \cdot p(t_0) \cdot p(\tau)$$

### 4.4 Credible Intervals

Bayesian credible intervals are derived directly from the posterior distributions of parameters, typically using the 2.5th and 97.5th percentiles for 95% credible intervals.

## 5. Sex-Specific Growth Models

### 5.1 Separate Models Approach

For sex-specific growth analyses, separate models are fitted for each sex:

**Males**: $L_t^M = L_\infty^M \cdot (1 - e^{-k^M(t-t_0^M)})$

**Females**: $L_t^F = L_\infty^F \cdot (1 - e^{-k^F(t-t_0^F)})$

### 5.2 Combined Model with Sex Factor

Alternatively, a combined model can be specified with sex as a factor:

$$L_t = L_\infty^{sex} \cdot (1 - e^{-k^{sex}(t-t_0^{sex})})$$

Where $sex \in \{M, F\}$ and parameters are estimated separately for each sex level.

## 6. Model Comparison and Selection

### 6.1 Information Criteria

Models can be compared using:

- **Akaike Information Criterion (AIC)**: $AIC = -2\ln(L) + 2p$
- **Bayesian Information Criterion (BIC)**: $BIC = -2\ln(L) + p\ln(n)$

Where $L$ is the likelihood, $p$ is the number of parameters, and $n$ is the sample size.

### 6.2 Likelihood Ratio Tests

For nested models, likelihood ratio tests can be used:

$$LR = -2\ln\left(\frac{L_{\text{reduced}}}{L_{\text{full}}}\right)$$

Under the null hypothesis, LR follows a chi-squared distribution with degrees of freedom equal to the difference in the number of parameters.

### 6.3 Bayesian Model Comparison

For Bayesian models, comparisons can use:

- **WAIC (Widely Applicable Information Criterion)**
- **LOO-CV (Leave-One-Out Cross-Validation)**
- **Bayes Factors**

## 7. Posterior Predictive Checks

### 7.1 Posterior Predictive Distribution

The posterior predictive distribution for a new observation $\tilde{L}$ at age $\tilde{t}$ is:

$$p(\tilde{L}|\tilde{t},\text{data}) = \int p(\tilde{L}|\tilde{t},\theta) \cdot p(\theta|\text{data}) \, d\theta$$

Where $\theta = (L_\infty, k, t_0, \tau)$ represents all model parameters.

### 7.2 Graphical Checks

Posterior predictive checks include:

1. **Density overlays**: Comparing observed data density to simulated posterior predictive draws
2. **Quantile-quantile plots**: Comparing observed data quantiles to posterior predictive quantiles
3. **Residual plots**: Examining systematic deviations between observed data and posterior predictions

## 8. Implementation in the growthVB Package

The `growthVB` package implements these methods with the following core functions:

1. `fit_vb_mle()`: Frequentist MLE estimation with optional heteroscedasticity modeling
2. `fit_vb_brms()`: Bayesian estimation using the brms package with Stan backend
3. `plot_vb_predictions()`: Visualizing model fit with confidence/credible intervals
4. `plot_vb_posteriors()`: Diagnostic plots for Bayesian models including posterior distributions and correlations

## 9. References

1. von Bertalanffy, L. (1938). A quantitative theory of organic growth. Human Biology, 10(2), 181-213.
2. Katsanevakis, S. (2006). Modelling fish growth: Model selection, multi-model inference and model selection uncertainty. Fisheries Research, 81(2-3), 229-235.
3. Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., & Rubin, D. B. (2013). Bayesian Data Analysis (3rd ed.). CRC Press.
4. Bürkner, P. C. (2017). brms: An R package for Bayesian multilevel models using Stan. Journal of Statistical Software, 80(1), 1-28.
