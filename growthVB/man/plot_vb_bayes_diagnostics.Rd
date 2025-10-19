\name{plot_vb_bayes_diagnostics}
\alias{plot_vb_bayes_diagnostics}
\title{Comprehensive Bayesian Diagnostics for von Bertalanffy Models}
\usage{
plot_vb_bayes_diagnostics(
  model,
  ndraws = 50,
  diagnostic_types = "all",
  pars = NULL
)
}
\arguments{
\item{model}{A model object returned by fit_vb_brms()}

\item{ndraws}{Number of posterior draws to include in plots (default 50)}

\item{diagnostic_types}{Character vector specifying which diagnostics to include.
Options: "convergence", "pp_checks", "residuals", "parameters", "growth_specific", "all"}

\item{pars}{Optional character vector of parameter names to focus on}
}
\value{
A list of ggplot2 objects organised by diagnostic type
}
\description{
This function creates comprehensive MCMC and posterior predictive diagnostic plots 
for Bayesian von Bertalanffy growth models, including convergence diagnostics,
enhanced posterior predictive checks, and model-specific growth diagnostics.
}
\examples{
\dontrun{
# Simple example with simulated data
age <- 1:15
length <- 100 * (1 - exp(-0.2 * (age - (-0.5)))) + rnorm(15, 0, 5)
fit <- fit_vb_brms(age = age, length = length)

# Full diagnostics
diagnostics <- plot_vb_bayes_diagnostics(fit)

# Specific diagnostic types
convergence <- plot_vb_bayes_diagnostics(fit, diagnostic_types = "convergence")
pp_checks <- plot_vb_bayes_diagnostics(fit, diagnostic_types = "pp_checks")
}
}
