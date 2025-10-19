\name{plot_vb_growth_pp_checks}
\alias{plot_vb_growth_pp_checks}
\title{Growth-Specific Posterior Predictive Checks}
\usage{
plot_vb_growth_pp_checks(model, ndraws = 100, check_types = "all")
}
\arguments{
\item{model}{A model object returned by fit_vb_brms()}

\item{ndraws}{Number of posterior draws to include (default 100)}

\item{check_types}{Character vector specifying which checks to perform.
Options: "growth_patterns", "age_effects", "length_distribution", "residual_patterns", "all"}
}
\value{
A list of ggplot2 objects with growth-specific posterior predictive checks
}
\description{
This function creates posterior predictive checks specifically designed for 
von Bertalanffy growth models, focusing on growth-relevant patterns and diagnostics.
}
\examples{
\dontrun{
# Growth model with posterior predictive checks
age <- 1:15
length <- 100 * (1 - exp(-0.2 * (age - (-0.5)))) + rnorm(15, 0, 5)
fit <- fit_vb_brms(age = age, length = length)

# Growth-specific checks
growth_checks <- plot_vb_growth_pp_checks(fit)
}
}
