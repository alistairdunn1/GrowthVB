library(tidyverse)
library(brms)
library(loo)
library(bayesplot)

library(ggplot2)
library(dplyr)
require(viridis)
require(ggspatial)
require(raster)
require(lwgeom)
require(tidyr)

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
# Sys.setenv(LOCAL_CPPFLAGS = '-march=native')

load(file = make.filename("Age data/age_data.rdata", DIR$Data), verbose = TRUE)
data <- age.data
bad <- c(which(data$nage < 1))
#data <- data[-bad, ]
data <- data %>% filter(CCAMLR_year %in% 2023:2024)


data <- data %>% filter(!is.na(length), !is.na(age), !is.na(sex), sex %in% c("M", "F")) %>%
        mutate(age = floor(as.numeric(as.character(age))))
Table(data$age, data$sex)


#########################################################################################################
# brms von Bertalanffy
#########################################################################################################

age_F <- data %>%
  filter(sex == "F")

age_M <- data %>%
  filter(sex == "M")

# Run von Bertalanffy models in brms ----
# Strong k prior as the model for males has poor estimates
priors <- c(
  prior(normal(150, 100), nlpar = "Linf", lb = 0),
  prior(normal(0.0, 100), nlpar = "k", lb = 0),
  prior(normal(0, 1), nlpar = "t0"),
  prior(normal(0, 100), nlpar = "tau", lb = 0)
)

vb_F_1 <- brm(bf(length ~ eta, nl = TRUE) + nlf(eta ~ 1 + Linf * (1.0 - exp(-k * (age - t0)))) + nlf(sigma ~ eta * tau) + lf(Linf ~ 1, k ~ 1, t0 ~ 1, tau ~ 1),
  data = age_F, chains = 4, prior = priors, family = brmsfamily("gaussian", link_sigma = "identity"), iter = 4000
)
loo_F_1_vb <- loo(vb_F_1)

vb_M_1 <- brm(bf(length ~ eta, nl = TRUE) + nlf(eta ~ 1 + Linf * (1.0 - exp(-k * (age - t0)))) + nlf(sigma ~ eta * tau) + lf(Linf ~ 1, k ~ 1, t0 ~ 1, tau ~ 1),
  data = age_M, chains = 4, prior = priors, family = brmsfamily("gaussian", link_sigma = "identity"), iter = 4000
)
loo_M_1_vb <- loo(vb_M_1)

vb_M_1
vb_F_1

csv(rbind(fixef(vb_M_1),fixef(vb_F_1)), "AgeLength/vonB all ages (brms)", path = DIR$Tables)

dev.off()
new.graph(0.5)
mod <- vb_F_1
ce1 <- conditional_effects(mod, spaghetti = FALSE, method = "posterior_predict")
plot(ce1, ask = FALSE, points = TRUE)
mod <- vb_M_1
ce1 <- conditional_effects(mod, spaghetti = FALSE, method = "posterior_predict")
plot(ce1, ask = FALSE, points = TRUE)

df1 <- data.frame(age = sort(unique(vb_F_1$data$age)), Sex = "Female")
pp1 <- predict(vb_F_1, newdata = df1, probs = c(0.025, 0.975)) %>%
  data.frame() %>%
  bind_cols(df1) %>%
  mutate(Model = "von Bertalanffy")
df2 <- data.frame(age = sort(unique(vb_M_1$data$age)), Sex = "Male")
pp2 <- predict(vb_M_1, newdata = df2, probs = c(0.025, 0.975)) %>%
  data.frame() %>%
  bind_cols(df2) %>%
  mutate(Model = "von Bertalanffy")
df1 <- data.frame(vb_F_1$data, Sex = "Female", Model = "von Bertalanffy")
df2 <- data.frame(vb_M_1$data, Sex = "Male", Model = "von Bertalanffy")

dev.off()
new.graph(0.5)
p1 <- ggplot(data = bind_rows(pp1, pp2)) +
  geom_ribbon(aes(x = age, ymin = Q2.5, ymax = Q97.5, fill = Sex), alpha = 0.3) +
  geom_point(data = bind_rows(df1, df2), aes(x = age, y = length, colour = Sex), alpha = 0.25) +
  geom_line(aes(x = age, y = Estimate, colour = Sex), size = 1) +
  labs(x = "Age", y = "Length (cm)") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  st1 +
  theme(legend.position = c(0.8, 0.2))
print(p1)
SavePlot("AgeLength/brms_growth_vb", DIR$Figures)

dev.off()
new.graph(0.5)
yrep1 <- posterior_predict(vb_F_1, ndraws = 500)
p1 <- ppc_dens_overlay(y = vb_F_1$data[, 1], yrep = yrep1) +
  # coord_cartesian(xlim = range(exp(fit_mls2$data[,1]))) +
  labs(title = "(a) Female", x = "Length (cm)", y = "Density") +
  st1 +
  theme(legend.position = c(0.8, 0.2))
print(p1)
SavePlot("AgeLength/brms_growth_vb_F_ppdens", DIR$Figures)

yrep2 <- posterior_predict(vb_M_1, ndraws = 500)
p2 <- ppc_dens_overlay(y = vb_M_1$data[, 1], yrep = yrep2) +
  # coord_cartesian(xlim = range(exp(fit_mls2$data[,1]))) +
  labs(title = "(b) Male", x = "Length (cm)", y = "Density") +
  st1 +
  theme(legend.position = c(0.8, 0.2))
print(p2)
SavePlot("AgeLength/brms_growth_vb_M_ppdens", DIR$Figures)

dev.off()
new.graph(0.45)
library(ggpubr)
ggarrange(p1, p2, common.legend = T, legend = "bottom")
SavePlot("AgeLength/brms_growth_vb_ppdens", DIR$Figures)

dev.off()
new.graph(0.5)
res1 <- residuals(vb_F_1, type = "pearson") %>%
  as.data.frame() %>%
  bind_cols(age_F) %>%
  mutate(Sex = "Female")
res2 <- residuals(vb_M_1, type = "pearson") %>%
  as.data.frame() %>%
  bind_cols(age_M) %>%
  mutate(Sex = "Male")
p1 <- ggplot(data = bind_rows(res1, res2)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_boxplot(aes(x = factor(age), y = Estimate, fill = Sex), alpha = 0.65) +
  labs(x = "Age", y = "Pearson residual") +
  scale_y_continuous(limits = c(-4, 4)) +
  facet_wrap(Sex ~ .) +
  theme(legend.position = "none") 
print(p1)
SavePlot("AgeLength/brms_growth_vb_resid", DIR$Figures)

sink(make.filename("AgeLength/brms vonB.txt", DIR$Tables))
cat("====================================================================\n")
cat("vonB Females: model = vb_F_1\n")
cat("====================================================================\n")
vb_F_1
cat("\n\n")
cat("====================================================================\n")
cat("vonB Males: model = vb_M_1\n")
cat("====================================================================\n")
vb_M_1
sink()
