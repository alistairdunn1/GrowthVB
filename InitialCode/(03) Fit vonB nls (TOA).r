library(ggplot2)
library(viridis)
library(tidyverse)

theme_set(t1())

source(make.filename("Table functions.r", DIR$Functions))
source(make.filename("08-LengthAge/Functions/Functions.r", DIR$R))

# use data from "(01) Length age (TOA).r"
load(file = make.filename("Age data/age_data.rdata", DIR$Data), verbose = TRUE)
data <- age.data
bad <- c(which(data$nage < 1))
#data <- data[-bad, ]
data <- data %>% filter(CCAMLR_year %in% 2023:2024)

a <- fit_vb(sex = data$sex, length = data$length, age = data$nage, length_bins = NULL, sampling_prob = 1)
a$parameters
csv(a$parameters, "AgeLength/NLS vonB", DIR$Tables)

# Combine original data with model fits
data <- bind_cols(data, a$data)

# Plot fit
dev.off()
new.graph(0.5)
p1 <- ggplot(a$data, aes(x = age, y = length, colour = sex)) +
  geom_point(alpha = 0.5) +
  geom_line(data = a$fits, aes(x = age, y = mean), colour = "black") +
  geom_line(data = a$fits, aes(x = age, y = lowerCI), colour = "black", linetype = "dashed") +
  geom_line(data = a$fits, aes(x = age, y = upperCI), colour = "black", linetype = "dashed") +
  xlim(0, NA) +
  ylim(0, NA) +
  ylab("Length (cm)") +
  xlab("Age") +
  facet_wrap(~sex) +
  st1
print(p1)
SavePlot("Length age/VonB fits")

# Add diagnostic plots

plot(data$student)
plot(data$residual)

ggplot(a$data, aes(x = sex, y = residual, fill = sex)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Sex",y = "Residuals")

