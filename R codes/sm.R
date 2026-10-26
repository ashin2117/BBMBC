### Supplementary Materials
# You must restart your R session to run this code if you have already install BBMBC package.
install.packages(c("DHARMa", "glmmTMB"))

library(DHARMa)
library(glmmTMB)
library(dplyr)

df_app <- read.csv("cd4b_igg_wk24.csv")

analysis_df <- df_app %>%
  # Remove any rows with missing data (e.g., the blank row at the end of the file)
  filter(!is.na(group)) %>%
  # Rename columns to match the existing model logic
  rename(
    parent_count = "iggb",
    target_count = "cd4b_iggb"
  ) %>%
  mutate(
    group = as.factor(group), 
    failures = parent_count - target_count,
    prop = target_count / parent_count
  )

resp_matrix <- as.matrix(analysis_df[, c("target_count", "failures")])


fit_bb <- glmmTMB(resp_matrix ~ group,
                  dispformula = ~ group,
                  family = glmmTMB::betabinomial(link = "logit"),
                  data = analysis_df)

res_bb <- simulateResiduals(fittedModel = fit_bb)

# plot the standard diagnostic dashboard for Beta Binomial model
plot(res_bb)

fit_b <- glmmTMB(resp_matrix ~ group,
                 family = binomial(link = "logit"),
                 data = analysis_df)

res_b <- simulateResiduals(fittedModel = fit_b)

# plot the standard diagnostic dashboard for Binomial model
plot(res_b)