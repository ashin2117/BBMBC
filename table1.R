# Install required packages
install.packages(c("dplyr", "purrr", "tidyr", "foreach", "doParallel", "BBMBC", "here"))

# Load the libraries
library(dplyr)
library(purrr)
library(tidyr)
library(foreach)
library(doParallel)
library(BBMBC)

source("calc_t1_wilcox.R")

### Data Manipulation
df <- read.csv("BG505_MD_bcp.csv")
df$MAGNITUDE <- df$MAGNITUDE/100

# group and summarize directly
summary_results <- df %>%
  group_by(RX_LABEL, MONTH) %>%
  summarise(
    mean_magnitude = mean(MAGNITUDE, na.rm = TRUE),
    rho_mom = calc_rho_mom(MAGNITUDE),
    sample_sizes = n(),
    .groups = "drop"
  ) %>% 
  mutate(dataset = paste(RX_LABEL, MONTH, sep = "_")) %>%
  select(dataset, mean_magnitude, rho_mom, sample_sizes)

# write.csv(summary_results, "dataset.csv", row.names = FALSE)



### Simulation Data Generation
# Sample size scenarios bounded by observed [12, 18]
n_pairs <- bind_rows(
  data.frame(n0 = 15, n1 = 15, n_design = "Balanced"),
  data.frame(n0 = 12, n1 = 18, n_design = "Imbalanced")
)

# Rho scenarios bounded by observed [~0.002, ~0.032]
rho_pairs <- bind_rows(
  data.frame(rho0 = c(0.005, 0.01, 0.03), rho1 = c(0.005, 0.01, 0.03), 
             rho_design = "Equal Dispersion"),
  data.frame(rho0 = 0.005, rho1 = 0.03, 
             rho_design = "Variance Inflation"),
  data.frame(rho0 = 0.03, rho1 = 0.005, 
             rho_design = "Variance Reduction")
)

# Mean magnitudes bounded by observed [~0.006, ~0.047]
prop_values <- data.frame(prop = c(0.01, 0.03, 0.05))

empirical_sim_grid <- rho_pairs %>%
  cross_join(n_pairs) %>%
  cross_join(prop_values) %>%
  mutate(sim_id = row_number()) %>%
  select(sim_id, rho0, rho1, rho_design, n0, n1, n_design, prop)



### Simulation (Type I error rate)
# Set global parameters for the run
num_sims <- 5000 
cores_per_sim <- parallel::detectCores() - 1

# Iterate over the grid
full_simulation_results <- empirical_sim_grid %>%
  mutate(
    # pmap iterates over lists of arguments row by row
    metrics = pmap(
      list(n0, n1, prop, rho0, rho1, sim_id),
      function(n0_val, n1_val, prop_val, rho0_val, rho1_val, id) {
        
        message(sprintf("Evaluating Scenario %d / %d...", id, nrow(empirical_sim_grid)))
        
        # Compute Beta-Binomial Type-I Error
        bb_res <- calc_t1_bb(
          n_0 = n0_val, n_1 = n1_val, p_base = prop_val,
          rho_0 = rho0_val, rho_1 = rho1_val,
          n_sims = num_sims, n_cores = cores_per_sim
        )
        
        # Compute Wilcoxon Type-I Error
        wilcox_res <- calc_t1_wilcox(
          n_0 = n0_val, n_1 = n1_val, p_base = prop_val,
          rho_0 = rho0_val, rho_1 = rho1_val,
          n_sims = num_sims, n_cores = cores_per_sim
        )
        
        # Return the extracted metrics as a temporary data frame
        data.frame(
          t1_error_bb = bb_res$type_1_error_rate,
          t1_error_wilcox = wilcox_res$type_1_error_rate,
          conv_bb = bb_res$convergences,
          conv_wilcox = wilcox_res$convergences
        )
      }
    )
  ) %>%
  # Unnest expands the nested dataframes into standard columns
  unnest(metrics)

# write.csv(full_simulation_results, "t1_manuscript.csv")