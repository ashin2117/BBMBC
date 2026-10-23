# Installation
install.packages(c("dplyr", "tidyr", "purrr", "ggplot2", "BBMBC"))

# Loading (run at the top of your script)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(BBMBC)

source("calc_power_wilcox.R")

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



### Simulation (Power)
num_sims <- 500 
cores_per_sim <- parallel::detectCores() - 1

target_fold_changes <- c(1, 1.25, 1.5, 1.75, 2, 2.5, 3)

# Expand the simulation grid
expanded_sim_grid <- empirical_sim_grid %>%
  tidyr::crossing(fold_change = target_fold_changes)

# Iterate over the expanded grid
full_simulation_results <- expanded_sim_grid %>%
  mutate(
    # Add fold_change to the pmap argument list
    metrics = pmap(
      list(n0, n1, prop, rho0, rho1, fold_change, sim_id),
      function(n0_val, n1_val, prop_val, rho0_val, rho1_val, fc_val, id) {
        
        # Optional: Print progress with the specific fold change
        message(sprintf("Evaluating Scenario %d with Fold Change %.1f...", id, fc_val))
        
        # Compute Beta-Binomial Power
        bb_res <- calc_power_bb(
          n_0 = n0_val, 
          n_1 = n1_val, 
          p_base = prop_val,
          fold_change = fc_val, # Passed dynamically
          rho_0 = rho0_val, 
          rho_1 = rho1_val,
          n_sims = num_sims, 
          n_cores = cores_per_sim
        )
        
        # Compute Wilcoxon Power
        wilcox_res <- calc_power_wilcox(
          n_0 = n0_val, 
          n_1 = n1_val, 
          p_base = prop_val,
          fold_change = fc_val, # Passed dynamically
          rho_0 = rho0_val, 
          rho_1 = rho1_val,
          n_sims = num_sims, 
          n_cores = cores_per_sim
        )
        
        # Return the extracted metrics as a temporary data frame
        data.frame(
          power_bb = bb_res$power,
          power_wilcox = wilcox_res$power,
          conv_bb = bb_res$convergences,
          conv_wilcox = wilcox_res$convergences
        )
      }
    )
  ) %>%
  # Unnest expands the nested dataframes into standard columns
  unnest(metrics)

# write.csv(full_simulation_results, "power_manuscript_multiple_folds.csv", row.names = FALSE)


### Plotting
power_df <- full_simulation_results

## Figure 1
figure1_data <- power_df %>%
  filter(
    prop == 0.01,
    (rho0 == 0.03 & rho1 == 0.03) |
      (rho0 == 0.005 & rho1 == 0.03) |
      (rho0 == 0.03 & rho1 == 0.005)
  ) %>%
  pivot_longer(
    cols = c(power_bb, power_wilcox),
    names_to = "Method",
    values_to = "Power"
  ) %>%
  mutate(
    Method = if_else(Method == "power_bb", "Beta-Binomial", "Wilcoxon Rank-Sum"),
    Dispersion = case_when(
      rho0 == 0.03 & rho1 == 0.03 ~ "(0.03, 0.03)",
      rho0 == 0.005 & rho1 == 0.03 ~ "(0.005, 0.03)",
      rho0 == 0.03 & rho1 == 0.005 ~ "(0.03, 0.005)"
    ),
    Dispersion = factor(Dispersion, levels = c(
      "(0.03, 0.03)", 
      "(0.005, 0.03)", 
      "(0.03, 0.005)"
    )),
    Sample_Size = sprintf("n0 = %d, n1 = %d", n0, n1),
    Sample_Size = factor(Sample_Size, levels = c(
      "n0 = 15, n1 = 15", 
      "n0 = 12, n1 = 18"
    ))
  ) %>%
  # Keep all Beta-Binomial data, but only keep Wilcoxon if dispersion is Equal
  filter(Method == "Beta-Binomial" | Dispersion == "(0.03, 0.03)")

# Build the ggplot object
figure1 <- ggplot(figure1_data, aes(x = fold_change, y = Power, color = Method, linetype = Method, shape = Method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  facet_grid(Dispersion ~ Sample_Size) +
  scale_color_manual(values = c("Beta-Binomial" = "#0072B2", "Wilcoxon Rank-Sum" = "#D55E00")) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
  labs(
    x = "Fold Change",
    y = "Statistical Power") +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", color = "black"),
    strip.text = element_text(face = "bold")
  )

# ggsave("power_curves_samples(p=0.01).pdf", plot = p, device = grDevices::cairo_pdf, width = 10, height = 7)

## Figure 2
figure2_data <- power_df %>%
  filter(
    prop == 0.03,
    (rho0 == 0.03 & rho1 == 0.03) |
      (rho0 == 0.005 & rho1 == 0.03) |
      (rho0 == 0.03 & rho1 == 0.005)
  ) %>%
  pivot_longer(
    cols = c(power_bb, power_wilcox),
    names_to = "Method",
    values_to = "Power"
  ) %>%
  mutate(
    Method = if_else(Method == "power_bb", "Beta-Binomial", "Wilcoxon Rank-Sum"),
    Dispersion = case_when(
      rho0 == 0.03 & rho1 == 0.03 ~ "(0.03, 0.03)",
      rho0 == 0.005 & rho1 == 0.03 ~ "(0.005, 0.03)",
      rho0 == 0.03 & rho1 == 0.005 ~ "(0.03, 0.005)"
    ),
    Dispersion = factor(Dispersion, levels = c(
      "(0.03, 0.03)", 
      "(0.005, 0.03)", 
      "(0.03, 0.005)"
    )),
    Sample_Size = sprintf("n0 = %d, n1 = %d", n0, n1),
    Sample_Size = factor(Sample_Size, levels = c(
      "n0 = 15, n1 = 15", 
      "n0 = 12, n1 = 18"
    ))
  ) %>%
  # Keep all Beta-Binomial data, but only keep Wilcoxon if dispersion is Equal
  filter(Method == "Beta-Binomial" | Dispersion == "(0.03, 0.03)")

# Build the ggplot object
figure2 <- ggplot(figure2_data, aes(x = fold_change, y = Power, color = Method, linetype = Method, shape = Method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  facet_grid(Dispersion ~ Sample_Size) +
  scale_color_manual(values = c("Beta-Binomial" = "#0072B2", "Wilcoxon Rank-Sum" = "#D55E00")) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
  labs(
    x = "Fold Change",
    y = "Statistical Power") +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", color = "black"),
    strip.text = element_text(face = "bold")
  )

# ggsave("power_curves_samples(p=0.03).pdf", plot = p, device = grDevices::cairo_pdf, width = 10, height = 7)




