# Installation
install.packages(c("dplyr", "tidyr", "purrr", "ggplot2", "BBMBC"))

# Loading (run at the top of your script)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(BBMBC)

source("calc_power_wilcox.R")

fold_changes <- seq(1.0, 2, by = 0.1)

rho_0_emp  <- 0.0007
rho_1_emp  <- 0.001
p_emp      <- 0.00152

# define the combinations of sample sizes
sample_sizes <- list(
  c(n0 = 13, n1 = 17),
  c(n0 = 17, n1 = 17),
  c(n0 = 17, n1 = 20),
  c(n0 = 20, n1 = 20)
)

### Simulation (Power)
all_power_results <- lapply(sample_sizes, function(ss) {
  n0 <- ss["n0"]
  n1 <- ss["n1"]
  
  # create a label for faceting
  sample_label <- paste0("n_0 = ", n0, ", n_1 = ", n1)
  
  cat("\nRunning simulations for:", sample_label, "\n")
  
  # iterate through fold changes for this specific sample size pair
  power_results <- lapply(fold_changes, function(fc) {
    cat("  Fold change:", fc, "\n")
    
    # beta-binomial simulation
    res_bb <- calc_power_bb(
      n_0 = n0,
      n_1 = n1,
      p_base = p_emp,
      fold_change = fc,
      rho_0 = rho_0_emp,
      rho_1 = rho_1_emp,
      seed = 2026,
      n_sims = 500  
    )
    
    # Wilcoxon Rank Sum Simulation
    res_wilcox <- calc_power_wilcox(
      n_0 = n0,
      n_1 = n1,
      p_base = p_emp,
      fold_change = fc,
      rho_0 = rho_0_emp,
      rho_1 = rho_1_emp,
      seed = 2026,
      n_sims = 500  
    )
    
    rbind(
      data.frame(Sample_Size = sample_label, fold_change = fc, Power = res_bb$power, Method = "Beta-Binomial"),
      data.frame(Sample_Size = sample_label, fold_change = fc, Power = res_wilcox$power, Method = "Wilcoxon Rank-Sum")
    )
  })
  
  # Combine the fold change list into one dataframe for this sample size pair
  do.call(rbind, power_results)
})

# Combine all sample size results into one master data frame
plot_data <- do.call(rbind, all_power_results)

# Ensure Sample_Size is a factor so it orders correctly in the facet grid
plot_data$Sample_Size <- factor(plot_data$Sample_Size, 
                                levels = c("n_0 = 13, n_1 = 17", 
                                           "n_0 = 17, n_1 = 17", 
                                           "n_0 = 17, n_1 = 20", 
                                           "n_0 = 20, n_1 = 20"))

### Plotting
# Generate the faceted comparative power plot
figure3 <- ggplot(plot_data, aes(x = fold_change, y = Power, color = Method, linetype = Method, shape = Method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  # Create a 2x2 grid based on the sample sizes
  facet_wrap(~ Sample_Size, nrow = 2) +
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

print(figure3)