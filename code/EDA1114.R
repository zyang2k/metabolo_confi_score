# 1114 EDA


# Group by mol_name and calculate summary statistics for delta_mz
delta_mz_summary <- combined_search_combined %>%
  group_by(mol_name) %>%
  summarise(
    mean_delta_mz = mean(delta_mz, na.rm = TRUE),
    median_delta_mz = median(delta_mz, na.rm = TRUE),
    sd_delta_mz = sd(delta_mz, na.rm = TRUE),
    min_delta_mz = min(delta_mz, na.rm = TRUE),
    max_delta_mz = max(delta_mz, na.rm = TRUE),
    count = n()
  )

head(delta_mz_summary)

# select weird values of delta m/z
# Set threshold for what you consider as "weird"
sd_threshold <- 0.005  # Example: standard deviation greater than 0.005
mean_threshold <- 0.003  # Example: mean_delta_mz greater than 0.003 or less than -0.003

# Filter the rows based on the criteria
weird_delta_mz <- delta_mz_summary %>%
  filter(
    abs(mean_delta_mz) > mean_threshold |  # Absolute value of mean greater than the threshold
      sd_delta_mz > sd_threshold |           # Standard deviation greater than the threshold
      is.na(sd_delta_mz)                     # Include rows where sd is NA (e.g., only 1 observation)
  )

# View the filtered results
head(weird_delta_mz)

##############################################################################
############## Retention time Bayesian
##############################################################################

citric_acid_rt_data <- combined_search_combined %>%
  filter(mol_name == "citric acid" | grepl("citric acid", name, ignore.case = TRUE)) %>%
  pull(rt)

# Remove NAs
citric_acid_rt_data <- na.omit(citric_acid_rt_data)

summary(citric_acid_rt_data)

# Fit KDE
rt_kde_fit <- density(citric_acid_rt_data)

# Plot the KDE
ggplot(data.frame(rt = citric_acid_rt_data), aes(x = rt)) +
  geom_density(fill = "blue", alpha = 0.5) +
  labs(title = "Density Plot of Retention Time for Citric Acid",
       x = "Retention Time",
       y = "Density") +
  theme_minimal()


new_rt <- 137.0  # Example value

# Find the likelihood using KDE
rt_likelihood <- approx(rt_kde_fit$x, rt_kde_fit$y, xout = new_rt, rule = 2)$y

total_counts <- sum(mol_name_counts)

# Compute the evidence
rt_evidence <- sum(sapply(names(mol_name_counts), function(name) {
  # Filter data for the compound
  compound_data <- combined_search_combined %>%
    filter(mol_name == name | grepl(name, mol_name, ignore.case = TRUE)) %>%
    pull(rt)
  
  # Proceed only if we have at least two points
  if (length(compound_data) >= 2) {
    # Fit KDE for each compound
    kde <- density(compound_data, na.rm = TRUE)
    rt_likelihood_compound <- approx(kde$x, kde$y, xout = new_rt)$y
    prior_compound <- mol_name_counts[name] / total_counts
    return(rt_likelihood_compound * prior_compound)
  } else {
    # If not enough data points, skip this compound
    return(0)
  }
}), na.rm = TRUE)

# Use mol_name_counts as prior
prior_citric_acid <- mol_name_counts["citric acid"] / total_counts

# Calculate the posterior probability
if (!is.na(rt_likelihood) && !is.na(rt_evidence) && rt_evidence != 0) {
  rt_posterior <- (rt_likelihood * prior_citric_acid) / rt_evidence
  cat("Posterior probability for citric acid given new RT:", rt_posterior, "\n")
} else {
  cat("Unable to calculate posterior probability due to missing or zero values.\n")
}


