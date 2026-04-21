##############################################################################
##############################################################################
# Set up the buffer zone
# basic idea: have a minimum variance threshold so that only compounds 
# with relatively large delta m/z variance are considered
# if within this buffer zone, then likelihood = 1
##############################################################################
##############################################################################
# Load necessary libraries
library(dplyr)
library(ggplot2)

# Define the variance threshold for filtering
var_thresh <- 1.263e-05  # Used mean of all variances as the threshold

# Step 1: Clean the data to ensure proper handling of compound names
cleaned_data <- combined_search_combined %>%
  filter(!is.na(delta_mz) & 
           (!is.na(mol_name) & mol_name != "" | !is.na(name) & name != "")) %>%
  mutate(compound_name = ifelse(!is.na(mol_name) & mol_name != "", as.character(mol_name), as.character(name)))

# Step 2: Calculate variance for each compound
compound_variances <- cleaned_data %>%
  group_by(compound_name) %>%  # Group by compound name
  summarise(variance = var(delta_mz, na.rm = TRUE)) %>%  # Calculate variance of delta m/z with NA removed
  filter(!is.na(variance) & variance > 0)  # Remove any NA variances and ensure variance is positive

# Step 3: Filter compounds that exceed the variance threshold
filtered_variances <- compound_variances %>%
  filter(variance > var_thresh)

# Step 4: Extract the filtered data for further analysis
filtered_compounds <- cleaned_data %>%
  filter(compound_name %in% filtered_variances$compound_name)

# Step 5: Assign likelihood based on variance threshold
filtered_compounds <- filtered_compounds %>%
  mutate(likelihood = ifelse(compound_name %in% filtered_variances$compound_name, NA, 1))  # Assign likelihood = 1 if variance is below threshold

# Step 6: Identify high, medium, and low variance compounds for demonstration
high_variance_compound <- compound_variances %>%
  slice_max(order_by = variance, n = 1)

medium_variance_compound <- compound_variances %>%
  slice_min(order_by = abs(variance - median(variance, na.rm = TRUE)), n = 1)

low_variance_compound <- compound_variances %>%
  slice_min(order_by = variance, n = 1)

# Step 7: Filter the dataset to include only high, medium, and low variance compounds
demo_compounds <- cleaned_data %>%
  filter(compound_name %in% c(high_variance_compound$compound_name,
                              medium_variance_compound$compound_name[1],
                              low_variance_compound$compound_name))

# Step 8: Create density plots to visualize the differences in delta m/z variance
# Extract compound names for high, medium, and low variance compounds
high_compound_name <- as.character(high_variance_compound$compound_name[1])
medium_compound_name <- as.character(medium_variance_compound$compound_name[1])
low_compound_name <- as.character(low_variance_compound$compound_name[1])

# Create the plot with specific colors for high, medium, and low variance compounds
ggplot(demo_compounds, aes(x = delta_mz, fill = compound_name)) +
  geom_density(alpha = 0.6) +
  facet_wrap(~ compound_name, scales = "free_y") +
  scale_fill_manual(values = setNames(
    c("red", "yellow", "green"),
    c(high_compound_name, medium_compound_name, low_compound_name)
  )) +
  labs(title = "Delta m/z Density Distribution for High, Medium, and Low Variance Compounds",
       x = "Delta m/z",
       y = "Density") +
  theme_minimal() +
  theme(legend.position = "none")

##############################################################################
##############################################################################
# Variance Entropy
# Intuition: with smaller sample size, the variance distribution is less reliable
# The idea is to leverage the overall distribution of variances across all compounds 
# to better estimate the variance for each individual compound, 
# especially those with limited sample sizes.

##############################################################################
##############################################################################

# Load necessary libraries
library(dplyr)
library(fitdistrplus)  # For fitting distributions
library(ggplot2)

# Step 1: Calculate Sample Variance for Each Compound
compound_variances <- cleaned_data %>%
  group_by(compound_name) %>%
  summarise(
    variance = var(delta_mz, na.rm = TRUE),
    n = n()
  ) %>%
  filter(!is.na(variance) & variance > 0)  # Remove compounds with NA or zero variance

# Step 2: Fit a Distribution to Model the Overall Variance Distribution
# Fit an inverse-gamma or log-normal distribution to the variances across all compounds
# Here, we'll use log-normal as an example
fit_gamma <- fitdist(compound_variances$variance, "gamma")
# Extract the estimated parameters for the gamma distribution
shape_est <- fit_gamma$estimate["shape"]
rate_est <- fit_gamma$estimate["rate"]

# Step 3: Calculate Shrunk Variance Using Empirical Bayes
# The empirical Bayes shrinkage is applied by combining the sample variance with the global estimate
# Calculate the posterior variance estimate for each compound
compound_variances <- compound_variances %>%
  mutate(
    global_variance = shape_est / rate_est,  # Mean of gamma distribution (shape / rate)
    shrinkage_factor = n / (n + shape_est),  # Shrinkage depends on the sample size and gamma shape parameter
    eb_variance = (shrinkage_factor * variance) + ((1 - shrinkage_factor) * global_variance)  # Empirical Bayes variance
  )

# Step 4: Compare Sample Variance and Empirical Bayes Variance
ggplot(compound_variances, aes(x = variance, y = eb_variance)) +
  geom_point(alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Sample Variance vs. Empirical Bayes Variance (Gamma Prior)",
       x = "Sample Variance",
       y = "Empirical Bayes Variance") +
  theme_minimal()

# Step 5: Visualize the Shrunk Variance Distribution
ggplot(compound_variances, aes(x = variance)) +
  geom_density(alpha = 0.6, fill = "lightgreen") +
  labs(title = "Density of Sample Variance",
       x = "Variance",
       y = "Density") +
  theme_minimal()



##############################################################################
##############################################################################
# Explore RT distribution
# Plot RT density distribution for each compound (counts in top 5 and last 5 in this case)
##############################################################################
##############################################################################
# Load necessary libraries
library(dplyr)
library(ggplot2)

# Step 1: Filter the dataset to keep only entries from annotation library with valid RT values, excluding names starting with "zz" or "yy"
filtered_annotation_data <- annotation_search_combined %>%
  mutate(compound_name = ifelse(!is.na(mol_name) & mol_name != "", as.character(mol_name), as.character(name))) %>%
  filter(!grepl("^(zz|yy)", compound_name, ignore.case = TRUE))  # Exclude names starting with "zz" or "yy"

# Step 2: Calculate counts for each compound within the filtered annotation dataset
compound_counts <- filtered_annotation_data %>%
  group_by(compound_name) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

# Step 3: Select the top 10 compounds based on counts
top_10_compounds <- compound_counts %>%
  slice_max(order_by = count, n = 9, with_ties = FALSE)

# Step 4: Filter the dataset to keep only the top 10 selected compounds
filtered_top_10_data <- filtered_annotation_data %>%
  filter(compound_name %in% top_10_compounds$compound_name)

# Step 4: Plot RT density distribution for each selected compound using ggplot2
ggplot(filtered_top_10_data, aes(x = rt, fill = compound_name)) +
  geom_density(alpha = 0.6) +
  facet_wrap(~ compound_name, scales = "free_y") +
  labs(title = "RT Density Distribution for Top 10 Compounds",
       x = "Retention Time (RT)",
       y = "Density") +
  theme_minimal() +
  theme(legend.position = "none")


# Step 4.1: Plot RT histogram for each selected compound using ggplot2
ggplot(filtered_top_10_data, aes(x = rt, fill = compound_name)) +
  geom_histogram(bins = 10, alpha = 0.6, position = "identity") +
  facet_wrap(~ compound_name, scales = "free_y") +
  labs(title = "RT Histogram for Top Compounds (Excluding 'zz_no peak')",
       x = "Retention Time (RT)",
       y = "Count") +
  theme_minimal() +
  theme(legend.position = "none")

# Step 4.2: Plot RT boxplot for each selected compound using ggplot2
ggplot(filtered_top_10_data, aes(x = compound_name, y = rt, fill = compound_name)) +
  geom_boxplot(alpha = 0.6) +
  labs(title = "RT Boxplot for Top Compounds",
       x = "Compound Name",
       y = "Retention Time (RT)") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability

## doesnt seem to be super doable because library entries does not RT entries. 
## just proceed with annotation data