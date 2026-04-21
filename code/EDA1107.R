# wiki_data structure:
# a large list
# each sublist corresponds to a wiki_id
# each sublist has two DataFrames: 
# annotation_search: contains annotation search results.
# library_search: contains library search results.
# each annotation df has 12 columns, each library df has 10 columns


##############################################################################
# Extract annotation_search and library_search data from wiki_data
##############################################################################

annotation_search_list <- list()
library_search_list <- list()

# Iterate over each wiki_id in wiki_data
for (wiki_id in names(wiki_data)) {
  annotation_search_df <- wiki_data[[wiki_id]]$annotation_search
  library_search_df <- wiki_data[[wiki_id]]$library_search
  
  # Only add non-empty DataFrames
  if (nrow(annotation_search_df) > 0) {
    annotation_search_list[[wiki_id]] <- annotation_search_df
  }
  
  if (nrow(library_search_df) > 0) {
    library_search_list[[wiki_id]] <- library_search_df
  }
}



# List to store annotation search DataFrames
annotation_search_list <- lapply(wiki_data, function(x) x$annotation_search)

# Ensure all DataFrames have the same columns before combining
annotation_search_list <- lapply(annotation_search_list, function(df) {
  if (ncol(df) == 0) {
    df <- setNames(data.frame(matrix(ncol = length(target_columns), nrow = 0)), target_columns)
  }
  df
})

# Combine the DataFrames into one
annotation_search_combined <- bind_rows(annotation_search_list)
library_search_combined <- do.call(rbind, library_search_list)



##############################################################################
# summary stats
##############################################################################
# Summary statistics for numerical columns in the combined DataFrames
summary(annotation_search_combined$precursor_mz)
summary(library_search_combined$precursor_mz)

summary(annotation_search_combined$delta_mz)
summary(library_search_combined$delta_mz)

# Calculate molecular name counts
combined_search_combined <- bind_rows(annotation_search_combined, library_search_combined)

mol_name_counts <- table(combined_search_combined$mol_name)
mol_name_counts <- sort(mol_name_counts, decreasing = TRUE)




##############################################################################
# Visualize delta_mz distribution
##############################################################################

annotation_search_combined$source <- "annotation"
library_search_combined$source <- "library"

combined_delta_mz <- rbind(annotation_search_combined[, c("delta_mz", "source")], library_search_combined[, c("delta_mz", "source")])

# Plot the distribution of delta_mz values
ggplot(combined_delta_mz, aes(x = delta_mz, fill = source)) +
  geom_histogram(binwidth = 0.001, alpha = 0.5, position = "identity") +
  labs(title = "Comparison of Delta m/z Values by Source", x = "Delta m/z", y = "Frequency") +
  scale_fill_manual(values = c("red", "green"))



# Plot density for delta_mz in annotation and library data
ggplot() +
  geom_density(data = annotation_search_combined, aes(x = delta_mz), fill = "blue", alpha = 0.5) +
  geom_density(data = library_search_combined, aes(x = delta_mz), fill = "red", alpha = 0.5) +
  labs(title = "Density Plot of Delta m/z Values", x = "Delta m/z", y = "Density") +
  theme_minimal()


##############################################################################
# Explore zz entries
##############################################################################

zz_annotation_list <- lapply(wiki_data, function(x) x$annotation_search)
zz_combined_df <- bind_rows(zz_annotation_list) %>% filter(grepl("zz", tolower(mol_name)))

# Plot histogram of delta m/z for zz entries
ggplot(zz_combined_df, aes(x = delta_mz)) +
  geom_histogram(binwidth = 0.001, alpha = 0.5, position = "identity", fill = "blue") +
  labs(title = "Frequency Plot of Delta m/z Values of zz entries", x = "Delta m/z", y = "Frequency") +
  theme_minimal()

##############################################################################
# Explore Bayesian
##############################################################################

# Step 1: Extract delta m/z values for 'citric acid'
citric_acid_data <- combined_search_combined %>%
  filter(mol_name == "citric acid" | grepl("citric acid", name, ignore.case = TRUE)) %>%
  pull(delta_mz)


# Remove NAs from the data
citric_acid_data <- na.omit(citric_acid_data)

# Step 2: Fit Gaussian distribution
mean_citric_acid <- mean(citric_acid_data)
sd_citric_acid <- sd(citric_acid_data)

# Step 3: Create a density plot and overlay the Gaussian fit
## The Gaussian fit is not okay!
ggplot(data.frame(delta_mz = citric_acid_data), aes(x = delta_mz)) +
  geom_density(aes(y = ..density..), fill = "blue", alpha = 0.5) +
  stat_function(fun = dnorm, args = list(mean = mean_citric_acid, sd = sd_citric_acid), color = "red", size = 1) +
  labs(title = "Density Plot of Delta m/z for Citric Acid with Gaussian Fit",
       x = "Delta m/z",
       y = "Density") +
  theme_minimal()

# Step 3.1: Use non-parametric
ggplot(data.frame(delta_mz = citric_acid_data), aes(x = delta_mz)) +
  geom_density(fill = "blue", alpha = 0.5) +
  labs(title = paste("Density Plot of Delta m/z for", specific_compound),
       x = "Delta m/z", y = "Density") +
  theme_minimal()

# Fit KDE
kde_fit <- density(citric_acid_data)

# Step 3: Calculate the likelihood for a new delta m/z value
new_delta_mz <- 0.000001

# Interpolate to find the likelihood for the new delta m/z value
likelihood <- approx(kde_fit$x, kde_fit$y, xout = new_delta_mz, rule = 2)$y

# Print the likelihood value
if (is.na(likelihood)) {
  cat("The new delta m/z value is outside the range of the KDE fit.\n")
} else {
  cat("Likelihood of new delta m/z (annotation):", likelihood, "\n")
}


# Step 4: Calculate the evidence (P(D))
# The evidence is the sum of all possible outcomes weighted by their priors.
# For simplicity, we assume the density function is normalized, so P(D) = sum(likelihood * prior)
# Here, we'll approximate the evidence as the likelihood from KDE weighted by the prior.
# Note: This is a simplified way to approximate the evidence.

evidence <- sum(sapply(names(mol_name_counts), function(name) {
  compound_data <- combined_search_combined %>%
    filter(mol_name == name) %>%
    pull(delta_mz)
  kde <- density(compound_data)
  likelihood_compound <- approx(kde$x, kde$y, xout = new_delta_mz)$y
  prior_compound <- mol_name_counts[name] / total_counts
  return(likelihood_compound * prior_compound)
}), na.rm = TRUE)

# Step 5: Calculate the posterior probability
posterior <- (likelihood * prior_citric_acid) / evidence

##############################################################################
# Distribution by name
##############################################################################
# problem: so many NAs! a lot of entries that has no name. 

