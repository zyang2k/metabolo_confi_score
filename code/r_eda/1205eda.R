# Load necessary libraries
library(ggplot2)
library(dplyr)
library(purrr)
library(car)


################# Library data EDA ###############################################

# Basic Statistics for Key Features
summary_stats <- summary(library_search_combined[, c("delta_mz", "precursor_mz", "charge")])
print("Basic Statistics for Key Features:")
print(summary_stats)

# Identify Common Adducts
common_adducts <- table(library_search_combined$adduct)
print("Common Adducts:")
print(sort(common_adducts, decreasing = TRUE))

# Identify Common Databases
common_databases <- table(library_search_combined$db)
print("Common Databases:")
print(sort(common_databases, decreasing = TRUE))

################# sample 5000 for further analysis ###############################################

# Set the sample size (e.g., 5000 points)
sample_size <- 5000

# Randomly sample indices from the dataset
set.seed(123) # Set seed for reproducibility
sample_indices <- sample(nrow(library_search_combined), sample_size)

# Subset the dataset with the sampled indices
sampled_data <- library_search_combined[sample_indices, ]


################# delta-mz histogram ###############################################
# to look a total distribution of delta m/z
# identifiy things that worth looking at

# Create a histogram in base R
hist_data <- hist(
  sampled_data$delta_mz, 
  breaks = 50, 
)

# Extract bin edges and frequencies
bin_edges <- hist_data$breaks
frequencies <- hist_data$counts

# Combine bin edges and frequencies into a data frame
bin_summary <- data.frame(
  bin_start = bin_edges[-length(bin_edges)], # Start of each bin
  bin_end = bin_edges[-1],                  # End of each bin
  frequency = frequencies
)

# Sort the data frame by frequency in descending order
bin_summary_sorted <- bin_summary[order(-bin_summary$frequency), ]

# Display the top bins
head(bin_summary_sorted)



# weird bin 1: [0.0040, 0.0045]

# Subset data for the bin [0.0040, 0.0045]
bin_004_0045_data <- library_search_combined[
  library_search_combined$delta_mz >= 0.0040 & library_search_combined$delta_mz < 0.0045, 
]

# Check the structure and summary of the subset
str(bin_004_0045_data)
summary(bin_004_0045_data)


# Frequency tables for metadata
table(bin_004_0045_data$adduct)  # Adduct types
table(bin_004_0045_data$db)      # Databases contributing to this bin


plot(
  bin_004_0045_data$precursor_mz, bin_004_0045_data$delta_mz,
  main = "Precursor m/z vs Delta m/z in Bin [0.0040, 0.0045]",
  xlab = "Precursor m/z",
  ylab = "Delta m/z",
  col = "blue"
)
grid()


sort(table(bin_004_0045_data$name), decreasing = TRUE)[1:10]


################# delta-mz by adduct ###############################################

# group by adducts

# Identify the top 10 adducts by count
top_adducts <- sampled_data %>%
  count(adduct) %>%
  arrange(desc(n)) %>%
  slice_head(n = 10) %>%
  pull(adduct)

# Filter the dataset for only the top 6 adducts
filtered_data <- sampled_data %>%
  filter(adduct %in% top_adducts)

# Create a boxplot for delta m/z of the top 6 adducts
ggplot(filtered_data, aes(x = adduct, y = delta_mz, fill = adduct)) +
  geom_boxplot(alpha = 0.7) +
  labs(
    title = "Delta m/z Distribution for Top 10 Adducts",
    x = "Adduct",
    y = "Delta m/z"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# [M-H-CH3]-" [M-2H]2-" looks weird

# stat test to see if they are v-different
# Load necessary library

# Ensure adduct is a factor
sampled_data <- sampled_data %>%
  mutate(adduct = as.factor(adduct))

# Perform one-way ANOVA
anova_result <- aov(delta_mz ~ adduct, data = sampled_data)

# Display the ANOVA summary
summary(anova_result)
# significant difference between adudct group


# Check for assumptions: Homogeneity of variance (Levene's Test)
leveneTest(delta_mz ~ adduct, data = sampled_data)

# violated, but why?
# does adduct affect instrumental accuracy?


# kruskal test
kruskal.test(delta_mz ~ adduct, data = sampled_data)


# shrinked pairwise wilcox to identify the difference. 
filtered_data <- sampled_data %>%
  filter(adduct %in% top_adducts)

pairwise.wilcox.test(filtered_data$delta_mz, filtered_data$adduct, p.adjust.method = "bonferroni")


################# accurate mass vs. delta mz? ###############################################

# Scatter plot using the sampled data
plot(
  sampled_data$precursor_mz,
  sampled_data$delta_mz,
  main = "Relationship Between Precursor_mz and Delta_mz (Sampled Data)",
  xlab = "Precursor m/z",
  ylab = "Delta m/z",
  pch = 19,            # Solid circle points
  col = rgb(0, 0, 1, 0.6), # Semi-transparent blue
  cex = 0.7            # Adjust point size
)



################# # QUESTIONS ###############################################

# 1. a lot of entries does not have a msms match that passes the threshold, 
# what to do then?
#   
# 2. Rightnow the workflow is to start from msms match directly. 
# Which is kinda different from the proposed 5-1 workflow.
# Shall i start de novo?
#   
# 3. in the 5-1 workflow, what feature is considered to be orthogonal?
# 


