# analysis of parent-child relationship
# functions stored in get_parent_child.R

# data comes from compound, 5m hilic premier | orbitrap | beh amide | negative
# > dim(data)
# [1] 235940     52


# tree_data <- build_family_trees(data)
tree_data <- tree_data %>%
  left_join(sample_metadata %>% select(sample_id, content),
            by = c("sample" = "sample_id"))

# > dim(tree_data)
# [1] 235940     59

direct_similarities <- calculate_direct_similarities(tree_data)

# visualization

# Histogram of entropy similarity
ggplot(direct_similarities, aes(x = entropy_similarity)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black") +
  labs(title = "Distribution of Entropy Similarity", x = "Entropy Similarity", y = "Count")

# Boxplot comparing entropy similarity by same_sample
ggplot(direct_similarities, aes(x = factor(same_sample), y = entropy_similarity)) +
  geom_boxplot(fill = c("lightgreen", "lightcoral")) +
  labs(title = "Entropy Similarity by Same Sample", x = "Same Sample (FALSE = No, TRUE = Yes)", y = "Entropy Similarity")



# Remove NA values and create the boxplot
boxplot(direct_similarities$entropy_similarity ~ direct_similarities$same_species, 
        data = direct_similarities, 
        main = "Entropy Similarity by Species Relationship",
        xlab = "Same Species", 
        ylab = "Entropy Similarity")

counts <- table(direct_similarities$same_species)
text(1:2, par("usr")[3], 
     labels = paste("n =", counts), 
     pos = 1, 
     xpd = TRUE)


boxplot(direct_similarities$entropy_similarity ~ direct_similarities$same_organ, 
        data = direct_similarities, 
        main = "Entropy Similarity by Organ Relationship",
        xlab = "Same Organ", 
        ylab = "Entropy Similarity")

counts <- table(direct_similarities$same_organ)
text(1:2, par("usr")[3], 
     labels = paste("n =", counts), 
     pos = 1, 
     xpd = TRUE)


# Filter for mouse species parent-children pairs
mouse_pairs <- direct_similarities[
  direct_similarities$same_species == TRUE & 
    direct_similarities$parent_species == "mouse",
  "entropy_similarity"
]

# Count the number of pairs
num_mouse_pairs <- length(mouse_pairs$entropy_similarity)

# Create histogram
hist(mouse_pairs$entropy_similarity, 
     main = paste("Entropy Similarity Distribution for Mouse Species\n(n =", num_mouse_pairs, "pairs)"), 
     xlab = "Entropy Similarity", 
     ylab = "Frequency",
     col = "skyblue",
     border = "black")


# Filter for human species, same species, and non-NA entropy similarities
human_similarities <- direct_similarities[
  direct_similarities$same_species == TRUE & 
    !is.na(direct_similarities$entropy_similarity) & 
    direct_similarities$parent_species == "human",
  "entropy_similarity"
]

# Create histogram
hist(human_similarities$entropy_similarity, 
     main = "Entropy Similarity Distribution\nfor Human Species Metabolites", 
     xlab = "Entropy Similarity", 
     ylab = "Frequency",
     col = "lightgreen",
     border = "black")



# Count and display number of human pairs
num_human_pairs <- length(human_similarities$entropy_similarity)
text(x = par("usr")[2] * 0.7, y = par("usr")[4] * 0.9, 
     labels = paste("n =", num_human_pairs, "pairs"), 
     adj = c(0, 0))



# Filter for Plant species, same species, and non-NA entropy similarities
plant_similarities <- direct_similarities[
  direct_similarities$same_species == TRUE & 
    !is.na(direct_similarities$entropy_similarity) & 
    direct_similarities$parent_species == "Plant",
  "entropy_similarity"
]

# Count the number of plant pairs
num_plant_pairs <- length(plant_similarities$entropy_similarity)

# Create histogram
hist(plant_similarities$entropy_similarity, 
     main = paste("Entropy Similarity Distribution for Plant Species\n(n =", num_plant_pairs, "pairs)"), 
     xlab = "Entropy Similarity", 
     ylab = "Frequency",
     col = "lightpink",
     border = "black")
