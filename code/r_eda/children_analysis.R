# children analysis
library(ggplot2)
library(dplyr)
library(tidyr)
library(msentropy)


#### build family tree ######
# fxn to build family trees
build_family_trees <- function(data) {
  # Create a copy of original data with all columns
  tree_data <- data %>%
    mutate(
      # Add parent column, which is the duplicate_of value
      parent_id = duplicate_of,
      # For rows that are not duplicates (parents), parent_id will be NA
      is_parent = is.na(duplicate_of)
    )
  
  # find root parent
  find_root_parent <- function(node_id, tree_data) {
    current_id <- node_id
    visited <- c()
    
    # Follow the chain of parents until we hit NA or a cycle
    while (!is.na(current_id) && !(current_id %in% visited)) {
      parent_id <- tree_data$parent_id[tree_data$id == current_id]
      if (length(parent_id) == 0 || is.na(parent_id)) break
      visited <- c(visited, current_id)
      current_id <- parent_id
    }
    
    return(current_id)
  }
  
  # Add root parent information
  tree_data <- tree_data %>%
    rowwise() %>%
    mutate(
      root_parent_id = find_root_parent(id, tree_data),
      # Calculate depth in tree (0 for roots, 1 for direct children, etc.)
      tree_depth = if_else(is_parent, 0L,
                           length(unique(na.omit(c(
                             id,
                             parent_id,
                             root_parent_id
                           )))) - 1L
      )
    ) %>%
    ungroup()
  
  # Add family size information
  tree_data <- tree_data %>%
    group_by(root_parent_id) %>%
    mutate(
      family_size = n(),
      is_root = id == root_parent_id
    ) %>%
    ungroup()
  
  return(tree_data)
}

### get children based on parent id #######

# get children data from a parent id
get_children <- function(parent_id, tree_data) {
  tree_data %>% filter(duplicate_of == !!parent_id)
}

### retrieve all tree data #######

# First build tree data with all the extra information
tree_data <- build_family_trees(data)


### check how many parents ######
# Count unique parents (unique values in duplicate_of column, excluding NA)
unique_parents <- unique(na.omit(orb_neg_data$duplicate_of))
length(unique_parents)


### check if a child has child ###
# Function to analyze child relationships
# Function to analyze child relationships with correct parent identification
check_child_levels <- function(tree_data) {
  # Get all children (those with a parent_id)
  children <- tree_data %>%
    filter(!is.na(parent_id)) %>%
    select(id, parent_id, tree_depth, family_size)
  
  # Identify children who have their own children
  # by checking if their id appears as someone else's parent_id
  children_as_parents <- children %>%
    filter(id %in% tree_data$parent_id)
  
  # Summary statistics
  summary_stats <- list(
    total_children = nrow(children),
    children_as_parents = nrow(children_as_parents),
    percentage_child_parents = round(100 * nrow(children_as_parents) / nrow(children), 2)
  )
  
  # Distribution of tree depths
  depth_dist <- children %>%
    count(tree_depth) %>%
    arrange(tree_depth)
  
  # Get details about children who are also parents
  parent_child_details <- children_as_parents %>%
    mutate(
      num_own_children = sapply(id, function(pid) sum(tree_data$parent_id == pid, na.rm = TRUE))
    ) %>%
    arrange(desc(num_own_children))
  
  return(list(
    summary = summary_stats,
    depth_distribution = depth_dist,
    children_who_are_parents = parent_child_details
  ))
}

child_results <- check_child_levels(tree_data)

# Print results
print("Child Summary:")
print(child_results$summary)

print("\nTree Depth Distribution:")
print(child_results$depth_distribution)

# To check specific child
check_specific_child <- function(tree_data, child_id) {
  # Get the child's info
  child <- tree_data %>%
    filter(id == child_id) %>%
    select(id, parent_id, is_parent, tree_depth, family_size)
  
  # Get children of this child (if any)
  grandchildren <- tree_data %>%
    filter(parent_id == child_id) %>%
    select(id, parent_id, is_parent, tree_depth, family_size)
  
  list(
    child_info = child,
    has_children = nrow(grandchildren) > 0,
    number_of_children = nrow(grandchildren),
    grandchildren = grandchildren
  )
}

# Get overall statistics
child_results <- check_child_levels(tree_data)

child_results

##### seems like some child has their own child ######

# Investigate tree structure discrepancy
investigate_tree <- function(tree_data) {
  # Look at some depth 2 examples
  depth_2_examples <- tree_data %>%
    filter(tree_depth == 2) %>%
    select(id, parent_id, is_parent, tree_depth, family_size, root_parent_id) %>%
    head(5)
  
  # Get their immediate parents
  parent_ids <- unique(depth_2_examples$parent_id)
  immediate_parents <- tree_data %>%
    filter(id %in% parent_ids) %>%
    select(id, parent_id, is_parent, tree_depth, family_size, root_parent_id)
  
  # Check root parents
  root_parent_ids <- unique(depth_2_examples$root_parent_id)
  root_parents <- tree_data %>%
    filter(id %in% root_parent_ids) %>%
    select(id, parent_id, is_parent, tree_depth, family_size)
  
  # Get counts for different relationships
  relationship_summary <- tree_data %>%
    summarise(
      total_rows = n(),
      has_parent_id = sum(!is.na(parent_id)),
      is_marked_parent = sum(is_parent == TRUE),
      has_children = n_distinct(parent_id[!is.na(parent_id)]),
      depth_0 = sum(tree_depth == 0),
      depth_1 = sum(tree_depth == 1),
      depth_2 = sum(tree_depth == 2),
      unique_root_parents = n_distinct(root_parent_id[!is.na(root_parent_id)])
    )
  
  return(list(
    depth_2_samples = depth_2_examples,
    their_parents = immediate_parents,
    root_parents = root_parents,
    relationship_counts = relationship_summary
  ))
}

# Run investigation
tree_investigation <- investigate_tree(tree_data)

# Print results
print("Relationship Summary:")
print(tree_investigation$relationship_counts)

print("\nSample Depth 2 Nodes:")
print(tree_investigation$depth_2_samples)

print("\nTheir Immediate Parents:")
print(tree_investigation$their_parents)

print("\nRoot Parents:")
print(tree_investigation$root_parents)




# Function to trace a complete family tree branch
trace_tree_branch <- function(tree_data, root_id) {
  # Get root node
  root <- tree_data %>%
    filter(id == root_id) %>%
    select(id, parent_id, is_parent, tree_depth, family_size, sample)
  
  # Get all descendants
  descendants <- tree_data %>%
    filter(root_parent_id == root_id) %>%
    select(id, parent_id, is_parent, tree_depth, family_size, sample) %>%
    arrange(tree_depth)
  
  # Get depth distribution in this family
  depth_dist <- descendants %>%
    count(tree_depth) %>%
    arrange(tree_depth)
  
  # Check intermediate nodes (depth 1 nodes that have children)
  intermediate_nodes <- descendants %>%
    filter(tree_depth == 1) %>%
    mutate(
      has_children = id %in% descendants$parent_id
    )
  
  list(
    root = root,
    descendants = descendants,
    depth_distribution = depth_dist,
    intermediate_summary = intermediate_nodes %>%
      summarise(
        total_depth_1 = n(),
        depth_1_with_children = sum(has_children),
        depth_1_marked_parent = sum(is_parent)
      ),
    intermediate_nodes = intermediate_nodes
  )
}

# Let's examine one of the larger families we saw
# Using root_id 23892019 which had family_size 355
example_trace <- trace_tree_branch(tree_data, 23892019)

# Print results
print("Root node:")
print(example_trace$root)

print("\nDepth distribution in this family:")
print(example_trace$depth_distribution)

print("\nIntermediate nodes summary:")
print(example_trace$intermediate_summary)

# Look at some intermediate nodes that have children
print("\nSample of depth 1 nodes that have children:")
print(example_trace$intermediate_nodes %>% 
        filter(has_children) %>% 
        head(5))


### Revised fxn for direct family analysis ####
# Function to analyze complete family hierarchies
analyze_family_hierarchies <- function(tree_data) {
  # Get all root parents (those who have grandchildren)
  grandparents <- tree_data %>%
    filter(is.na(parent_id)) %>%  # Root level
    filter(id %in% (tree_data %>% 
                      filter(tree_depth == 2) %>% 
                      pull(root_parent_id))) %>%
    select(id, family_size)
  
  # Analyze each level of the hierarchy
  hierarchy_analysis <- tree_data %>%
    mutate(
      is_grandparent = id %in% grandparents$id,
      is_parent_only = id %in% tree_data$parent_id & !is_grandparent,
      has_own_children = id %in% tree_data$parent_id
    ) %>%
    summarise(
      total_grandparents = sum(is_grandparent),
      total_parents = sum(has_own_children),
      parents_who_are_children = sum(has_own_children & !is.na(parent_id)),
      total_children = sum(!is.na(parent_id)),
      total_grandchildren = sum(tree_depth == 2)
    )
  
  # Analyze family sizes and structures
  family_structures <- tree_data %>%
    group_by(root_parent_id) %>%
    summarise(
      total_members = n(),
      direct_children = sum(tree_depth == 1),
      grandchildren = sum(tree_depth == 2),
      has_complete_family = any(tree_depth == 2)
    ) %>%
    filter(!is.na(root_parent_id)) %>%
    ungroup()
  
  # Get detailed stats about parents who are also children
  parent_children_details <- tree_data %>%
    filter(!is.na(parent_id)) %>%  # They are children
    filter(id %in% tree_data$parent_id) %>%  # They are also parents
    select(id, parent_id, root_parent_id, tree_depth) %>%
    mutate(
      num_own_children = sapply(id, function(pid) sum(tree_data$parent_id == pid, na.rm = TRUE))
    )
  
  return(list(
    summary = hierarchy_analysis,
    family_structures = family_structures,
    parent_children = parent_children_details
  ))
}

# Function to trace a complete family line
trace_family_line <- function(tree_data, root_id) {
  # Get grandparent info
  grandparent <- tree_data %>%
    filter(id == root_id) %>%
    select(id, family_size, sample)
  
  # Get all direct children (parents)
  parents <- tree_data %>%
    filter(parent_id == root_id) %>%
    select(id, parent_id, tree_depth, family_size, sample) %>%
    mutate(
      num_own_children = sapply(id, function(pid) sum(tree_data$parent_id == pid, na.rm = TRUE))
    )
  
  # Get all grandchildren
  grandchildren <- tree_data %>%
    filter(root_parent_id == root_id, tree_depth == 2) %>%
    select(id, parent_id, tree_depth, family_size, sample)
  
  # Create family summary
  family_summary <- list(
    total_family_size = nrow(parents) + nrow(grandchildren) + 1,
    num_children = nrow(parents),
    num_grandchildren = nrow(grandchildren),
    children_who_are_parents = sum(parents$num_own_children > 0),
    avg_grandchildren_per_parent = round(nrow(grandchildren) / sum(parents$num_own_children > 0), 2)
  )
  
  return(list(
    grandparent = grandparent,
    parents = parents,
    grandchildren = grandchildren,
    summary = family_summary
  ))
}

# Example usage:
hierarchy_results <- analyze_family_hierarchies(tree_data)
print("Complete Family Hierarchy Summary:")
print(hierarchy_results$summary)

print("\nFamily Structures Overview:")
print(head(hierarchy_results$family_structures))

print("\nDetails of Parents who are Children:")
print(head(hierarchy_results$parent_children))

# To trace a specific family line:
# family_line <- trace_family_line(tree_data, root_id = 23892019)
# print("Family Line Summary:")
# print(family_line$summary)
# print("\nParents in this family:")
# print(family_line$parents %>% arrange(desc(num_own_children)))

# To trace a specific family:
# family_example <- trace_direct_family(tree_data, parent_id = 23892019)
# print("Parent info:")
# print(family_example$parent)
# print("\nDirect children summary:")
# print(family_example$children_summary)

### direct family similarity distribution ####

# Function to calculate similarities for direct families
# Function to calculate similarities for direct families with optimization
analyze_direct_family_similarities <- function(tree_data, batch_size = 100) {
  # Get direct parents
  direct_parents <- tree_data %>%
    filter(is.na(parent_id)) %>%  # Root level parents
    filter(id %in% tree_data$parent_id) %>%  # Must have children
    distinct()
  
  print(paste("Found", nrow(direct_parents), "direct parents"))
  
  # Process in batches
  total_batches <- ceiling(nrow(direct_parents) / batch_size)
  
  all_similarities <- list()
  
  for(batch in 1:total_batches) {
    start_idx <- (batch - 1) * batch_size + 1
    end_idx <- min(batch * batch_size, nrow(direct_parents))
    
    print(paste("Processing batch", batch, "of", total_batches, 
                "(parents", start_idx, "to", end_idx, ")"))
    
    # Process this batch of parents
    batch_parents <- direct_parents[start_idx:end_idx, ]
    
    batch_similarities <- map_df(1:nrow(batch_parents), function(i) {
      parent_row <- batch_parents[i,]
      
      # Get direct children
      direct_children <- tree_data %>% 
        filter(parent_id == parent_row$id)
      
      if(nrow(direct_children) == 0 || is.na(parent_row$msms)) return(NULL)
      
      parent_peaks <- convert_msms_to_matrix(parent_row$msms)
      if(is.null(parent_peaks)) return(NULL)
      
      # Process all children for this parent
      children_similarities <- direct_children %>%
        rowwise() %>%
        mutate(
          parent_id = parent_row$id,
          parent_sample = parent_row$sample,
          entropy_similarity_to_parent = {
            child_peaks <- convert_msms_to_matrix(msms)
            if(is.null(child_peaks)) {
              NA_real_
            } else {
              tryCatch({
                calculate_entropy_similarity(
                  peaks_a = parent_peaks,
                  peaks_b = child_peaks,
                  ms2_tolerance_in_da = 0.02,
                  ms2_tolerance_in_ppm = -1,
                  clean_spectra = TRUE,
                  min_mz = 0,
                  max_mz = 1000,
                  noise_threshold = 0.01,
                  max_peak_num = 100
                )
              }, error = function(e) {
                NA_real_
              })
            }
          }
        ) %>%
        ungroup()
      
      children_similarities
    })
    
    all_similarities[[batch]] <- batch_similarities
    
    # Optional: save intermediate results
    if(batch %% 5 == 0) {
      saveRDS(do.call(rbind, all_similarities), 
              file = paste0("family_similarities_batch_", batch, ".rds"))
    }
  }
  
  # Combine all batches
  direct_family_similarities <- do.call(rbind, all_similarities)
  
  # Calculate family-level statistics
  family_stats <- direct_family_similarities %>%
    group_by(parent_id) %>%
    summarise(
      n_children = n(),
      mean_similarity = mean(entropy_similarity_to_parent, na.rm = TRUE),
      median_similarity = median(entropy_similarity_to_parent, na.rm = TRUE),
      min_similarity = min(entropy_similarity_to_parent, na.rm = TRUE),
      max_similarity = max(entropy_similarity_to_parent, na.rm = TRUE),
      similarity_range = max_similarity - min_similarity,
      same_sample_children = sum(sample == parent_sample, na.rm = TRUE)
    ) %>%
    ungroup()
  
  # Create visualizations
  # 1. Overall similarity distribution by family size
  p1 <- ggplot(direct_family_similarities, 
               aes(x = entropy_similarity_to_parent)) +
    geom_histogram(binwidth = 0.05) +
    facet_wrap(~parent_id, scales = "free_y") +
    theme_minimal() +
    labs(title = "Similarity Distribution Within Each Direct Family",
         x = "Entropy Similarity",
         y = "Count")
  
  # 2. Family size vs average similarity
  p2 <- ggplot(family_stats, 
               aes(x = n_children, y = mean_similarity)) +
    geom_point() +
    geom_smooth(method = "lm", se = TRUE) +
    theme_minimal() +
    labs(title = "Family Size vs Average Similarity",
         x = "Number of Direct Children",
         y = "Mean Similarity")
  
  # 3. Similarity range distribution
  p3 <- ggplot(family_stats, 
               aes(x = similarity_range)) +
    geom_histogram(binwidth = 0.05) +
    theme_minimal() +
    labs(title = "Distribution of Similarity Ranges within Families",
         x = "Similarity Range (Max - Min)",
         y = "Count")
  
  # Calculate summary statistics
  summary_stats <- list(
    total_direct_families = nrow(family_stats),
    mean_family_size = mean(family_stats$n_children),
    median_family_size = median(family_stats$n_children),
    overall_mean_similarity = mean(direct_family_similarities$entropy_similarity_to_parent, 
                                   na.rm = TRUE),
    mean_similarity_range = mean(family_stats$similarity_range, na.rm = TRUE),
    same_sample_pairs = sum(direct_family_similarities$sample == 
                              direct_family_similarities$parent_sample, 
                            na.rm = TRUE)
  )
  
  return(list(
    summary = summary_stats,
    family_stats = family_stats,
    similarity_data = direct_family_similarities,
    family_distributions = p1,
    size_vs_similarity = p2,
    range_distribution = p3
  ))
}



# Example usage with smaller batch size:
direct_family_results <- analyze_direct_family_similarities(tree_data, batch_size = 50)


# Create basic visualizations for similarity distributions
plot_similarity_distributions <- function(similarity_data) {
  # Set up a 2x2 plotting layout
  par(mfrow=c(2,2), mar=c(5,5,4,2))
  
  # 1. Overall histogram of similarity scores
  hist(similarity_data$entropy_similarity_to_parent, 
       breaks=30,
       main="Distribution of Similarity Scores",
       xlab="Entropy Similarity Score",
       ylab="Frequency",
       col="lightblue",
       border="white")
  abline(v=mean(similarity_data$entropy_similarity_to_parent, na.rm=TRUE), 
         col="red", lwd=2, lty=2)
  legend("topright", 
         legend=c("Mean Similarity"), 
         col=c("red"), 
         lty=2, 
         lwd=2)
  
  # 2. Box plot comparing same vs different sample similarities
  boxplot(entropy_similarity_to_parent ~ (sample == parent_sample), 
          data=similarity_data,
          main="Similarity by Sample Type",
          xlab="Same Sample",
          ylab="Entropy Similarity Score",
          col=c("lightblue", "lightgreen"),
          names=c("Different", "Same"))
  
  # 3. Density plot of similarities
  plot(density(similarity_data$entropy_similarity_to_parent, na.rm=TRUE),
       main="Density Plot of Similarity Scores",
       xlab="Entropy Similarity Score",
       ylab="Density",
       col="blue",
       lwd=2)
  
  # 4. QQ plot to check distribution
  qqnorm(similarity_data$entropy_similarity_to_parent,
         main="Normal Q-Q Plot of Similarities")
  qqline(similarity_data$entropy_similarity_to_parent, col="red")
  
  # Reset plotting parameters
  par(mfrow=c(1,1))
}

# Create summary plots for family-level statistics
plot_family_summaries <- function(family_stats) {
  # Set up a 2x2 plotting layout
  par(mfrow=c(2,2), mar=c(5,5,4,2))
  
  # 1. Histogram of family sizes
  hist(family_stats$n_children,
       breaks=30,
       main="Distribution of Family Sizes",
       xlab="Number of Children",
       ylab="Frequency",
       col="lightblue",
       border="white")
  
  # 2. Scatter plot of family size vs mean similarity
  plot(family_stats$n_children, 
       family_stats$mean_similarity,
       main="Family Size vs Mean Similarity",
       xlab="Number of Children",
       ylab="Mean Similarity Score",
       pch=19,
       col=rgb(0,0,1,0.2))
  abline(lm(mean_similarity ~ n_children, data=family_stats), 
         col="red", lwd=2)
  
  # 3. Histogram of similarity ranges within families
  hist(family_stats$similarity_range,
       breaks=30,
       main="Distribution of Within-Family\nSimilarity Ranges",
       xlab="Similarity Range (Max - Min)",
       ylab="Frequency",
       col="lightblue",
       border="white")
  
  # 4. Box plot of similarities by family size groups
  # Create size categories
  family_stats$size_category <- cut(family_stats$n_children, 
                                    breaks=c(0,2,5,10,20,max(family_stats$n_children)),
                                    labels=c("1-2","3-5","6-10","11-20","20+"))
  
  boxplot(mean_similarity ~ size_category,
          data=family_stats,
          main="Similarity by Family Size Category",
          xlab="Family Size",
          ylab="Mean Similarity Score",
          col="lightblue")
  
  # Reset plotting parameters
  par(mfrow=c(1,1))
}

# Example usage:
plot_similarity_distributions(direct_family_results$similarity_data)
plot_family_summaries(direct_family_results$family_stats)

# To save plots to PDF:
# pdf("similarity_distributions.pdf", width=12, height=10)
# plot_similarity_distributions(direct_family_results$similarity_data)
# dev.off()
# 
# pdf("family_summaries.pdf", width=12, height=10)
# plot_family_summaries(direct_family_results$family_stats)
# dev.off()




# Function to analyze similarity differences between same and different sample pairs
analyze_sample_similarities <- function(similarity_data) {
  # Create groups
  same_sample <- similarity_data$entropy_similarity_to_parent[
    similarity_data$sample == similarity_data$parent_sample]
  diff_sample <- similarity_data$entropy_similarity_to_parent[
    similarity_data$sample != similarity_data$parent_sample]
  
  # Basic summary statistics
  summary_stats <- list(
    same_sample = list(
      n = length(same_sample),
      mean = mean(same_sample, na.rm = TRUE),
      sd = sd(same_sample, na.rm = TRUE),
      median = median(same_sample, na.rm = TRUE)
    ),
    diff_sample = list(
      n = length(diff_sample),
      mean = mean(diff_sample, na.rm = TRUE),
      sd = sd(diff_sample, na.rm = TRUE),
      median = median(diff_sample, na.rm = TRUE)
    )
  )
  
  # Check normality (Shapiro-Wilk test)
  # Only test on a sample if n > 5000 due to Shapiro-Wilk limitations
  normality_test <- list(
    same_sample = shapiro.test(if(length(same_sample) > 5000) 
      sample(same_sample, 5000) else same_sample),
    diff_sample = shapiro.test(if(length(diff_sample) > 5000) 
      sample(diff_sample, 5000) else diff_sample)
  )
  
  # Statistical tests
  # 1. Welch's t-test (doesn't assume equal variances)
  t_test_result <- t.test(same_sample, diff_sample)
  
  # 2. Non-parametric test (Mann-Whitney U test)
  wilcox_result <- wilcox.test(same_sample, diff_sample)
  
  # Effect size (Cohen's d)
  cohens_d <- (mean(same_sample, na.rm = TRUE) - mean(diff_sample, na.rm = TRUE)) /
    sqrt((var(same_sample, na.rm = TRUE) + var(diff_sample, na.rm = TRUE)) / 2)
  
  # Create visualization
  par(mfrow = c(1,2))
  
  # Boxplot
  boxplot(list(
    "Different Sample" = diff_sample,
    "Same Sample" = same_sample
  ), 
  main = "Similarity Scores by Sample Type",
  ylab = "Entropy Similarity Score",
  col = c("lightblue", "lightgreen"))
  
  # Density plot
  plot(density(diff_sample, na.rm = TRUE), 
       main = "Density Plot of Similarity Scores",
       xlab = "Entropy Similarity Score",
       ylab = "Density",
       col = "blue",
       lwd = 2)
  lines(density(same_sample, na.rm = TRUE), 
        col = "green",
        lwd = 2)
  legend("topright", 
         legend = c("Different Sample", "Same Sample"),
         col = c("blue", "green"),
         lwd = 2)
  
  # Reset plotting parameters
  par(mfrow = c(1,1))
  
  # Return all results
  return(list(
    summary = summary_stats,
    normality = normality_test,
    t_test = t_test_result,
    wilcox_test = wilcox_result,
    effect_size = cohens_d
  ))
}

# Example usage:
sample_analysis <- analyze_sample_similarities(direct_family_results$similarity_data)
# 
# # Print results
print("Summary Statistics:")
print(sample_analysis$summary)
print("\nNormality Tests:")
print(sample_analysis$normality)
print("\nt-test Results:")
print(sample_analysis$t_test)
print("\nWilcoxon Test Results:")
print(sample_analysis$wilcox_test)
print("\nEffect Size (Cohen's d):")
print(sample_analysis$effect_size)

### Function to convert MSMS string to matrix #####
convert_msms_to_matrix <- function(msms_string) {
  if(is.na(msms_string)) return(NULL)
  
  pairs <- strsplit(msms_string, " ")[[1]]
  mz_values <- numeric()
  intensity_values <- numeric()
  
  for(pair in pairs) {
    values <- as.numeric(strsplit(pair, ":")[[1]])
    if(length(values) == 2) {  # Ensure we have both mz and intensity
      mz_values <- c(mz_values, values[1])
      intensity_values <- c(intensity_values, values[2])
    }
  }
  
  if(length(mz_values) == 0) return(NULL)
  matrix(c(mz_values, intensity_values), ncol = 2, byrow = FALSE)
}

# Function to calculate similarities for a parent and its children
calculate_family_similarities <- function(parent_row, tree_data) {
  # Get children
  children <- tree_data %>% 
    filter(parent_id == parent_row$id)
  
  if(nrow(children) == 0 || is.na(parent_row$msms)) return(NULL)
  
  # Convert parent MSMS
  parent_peaks <- convert_msms_to_matrix(parent_row$msms)
  if(is.null(parent_peaks)) return(NULL)
  
  # Calculate similarities for each child
  children_with_sim <- children %>%
    rowwise() %>%
    mutate(
      entropy_similarity_to_parent = {
        child_peaks <- convert_msms_to_matrix(msms)
        if(is.null(child_peaks)) {
          NA_real_
        } else {
          tryCatch({
            calculate_entropy_similarity(
              peaks_a = parent_peaks,
              peaks_b = child_peaks,
              ms2_tolerance_in_da = 0.02,
              ms2_tolerance_in_ppm = -1,
              clean_spectra = TRUE,
              min_mz = 0,
              max_mz = 1000,
              noise_threshold = 0.01,
              max_peak_num = 100
            )
          }, error = function(e) {
            NA_real_
          })
        }
      },
      parent_sample = parent_row$sample
    ) %>%
    ungroup()
  
  return(children_with_sim)
}

# Main analysis function
analyze_tree_similarities <- function(tree_data) {
  # Get all parents
  parents <- tree_data %>%
    filter(is_parent == TRUE)
  
  print(paste("Found", nrow(parents), "parents"))
  
  # Calculate similarities for all families
  all_children <- map_df(1:nrow(parents), function(i) {
    result <- calculate_family_similarities(parents[i,], tree_data)
    if(i %% 10 == 0) print(paste("Processed", i, "parents"))
    result
  })
  
  print("Finished calculating similarities")
  
  # 1. Overall distribution of similarities
  p1 <- ggplot(all_children, aes(x = entropy_similarity_to_parent)) +
    geom_histogram(binwidth = 0.05) +
    theme_minimal() +
    labs(title = "Distribution of Parent-Child Spectral Similarities",
         x = "Entropy Similarity",
         y = "Count")
  
  # 2. Same-sample distribution
  same_sample_pairs <- all_children %>%
    filter(!is.na(parent_sample) & sample == parent_sample)
  
  p2 <- ggplot(same_sample_pairs, aes(x = entropy_similarity_to_parent)) +
    geom_histogram(binwidth = 0.05) +
    theme_minimal() +
    labs(title = "Distribution of Parent-Child Similarities (Same Sample)",
         x = "Entropy Similarity",
         y = "Count")
  
  # Summary statistics
  summary_stats <- list(
    total_parents = nrow(parents),
    total_children = nrow(all_children),
    same_sample_pairs = nrow(same_sample_pairs),
    different_sample_pairs = nrow(all_children) - nrow(same_sample_pairs),
    mean_similarity = mean(all_children$entropy_similarity_to_parent, na.rm = TRUE),
    median_similarity = median(all_children$entropy_similarity_to_parent, na.rm = TRUE),
    same_sample_mean_similarity = mean(same_sample_pairs$entropy_similarity_to_parent, na.rm = TRUE)
  )
  
  # Similarity summary by sample type
  similarity_summary <- all_children %>%
    mutate(same_sample = !is.na(parent_sample) & parent_sample == sample) %>%
    group_by(same_sample) %>%
    summarise(
      mean_similarity = mean(entropy_similarity_to_parent, na.rm = TRUE),
      median_similarity = median(entropy_similarity_to_parent, na.rm = TRUE),
      n = n()
    )
  
  return(list(
    summary = summary_stats,
    overall_similarity_plot = p1,
    same_sample_similarity_plot = p2,
    similarity_data = all_children,
    by_sample_type = similarity_summary
  ))
}

# Run the analysis on subset of data
results <- analyze_tree_similarities(tree_data %>% slice(1:100))

# Print summary
print("Summary Statistics:")
print(results$summary)

# Print sample type breakdown
print("\nSimilarity by Sample Type:")
print(results$by_sample_type)

# Display plots
print(results$overall_similarity_plot)
print(results$same_sample_similarity_plot)

###### end of 1223 #####

# Function to calculate similarity with parent for a family
calculate_family_similarities <- function(parent_id, children_data) {
  # Get parent MSMS
  parent_msms <- data$msms[data$id == parent_id]
  parent_peaks <- convert_msms_to_matrix(parent_msms)
  
  # Calculate similarities for each child
  similarities <- numeric(nrow(children_data))
  
  for(i in 1:nrow(children_data)) {
    if(!is.na(children_data$msms[i])) {
      tryCatch({
        child_peaks <- convert_msms_to_matrix(children_data$msms[i])
        
        similarities[i] <- calculate_entropy_similarity(
          peaks_a = parent_peaks,
          peaks_b = child_peaks,
          ms2_tolerance_in_da = 0.02,
          ms2_tolerance_in_ppm = -1,
          clean_spectra = TRUE,
          min_mz = 0,
          max_mz = 1000,
          noise_threshold = 0.01,
          max_peak_num = 100
        )
      }, error = function(e) {
        print(paste("Error in family", parent_id, "child", i, ":", e$message))
        similarities[i] <- NA
      })
    } else {
      similarities[i] <- NA
    }
  }
  return(similarities)
}
# Calculate similarities for all sampled families
family_similarities <- list()
names_for_families <- character()  # To store parent IDs

for(i in seq_along(sampled_families)) {
  parent_id <- random_parents[i]
  children_data <- sampled_families[[i]]
  
  # Skip if no children data
  if(nrow(children_data) == 0) {
    print(paste("Skipping parent", parent_id, "- no children"))
    next
  }
  
  similarities <- calculate_family_similarities(parent_id, children_data)
  
  # Store similarities with proper indexing
  family_similarities[[i]] <- similarities
  names_for_families[i] <- parent_id
}

# Name the list elements after storing all data
names(family_similarities) <- names_for_families

# Verify the data
print(paste("Number of families processed:", length(family_similarities)))
print("Sample sizes of first few families:")
print(sapply(head(family_similarities), length))


# Combine all similarities
all_similarities <- unlist(family_similarities)

# Create summary statistics for each family
family_summaries <- lapply(names(family_similarities), function(parent_id) {
  sims <- family_similarities[[parent_id]]
  data.frame(
    parent_id = parent_id,
    mean_similarity = mean(sims, na.rm = TRUE),
    median_similarity = median(sims, na.rm = TRUE),
    sd_similarity = sd(sims, na.rm = TRUE),
    n_children = length(sims),
    q25 = quantile(sims, 0.25, na.rm = TRUE),
    q75 = quantile(sims, 0.75, na.rm = TRUE)
  )
}) %>% bind_rows()

# Overall distribution plot
p1 <- ggplot(data.frame(similarity = na.omit(all_similarities)), aes(x = similarity)) +
  geom_histogram(bins = 50, fill = "skyblue", color = "black", alpha = 0.7) +
  geom_vline(aes(xintercept = mean(na.omit(all_similarities))), 
             color = "red", linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "Overall Distribution of Parent-Child Similarity Scores",
    subtitle = paste("Mean =", round(mean(na.omit(all_similarities)), 3)),
    x = "Similarity Score",
    y = "Count"
  )

# Distribution of family median similarities
p2 <- ggplot(family_summaries, aes(x = median_similarity)) +
  geom_histogram(bins = 30, fill = "lightgreen", color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "Distribution of Median Similarities Across Families",
    x = "Median Family Similarity",
    y = "Count of Families"
  )

# Relationship between family size and similarity
p3 <- ggplot(family_summaries, aes(x = n_children, y = median_similarity)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = TRUE) +
  scale_x_log10() +  # log scale for better visualization
  theme_minimal() +
  labs(
    title = "Median Similarity vs Family Size",
    x = "Number of Children (log scale)",
    y = "Median Similarity"
  )

print(p1)
print(p2)
print(p3)

# Print summary statistics
print("Overall similarity statistics:")
print(summary(all_similarities))

print("\nFamily-level statistics:")
print(summary(family_summaries[, c("mean_similarity", "median_similarity", "sd_similarity", "n_children")]))




######### example, parent = 23359501 #########
children <- get_children("23359501", tree_data)
# 1. Distribution of retention index
hist(children$retention_index, 
     main="Distribution of Retention Index",
     xlab="Retention Index",
     col="skyblue",
     breaks=30)

# 2. Distribution of accurate mass
hist(children$accurate_mass,
     main="Distribution of Accurate Mass",
     xlab="Accurate Mass",
     col="lightgreen",
     breaks=30)

# 3. Distribution of p intensity
hist(children$pre_cursors_intensity,
     main="Distribution of Precursor Intensity",
     xlab="Precursor Intensity",
     col="lightgreen",
     breaks=30)

# 4. Distribution of pct standard found
hist(children$correction_standards_found_in_percent,
     main="Distribution of pct standard found",
     xlab="pct standard found",
     col="lightgreen",
     breaks=20)


# Basic summary statistics
summary_stats <- summary(children[c("retention_index", "accurate_mass", "normalized_entropy")])
print(summary_stats)

# Count of some interesting varibales
table(children$target_type)





# pairwise calculate spectrum similarity to parent msms

# convert MSMS string into a matrix
convert_msms_to_matrix <- function(msms_string) {
  pairs <- strsplit(msms_string, " ")[[1]]
  mz_values <- numeric()
  intensity_values <- numeric()
  
  for(pair in pairs) {
    values <- as.numeric(strsplit(pair, ":")[[1]])
    mz_values <- c(mz_values, values[1])
    intensity_values <- c(intensity_values, values[2])
  }
  
  matrix(c(mz_values, intensity_values), ncol = 2, byrow = FALSE)
}

# retreive parent info
parent_id <- "23359501"
parent_msms <- data$msms[data$id == parent_id]
parent_peaks <- convert_msms_to_matrix(parent_msms)

# calculate similarities
similarities <- numeric(nrow(children))
for(i in 1:nrow(children)) {
  if(!is.na(children$msms[i])) {
    tryCatch({
      child_peaks <- convert_msms_to_matrix(children$msms[i])
      
      # Calculate similarity with all parameters
      similarities[i] <- calculate_entropy_similarity(
        peaks_a = parent_peaks,
        peaks_b = child_peaks,
        ms2_tolerance_in_da = 0.02,
        ms2_tolerance_in_ppm = -1,  # disable ppm tolerance
        clean_spectra = TRUE,
        min_mz = 0,
        max_mz = 1000,
        noise_threshold = 0.01,
        max_peak_num = 100
      )
    }, error = function(e) {
      print(paste("Error in row", i, ":", e$message))
      similarities[i] <- NA
    })
  } else {
    similarities[i] <- NA
  }
}

# Add similarities to children dataframe
children$entropy_similarity_to_parent <- similarities

# Print summary statistics
print(summary(similarities))

# Plot distribution
hist(similarities[!is.na(similarities)], 
     main="Distribution of Entropy Similarities to Parent",
     xlab="Entropy Similarity",
     ylab="Count",
     col="lightblue",
     breaks=20)


################# db level parent-children analysis ##################

# Number of children per parent
children_counts <- table(data$duplicate_of)

median(children_counts)
mean(children_counts)
max(children_counts)
min(children_counts)
length(children_counts)

hist(children_counts)
# Create bins for the children counts
children_bins <- cut(as.numeric(children_counts), 
                     breaks = c(0, 1, 2, 5, 10, 50, 100, 500, 1000, 4000),
                     labels = c("1", "2", "3-5", "6-10", "11-50", 
                                "51-100", "101-500", "501-1000", ">1000"),
                     include.lowest = TRUE)

# Create summary dataframe
bin_summary <- data.frame(
  bin = children_bins
) %>%
  group_by(bin) %>%
  summarise(
    count = n(),
    pct = n()/length(children_counts) * 100
  )

# Print summary
print("Distribution of parents by number of children:")
print(bin_summary)

# Visualizations
# 1. Distribution plot (log scale)
p1 <- ggplot(data.frame(children = as.numeric(children_counts)), aes(x = children)) +
  geom_histogram(bins = 50) +
  scale_x_log10() +
  theme_minimal() +
  labs(
    title = "Distribution of Children per Parent (Log Scale)",
    x = "Number of Children (log scale)",
    y = "Count of Parents"
  )

# 2. Bar plot of binned data
p2 <- ggplot(bin_summary, aes(x = bin, y = count)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Number of Parents by Children Count Bins",
    x = "Number of Children",
    y = "Count of Parents"
  )



print(p1)
print(p2)


