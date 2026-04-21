library(dplyr)
library(readr)
library(httr)
library(jsonlite)

# 12062024
################## Load masswiki result.####################
masswiki_result <- read_csv("data/masswiki_result.csv")


################## filter masswiki result.####################
# filter out those without user annotated names
# "user_annotation-name"  == NA were filtered

filtered_masswiki_result <- masswiki_result %>%
  filter(!is.na(`user_annotation-name`))

# filter out those start with yy or zz
filtered_masswiki_result <- filtered_masswiki_result %>%
  filter(!grepl("^(yy|zz)", `user_annotation-name`))


################## extract library hits ####################
# Function to extract annotation and library lists with identity scores
# Function to extract annotation and library lists with identity scores
extract_annotation_and_library_lists <- function(spectrum_data) {
  annotation_list <- c()
  library_list <- c()
  annotation_scores <- c()
  library_scores <- c()
  
  if (!is.null(spectrum_data$analysis)) {
    analysis_data <- spectrum_data$analysis
    
    # Extract annotation search IDs and identity scores
    if (!is.null(analysis_data$annotation_search$identity_search)) {
      print("Annotation identity_search:")
      print(analysis_data$annotation_search$identity_search)
      
      # Extract annotation IDs and scores
      annotation_list <- sapply(analysis_data$annotation_search$identity_search, function(x) x[1])
      annotation_scores <- sapply(analysis_data$annotation_search$identity_search, function(x) as.numeric(x[2]))  # Ensure the score is numeric
      
      # Debugging: Check lengths before filtering
      print("Before filtering:")
      print(length(annotation_list))
      print(length(annotation_scores))
      
      # Filter out invalid IDs and their corresponding scores
      valid_annotation_indices <- nchar(annotation_list) >= 3 & !is.na(as.numeric(annotation_list))
      annotation_list <- annotation_list[valid_annotation_indices]
      annotation_scores <- annotation_scores[valid_annotation_indices]
      
      # Debugging: Check lengths after filtering
      print("After filtering:")
      print(length(annotation_list))
      print(length(annotation_scores))
    }
    
    # Extract library search IDs and identity scores
    if (!is.null(analysis_data$library_search$identity_search)) {
      print("Library identity_search:")
      print(analysis_data$library_search$identity_search)
      
      # Extract library IDs and scores
      library_list <- sapply(analysis_data$library_search$identity_search, function(x) x[1])
      library_scores <- sapply(analysis_data$library_search$identity_search, function(x) as.numeric(x[2]))  # Ensure the score is numeric
      
      # Filter out invalid IDs and their corresponding scores
      valid_library_indices <- nchar(library_list) >= 3 & !is.na(as.numeric(library_list))
      library_list <- library_list[valid_library_indices]
      library_scores <- library_scores[valid_library_indices]
      
      # Debugging: Check lengths after filtering
      print("After filtering library:")
      print(length(library_list))
      print(length(library_scores))
    }
  }
  
  return(list(
    annotation_list = annotation_list,
    annotation_scores = annotation_scores,
    library_list = library_list,
    library_scores = library_scores
  ))
}

# Retrieve and organize data into list of DataFrames
wiki_data_score <- retrieve_and_organize_data(filtered_masswiki_result$wiki_id[1])

# Check the result for a specific wiki_id
print(wiki_data$`aZ9CRE8/HF1KSXBS`$annotation_search)
print(wiki_data$`aZ9CRE8/HF1KSXBS`$library_search)

# Extract annotation_search and library_search data from wiki_data
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

# Combine the DataFrames into one
annotation_search_combined <- bind_rows(annotation_search_list, .id = "wiki_id")
library_search_combined <- bind_rows(library_search_list, .id = "wiki_id")

# Check the combined result
print(head(annotation_search_combined))
print(head(library_search_combined))

######################## each wiki id metadata####################


fetch_wiki_data <- function(wiki_id) {
  # Encode the wiki_id by replacing slashes with %2F
  encoded_id <- gsub("/", "%2F", wiki_id)
  
  # Construct the URL with the encoded wiki_id
  url <- paste0('https://masswiki.us-west-2.elasticbeanstalk.com/analysis/get_data?wiki_id=', encoded_id)
  
  # Send GET request
  response <- GET(url, add_headers('accept' = 'application/json'))
  
  # Check if the request was successful (status code 200)
  if (status_code(response) == 200) {
    # Parse and return the JSON content from the response
    return(content(response, "parsed"))
  } else {
    # Handle errors or non-200 responses
    message("Failed to fetch data for wiki_id: ", wiki_id)
    return(NULL)
  }
}

# fetch wiki data for each wiki id
results <- lapply(wiki_ids, fetch_wiki_data)

