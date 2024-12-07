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



################## use wiki_id to query lib and annotation####################
# Function to get spectrum data from MassWiki API
get_spectrum_data <- function(wiki_id) {
  url <- paste0('https://masswiki.us-west-2.elasticbeanstalk.com/analysis/get_data?wiki_id=', wiki_id)
  response <- GET(url, accept_json())
  
  if (status_code(response) == 200) {
    return(content(response, "parsed", simplifyVector = TRUE))
  } else {
    cat("Failed to get spectrum data for", wiki_id, ". Status code:", status_code(response), "\n")
    return(NULL)
  }
}

# Function to extract annotation and library lists from spectrum data
extract_annotation_and_library_lists <- function(spectrum_data) {
  annotation_list <- c()
  library_list <- c()
  
  if (!is.null(spectrum_data$analysis)) {
    analysis_data <- spectrum_data$analysis
    
    # Extract annotation search IDs
    if (!is.null(analysis_data$annotation_search$identity_search)) {
      annotation_list <- sapply(analysis_data$annotation_search$identity_search, function(x) x[1])
      # Filter out invalid IDs
      annotation_list <- annotation_list[nchar(annotation_list) >= 3 & !is.na(as.numeric(annotation_list)) == FALSE]
    }
    
    # Extract library search IDs
    if (!is.null(analysis_data$library_search$identity_search)) {
      library_list <- sapply(analysis_data$library_search$identity_search, function(x) x[1])
      # Filter out invalid IDs
      library_list <- library_list[nchar(library_list) >= 3 & !is.na(as.numeric(library_list)) == FALSE]
    }
  }
  
  return(list(annotation_list = annotation_list, library_list = library_list))
}

# Function to get annotation data from MassWiki API
get_annotation_data <- function(annotation_list) {
  url <- 'https://masswiki.us-west-2.elasticbeanstalk.com/get/masswiki_data'
  payload <- list(
    wiki_id_list = annotation_list,
    type = "annotation",
    get_details = FALSE
  )
  response <- POST(url, body = toJSON(payload, auto_unbox = TRUE), content_type_json(), accept_json())
  
  if (status_code(response) == 200) {
    return(content(response, "parsed", simplifyVector = TRUE))
  } else {
    cat("Failed to get annotation data. Status code:", status_code(response), "\n")
    return(NULL)
  }
}

# Function to get spectral library data from MassWiki API
get_spectral_library_data <- function(library_list) {
  url <- 'https://masswiki.us-west-2.elasticbeanstalk.com/get/masswiki_data'
  payload <- list(
    wiki_id_list = library_list,
    type = "spectral_library",
    get_details = FALSE
  )
  
  response <- POST(url, body = toJSON(payload, auto_unbox = TRUE), content_type_json(), accept_json())
  
  if (status_code(response) == 200) {
    content_data <- content(response, "parsed", simplifyVector = TRUE)
    if (length(content_data) == 0) {
      cat("No data found in library search response.\n")
      return(NULL)
    } else {
      return(content_data)
    }
  } else {
    cat("Failed to get spectral library data. Status code:", status_code(response), "\n")
    cat("Response content: ", content(response, "text"), "\n")
    return(NULL)
  }
}

# Main function to retrieve data and organize into list of DataFrames
retrieve_and_organize_data <- function(wiki_ids) {
  wiki_data <- list()
  
  for (wiki_id in wiki_ids) {
    # Fetch spectrum data from the API
    spectrum_data <- get_spectrum_data(wiki_id)
    
    # If valid spectrum data is returned, extract the annotation and library lists
    if (!is.null(spectrum_data)) {
      lists <- extract_annotation_and_library_lists(spectrum_data)
      annotation_list <- lists$annotation_list
      library_list <- lists$library_list
      
      # Fetch annotation data using the extracted annotation list
      if (length(annotation_list) > 0) {
        annotation_data <- get_annotation_data(annotation_list)
        if (!is.null(annotation_data)) {
          annotation_df <- as.data.frame(annotation_data)
        } else {
          annotation_df <- data.frame()
        }
      } else {
        annotation_df <- data.frame()
      }
      
      # Fetch spectral library data using the extracted library list
      if (length(library_list) > 0) {
        library_data <- get_spectral_library_data(library_list)
        if (!is.null(library_data)) {
          library_df <- as.data.frame(library_data)
          # Calculate delta m/z for both DataFrames
          if (!is.null(spectrum_data$spectrum$precursor_mz)) {
            precursor_mz <- spectrum_data$spectrum$precursor_mz
            if (!is.null(annotation_df$precursor_mz)) {
              annotation_df$delta_mz <- annotation_df$precursor_mz - precursor_mz
            }
            if (!is.null(library_df$precursor_mz)) {
              library_df$delta_mz <- library_df$precursor_mz - precursor_mz
            }
          }
        } else {
          library_df <- data.frame()
        }
      } else {
        library_df <- data.frame()
      }
      
      # Store the DataFrames in the list
      wiki_data[[wiki_id]] <- list(
        annotation_search = annotation_df,
        library_search = library_df
      )
    } else {
      cat("No valid spectrum data found for", wiki_id, "\n")
      wiki_data[[wiki_id]] <- list(
        annotation_search = data.frame(),
        library_search = data.frame()
      )
    }
  }
  
  return(wiki_data)
}


wiki_ids <- filtered_masswiki_result$wiki_id

# Retrieve and organize data into list of DataFrames
wiki_data <- retrieve_and_organize_data(wiki_ids)



######## Extract annotation_search and library_search data from wiki_data######## 

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

# Combine the DataFrames into one
annotation_search_combined <- bind_rows(annotation_search_list)
library_search_combined <- do.call(rbind, library_search_list)



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

