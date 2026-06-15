# new_get_wiki_data.R
#
# Description: This script retrieves and processes mass spectrometry data from MassWiki.
#              It includes functions for filtering results and querying spectral library hits.
#
# Author: Original by zyang2k
# Last Modified: 2025-01-14

#------------------------------------------------------------------------------
# Required Packages
#------------------------------------------------------------------------------
library(dplyr)     # For data manipulation
library(readr)     # For reading CSV files
library(httr)      # For HTTP requests
library(jsonlite)  # For JSON parsing


#------------------------------------------------------------------------------
# Function Definitions
#------------------------------------------------------------------------------

#' Filter MassWiki Results
#' 
#' This function filters mass spectrometry results based on specific criteria:
#' 1. Keeps only manually annotated entries
#' 2. Removes entries with names starting with 'yy' or 'zz'
#'
#' @param masswiki_result (data.frame): A data frame containing MassWiki results
#'        with columns including 'is_manual_annotated' and 'name'
#'
#' @return data.frame: Filtered results containing only manually annotated entries
#'         with valid names
#'
#' @examples
#' filtered_data <- filter_masswiki_results(raw_results)
filter_masswiki_results <- function(masswiki_result) {
  filtered_result <- masswiki_result %>%
    filter(is_manual_annotated == TRUE) %>%
    filter(!grepl("^(yy|zz)", name))
  
  return(filtered_result)
}

#' Get Spectrum Data from MassWiki API
#' 
#' Retrieves spectral data for given wiki IDs from the MassWiki API.
#' The function handles multiple wiki IDs and includes error handling for API requests.
#'
#' @param wiki_ids (vector): A vector of wiki IDs to query
#'
#' @return list: A nested list where each element corresponds to a wiki_id and contains:
#'        - reference_library: Results from reference library identity search
#'        - annotation_library: Results from annotation library identity search
#'        Returns NULL for entries where the API request failed
#'
#' @details
#' The function:
#' 1. Processes each wiki_id individually
#' 2. Skips empty or NA wiki_ids
#' 3. URL encodes the wiki_ids to handle special characters
#' 4. Makes GET requests to the MassWiki API
#' 5. Parses both reference and annotation library results
#'
#' @examples
#' results <- get_spectrum_data(c("wiki_id_1", "wiki_id_2"))
get_spectrum_data <- function(wiki_ids) {
  # Ensure wiki_ids is a vector
  wiki_ids <- as.vector(wiki_ids)
  # Initialize results list
  all_results <- list()
  
  # Process each wiki_id
  for (wiki_id in wiki_ids) {
    # Skip if wiki_id is NA or empty
    if (is.na(wiki_id) || wiki_id == "") {
      cat("Skipping empty or NA wiki_id\n")
      next
    }
    
    # URL encode the wiki_id to handle special characters like '/'
    encoded_wiki_id <- URLencode(wiki_id, reserved = TRUE)
    url <- paste0('https://masswiki.us-west-2.elasticbeanstalk.com/analysis/get_data?wiki_id=', encoded_wiki_id)
    
    tryCatch({
      response <- GET(url, accept("application/json"))
      
      if (status_code(response) == 200) {
        data <- content(response, "parsed", simplifyVector = TRUE)
        
        results <- list()
        
        # Extract reference library identity search results if they exist
        if (!is.null(data$analysis$reference_library$identity_search)) {
          results$reference_library <- data$analysis$reference_library$identity_search
        } else {
          results$reference_library <- NULL
        }
        
        # Extract annotation library identity search results if they exist
        if (!is.null(data$analysis$annotation_library$identity_search)) {
          results$annotation_library <- data$analysis$annotation_library$identity_search
        } else {
          results$annotation_library <- NULL
        }
        
        all_results[[wiki_id]] <- results
        
      } else {
        cat("Failed to get spectrum data for", wiki_id, ". Status code:", status_code(response), "\n")
        all_results[[wiki_id]] <- NULL
      }
    }, error = function(e) {
      cat("Error processing wiki_id:", wiki_id, "-", e$message, "\n")
      all_results[[wiki_id]] <- NULL
    })
  }
  
  return(all_results)
}

#------------------------------------------------------------------------------
# Usage Notes:
# 1. Ensure all required packages are installed
# 2. Verify that input CSV files exist in the data/ directory
# 3. Check internet connection for API requests
# 4. Consider implementing rate limiting for large numbers of API requests
# 5. Monitor API response times and implement appropriate timeout settings
#------------------------------------------------------------------------------