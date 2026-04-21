#' @title Mass Spectrometry Data Analysis Package
#' @description A package for analyzing mass spectrometry data and building hierarchical
#'   trees from duplicate relationships. Supports MS/MS spectral similarity calculations
#'   and tree-based organization of related spectra.
#' 
#' @section Package Dependencies:
#'   This package requires the following R packages:
#'   * dplyr - For data manipulation
#'   * msentropy - For calculating MS/MS spectral similarities
#'
#' @importFrom dplyr %>% mutate filter select rename inner_join
#' @importFrom msentropy calculate_entropy_similarity
#' @export

#' Build Hierarchical Trees from Duplicate Relationships
#' 
#' Takes a data frame containing duplicate relationships between entries and builds
#' complete hierarchical trees, including metadata like tree depth, root parents,
#' and family size.
#'
#' @param data A data frame containing at least:
#'   \describe{
#'     \item{id}{Unique identifier for each entry}
#'     \item{duplicate_of}{ID of the parent entry this is a duplicate of (NA for roots)}
#'   }
#'
#' @return A data frame with additional columns:
#'   \describe{
#'     \item{parent_id}{Same as duplicate_of, for clarity}
#'     \item{is_parent}{Logical indicating if entry is a root (has no parent)}
#'     \item{root_parent_id}{ID of the ultimate root parent}
#'     \item{tree_depth}{Integer depth in tree (0 for roots)}
#'     \item{family_size}{Total number of entries in this tree}
#'     \item{is_root}{Logical indicating if this is the root entry}
#'   }
#'
#' @examples
#' data <- data.frame(
#'   id = c(1, 2, 3),
#'   duplicate_of = c(NA, 1, 1)
#' )
#' tree_data <- build_family_trees(data)
#'
#' @export
build_family_trees <- function(data) {
  # Create copy of data with parent relationships
  tree_data <- data %>%
    mutate(
      parent_id = duplicate_of,
      is_parent = is.na(duplicate_of)
    )
  
  # Helper function to find root parent with cycle detection
  find_root_parent <- function(node_id, tree_data, cache = new.env()) {
    # Return cached result if available
    if (exists(as.character(node_id), cache)) {
      return(get(as.character(node_id), cache))
    }
    
    current_id <- node_id
    visited <- c()
    
    # Follow parent chain until root or cycle detected
    while (!is.na(current_id) && !(current_id %in% visited)) {
      parent_id <- tree_data$parent_id[tree_data$id == current_id]
      if (length(parent_id) == 0 || is.na(parent_id)) break
      visited <- c(visited, current_id) 
      current_id <- parent_id
    }
    
    # Cache and return result
    assign(as.character(node_id), current_id, cache)
    return(current_id)
  }
  
  # Add tree metadata
  tree_data <- tree_data %>%
    rowwise() %>%
    mutate(
      # Find root parent for each node
      root_parent_id = find_root_parent(id, tree_data),
      
      # Calculate depth in tree (0 for roots, 1+ for children)
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

#' Convert MS/MS Spectrum String to Matrix
#'
#' Converts a string representation of an MS/MS spectrum into a numeric matrix.
#' The input string should contain space-separated pairs of m/z values and 
#' intensities joined by colons (e.g., "100.5:1000 200.3:500").
#'
#' @param msms_string Character string containing MS/MS spectrum data in
#'   "mz1:int1 mz2:int2" format. NA values are allowed.
#'
#' @return Two-column numeric matrix with m/z values in first column and 
#'   intensities in second column. Returns NULL for NA input or if no valid
#'   pairs are found.
#'
#' @examples
#' msms <- "100.5:1000 200.3:500"
#' peaks <- convert_msms_to_matrix(msms)
#'
#' @export
convert_msms_to_matrix <- function(msms_string) {
  if(is.na(msms_string)) return(NULL)
  
  pairs <- strsplit(msms_string, " ")[[1]]
  mz_values <- numeric()
  intensity_values <- numeric()
  
  for(pair in pairs) {
    values <- as.numeric(strsplit(pair, ":")[[1]])
    if(length(values) == 2) {
      mz_values <- c(mz_values, values[1])
      intensity_values <- c(intensity_values, values[2])
    }
  }
  
  if(length(mz_values) == 0) return(NULL)
  matrix(c(mz_values, intensity_values), ncol = 2, byrow = FALSE)
}
#' Calculate MS/MS Spectral Similarities and Metadata Relationships
#'
#' This function analyzes metabolite relationships by calculating spectral 
#' similarities and extracting metadata from JSON content. It identifies 
#' parent-child pairs in metabolomic data, computing their MS/MS spectral 
#' similarity and checking for shared characteristics.
#'
#' @param tree_data A data frame containing metabolite tree relationships, 
#'   MS/MS spectra, and JSON metadata with the following columns:
#'   \describe{
#'     \item{id}{Unique metabolite identifier}
#'     \item{duplicate_of}{Parent metabolite ID for this entry}
#'     \item{msms}{MS/MS spectrum as a string representation}
#'     \item{sample}{Sample identifier}
#'     \item{content}{JSON string with detailed sample metadata}
#'   }
#'
#' @return A data frame with the following columns:
#'   \describe{
#'     \item{child_id}{Identifier of the child metabolite}
#'     \item{parent_id}{Identifier of the parent metabolite}
#'     \item{entropy_similarity}{Calculated spectral similarity score between 
#'           parent and child MS/MS spectra}
#'     \item{same_sample}{Logical indicating whether parent and child 
#'           originated from the same sample}
#'     \item{same_organ}{Logical indicating whether parent and child 
#'           are from the same organ}
#'     \item{same_species}{Logical indicating whether parent and child 
#'           are from the same species}
#'   }
#'
#' @details 
#' The function performs the following key operations:
#' \itemize{
#'   \item Parses JSON content to extract organ and species metadata
#'   \item Identifies direct parent-child relationships
#'   \item Calculates entropy-based similarity of MS/MS spectra
#'   \item Checks for shared sample, organ, and species characteristics
#' }
#'
#' @examples
#' \dontrun{
#' library(jsonlite)
#' 
#' # Example tree_data creation
#' tree_data <- data.frame(
#'   id = c("metabolite1", "metabolite2"),
#'   duplicate_of = c(NA, "metabolite1"),
#'   msms = c("100:1000", "100:990"),
#'   sample = c("sample_A", "sample_A"),
#'   content = c(
#'     '{"metadata": {"organ": "liver", "species": "mouse"}}',
#'     '{"metadata": {"organ": "liver", "species": "mouse"}}'
#'   )
#' )
#' 
#' # Calculate similarities
#' similarities <- calculate_direct_similarities(tree_data)
#' }
#'
#' @export
calculate_direct_similarities <- function(tree_data) {
  # First, add robust parsing with error handling
  tree_data <- tree_data %>%
    rowwise() %>%
    mutate(
      organ = tryCatch({
        # Use fromJSON with additional parsing options
        content_parsed <- fromJSON(content, simplifyVector = TRUE, flatten = TRUE)
        
        # Safely extract organ, with fallback to NA
        if(!is.null(content_parsed$metadata$organ)) {
          as.character(content_parsed$metadata$organ)
        } else {
          NA_character_
        }
      }, error = function(e) {
        # Log the error and the problematic content for debugging
        warning(paste("JSON parsing error for content:", content))
        warning(paste("Error details:", e$message))
        NA_character_
      }),
      species = tryCatch({
        content_parsed <- fromJSON(content, simplifyVector = TRUE, flatten = TRUE)
        
        # Safely extract species, with fallback to NA
        if(!is.null(content_parsed$metadata$species)) {
          as.character(content_parsed$metadata$species)
        } else {
          NA_character_
        }
      }, error = function(e) {
        # Log the error and the problematic content for debugging
        warning(paste("JSON parsing error for content:", content))
        warning(paste("Error details:", e$message))
        NA_character_
      })
    ) %>%
    ungroup()
  
  # Rest of the function remains the same...
  
  # Get direct relationships
  direct_pairs <- tree_data %>%
    filter(!is.na(duplicate_of)) %>%
    select(id, duplicate_of, msms, sample, organ, species) %>%
    rename(child_id = id, 
           parent_id = duplicate_of)
  
  # Get parent data
  parent_data <- tree_data %>%
    filter(id %in% direct_pairs$parent_id) %>%
    select(id, msms, sample, organ, species) %>%
    rename(parent_id = id,
           parent_msms = msms,
           parent_sample = sample,
           parent_organ = organ,
           parent_species = species)
  
  # Join parent and child data
  pairs_with_msms <- direct_pairs %>%
    inner_join(parent_data, by = "parent_id") %>%
    filter(!is.na(msms) & !is.na(parent_msms))
  
  # Calculate similarities
  similarity_results <- pairs_with_msms %>%
    rowwise() %>%
    mutate(
      entropy_similarity = {
        parent_peaks <- convert_msms_to_matrix(parent_msms)
        child_peaks <- convert_msms_to_matrix(msms)
        
        if(!is.null(parent_peaks) && !is.null(child_peaks)) {
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
          }, error = function(e) NA_real_)
        } else {
          NA_real_
        }
      },
      same_sample = tolower(trimws(sample)) == tolower(trimws(parent_sample)),
      same_organ = tolower(trimws(organ)) == tolower(trimws(parent_organ)),
      same_species = tolower(trimws(species)) == tolower(trimws(parent_species))
    ) %>%
    ungroup()
  
  return(similarity_results)
}
