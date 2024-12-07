masswiki_result <- read.csv("data/masswiki_result.csv")

# Create a named vector for splash mapping
splash_mapping <- setNames(masswiki_result$`binbase.splash`, masswiki_result$wiki_id)

# Add splash to each entry in wiki_data
wiki_data <- lapply(names(wiki_data), function(wiki_id) {
  # Copy existing data
  entry <- wiki_data[[wiki_id]]
  
  # Add splash value to the entry
  entry$splash <- splash_mapping[wiki_id]
  
  # Preserve wiki_id in the entry for reference
  entry$wiki_id <- wiki_id
  
  return(entry)
})

# Restore wiki_id as names of the list
names(wiki_data) <- names(wiki_data)


