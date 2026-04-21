# sql_retrieve_data.R
#
# Description: This script connects to the LCBinBase PostgreSQL database and 
#              retrieves compound data for specific HILIC methods.
#
# Author: Original by zyang2k
# Last Modified: 2025-01-14

#------------------------------------------------------------------------------
# Required Packages
#------------------------------------------------------------------------------
library(DBI)       # For database connection handling
library(RPostgres) # For PostgreSQL specific functionality

#------------------------------------------------------------------------------
# Function Definitions
#------------------------------------------------------------------------------

#' Create Database Connection
#' 
#' Establishes a connection to the LCBinBase PostgreSQL database with the specified
#' credentials and connection parameters.
#'
#' @param dbname (character): Name of the database
#' @param host (character): Database host address
#' @param port (numeric): Port number for the connection
#' @param user (character): Username for authentication
#' @param password (character): Password for authentication
#'
#' @return DBIConnection: A database connection object
#'
#' @examples
#' con <- create_db_connection(
#'   dbname = "carrot-prod",
#'   host = "lcb-standalone-cluster.example.com",
#'   port = 5432,
#'   user = "username",
#'   password = "password"
#' )
create_db_connection <- function(dbname, host, port = 5432, user, password) {
  tryCatch({
    con <- dbConnect(
      RPostgres::Postgres(),
      dbname = dbname,
      host = host,
      port = port,
      user = user,
      password = password
    )
    return(con)
  }, error = function(e) {
    stop("Failed to connect to database: ", e$message)
  })
}

#' Query Database Table by Method
#' 
#' Retrieves data from a specified database table filtered by method.
#'
#' @param con (DBIConnection): Active database connection
#' @param method (character): The specific method to query for
#' @param table (character): Name of the database table to query
#' @param save_path (character, optional): If provided, saves results to this file path
#'
#' @return data.frame: Results of the database query
#'
#' @examples
#' # Basic query
#' compound_data <- query_by_method(
#'   con = my_connection,
#'   method = "5m hilic premier | orbitrap | beh amide | negative",
#'   table = "compound"
#' )
#'
#' # Query with automatic saving
#' sample_annotations <- query_by_method(
#'   con = my_connection,
#'   method = "5m hilic premier | orbitrap | beh amide | negative",
#'   table = "sample_annotations",
#'   save_path = "results/annotations.RData"
#' )
query_by_method <- function(con, method, table, save_path = NULL) {
  tryCatch({
    # Construct and execute query
    query <- sprintf("
      SELECT *
      FROM %s
      WHERE method = '%s';
    ", table, method)
    
    result <- dbGetQuery(con, query)
    
    # Optionally save results
    if (!is.null(save_path)) {
      tryCatch({
        save(result, file = save_path)
        message(sprintf("Results saved to: %s", save_path))
      }, error = function(e) {
        warning("Failed to save results: ", e$message)
      })
    }
    
    return(result)
  }, error = function(e) {
    stop(sprintf("Failed to execute query on table '%s': %s", table, e$message))
  })
}

# Example usage:
# Create database connection
con <- create_db_connection(
  dbname = "carrot-prod",
  host = "lcb-standalone-cluster.cluster-czbqhgrlaqbf.us-west-2.rds.amazonaws.com",
  port = 5432,
  user = "zyang2k",
  password = "EzTwsGPzivZvpcN"
)

# Query different tables
# compounds <- query_by_method(
#   con = con,
#   method = "5m hilic premier | orbitrap | beh amide | negative",
#   table = "compound"
# )
# #
# annotations <- query_by_method(
#   con = con,
#   method = "5m hilic premier | orbitrap | beh amide | negative",
#   table = "sample_annotations"
# )

binbase_hilic_neg <- dbGetQuery(con, "SELECT * FROM compound WHERE method = '5m hilic premier | orbitrap | beh amide | negative' AND hidden = FALSE;")


# Close the database connection
dbDisconnect(con)

#------------------------------------------------------------------------------
# Usage Notes:
# 1. Ensure all required packages are installed:
#    install.packages(c("DBI", "RPostgres"))
# 2. Store database credentials securely (e.g., in environment variables)
# 3. Consider implementing connection pooling for multiple queries
# 4. Add appropriate error handling for your use case
# 5. For large datasets, consider using dbSendQuery() instead of dbGetQuery()
#------------------------------------------------------------------------------