## Use SQL to query data

# Install necessary packages
install.packages("DBI")
install.packages("RPostgres")

# Load the packages
library(DBI)
library(RPostgres)


# Establish a connection
con <- dbConnect(
  RPostgres::Postgres(),
  dbname = "carrot-prod",
  host = "lcb-standalone-cluster.cluster-czbqhgrlaqbf.us-west-2.rds.amazonaws.com",
  port = 5432,  # Default PostgreSQL port
  user = "zyang2k",
  password = "EzTwsGPzivZvpcN"
)

# change LIMIT if necessary
query <- "
SELECT *
FROM compound
WHERE method = '5m hilic premier | orbitrap | beh amide | negative';
"

# Execute the query and store the result in a data frame
data <- dbGetQuery(con, query)

