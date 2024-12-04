# Load necessary packages
if (!requireNamespace("httr", quietly = TRUE)) {
  install.packages("httr")
}
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite")
}

library(httr)
library(jsonlite)

# Replace with the ID of your published collection
collection_id <- "7569888"
# Figshare token
token <- "f3bf230bfc60b9907dd2a3ee6f2800714bf45178d9e8406ee4759e3fd392ed17c82ce5a643134afc4f365bad4ff6fd3dcb1e7467a7b36a1358d38c8989474719"

# API endpoint for the articles in the collection
articles_url <- paste0("https://api.figshare.com/v2/collections/", collection_id, "/articles")

# Directory to store the downloaded files
output_dir <- "/Users/dorsi/Desktop/github_tmp"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Fetch articles metadata
response <- httr::GET(articles_url, add_headers(Authorization = paste("token", token)))

if (response$status_code == 200) {
  # Parse the response to extract articles
  articles_data <- jsonlite::fromJSON(httr::content(response, "text"))
  
  # Loop through each article to fetch and download its files
  for (article_id in articles_data$id) {
    # Fetch article metadata
    article_url <- paste0("https://api.figshare.com/v2/articles/", article_id)
    article_response <- httr::GET(article_url, add_headers(Authorization = paste("token", token)))
    
    if (article_response$status_code == 200) {
      article_data <- jsonlite::fromJSON(httr::content(article_response, "text"))
      
      # Extract article title for subdirectory creation
      article_title <- gsub("[^A-Za-z0-9_]", "_", article_data$title)
      article_dir <- file.path(output_dir, article_title)
      if (!dir.exists(article_dir)) {
        dir.create(article_dir)
      }
      
      # Download each file in the article
      for (i in 1:nrow(article_data$files)) {
        file_url <- article_data$files$download_url[[i]]
        file_name <- article_data$files$name[[i]]
        download_path <- file.path(article_dir, file_name)
        
        cat("Downloading:", file_name, "from", file_url, "\n")
        download.file(file_url, destfile = download_path, mode = "wb")
      }
    } else {
      cat("Failed to fetch article metadata for ID:", article_id, "\n")
    }
  }
  cat("All files downloaded successfully!\n")
} else {
  cat("Error fetching articles. Status code:", response$status_code, "\n")
}


