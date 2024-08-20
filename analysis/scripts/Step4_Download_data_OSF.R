
## STEP 4. DOWNLOAD RAW DATA FROM OSF

# Load required packages
library(utils)
library(osfr)

###------------------------------------------------------------------------
## Core Data
# Create the directory to store core data if it doesn't exist
core_data <- file.path("analysis", "data", "raw_data", "core_data")
dir.create(core_data, recursive = TRUE, showWarnings = FALSE)

# Download core data from OSF (alternatively, use the link in a web browser to download the .zip file)
osf_retrieve_file("https://osf.io/mt4xv") %>%
  osf_download(path = "analysis/data/raw_data/core_data")

# Specify the path to the Zip file
zip_file <- "analysis/data/raw_data/core_data/core_data.zip"

# Extract files from the Zip archive
unzip(zip_file, exdir = "analysis/data/raw_data/core_data")



###------------------------------------------------------------------------
## Product Data
# Create the directory to store core data if it doesn't exist
product_data <- file.path("analysis", "data", "raw_data", "product_data")
dir.create(product_data, recursive = TRUE, showWarnings = FALSE)

# Download core data from OSF (alternatively, use the link in a web browser to download the .zip file)
osf_retrieve_file("") %>%
  osf_download(path = "analysis/data/raw_data/product_data")

# Specify the path to the Zip file
zip_file <- "analysis/data/raw_data/product_data.zip"

# Extract files from the Zip archive
unzip(zip_file, exdir = "analysis/data/raw_data/product_data")
