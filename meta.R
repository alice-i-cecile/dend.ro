# Libraries ####
library(ggplot2)
library(reshape2)
library(dplR)

# Set RNG seed ####
set.seed(seed=42)

# Generate a list of file names ####

# Find the root folder with the data
data_folder <- "./Meta-analysis/Extracted"

# Find the folders for each region
region_folders <- list.dirs(data_folder)

# Make a list of each file, indexed by which region it's from
data_files <- melt(sapply(region_folders, list.files))
colnames(data_files) <- c("file", "folder")

# Reading and truncating data ####
read_data <- function(data_file)
{
  data_path <- paste(data_file$folder, data_file$file, sep="/")
  rwl <- read.rwl()
  return(rwl)
}

# Descripitive statistics ####

# Find the number of measurements
count_entries <- function(rwl)
{
  sum(!is.na(rwl))
}

# Find the number of series (raw)
count_series <- function(rwl)
{
  ncol(rwl)
}


# Find the number of chronologies

# Find the sample depth by time

# Plot the graph of sample depth vs. time