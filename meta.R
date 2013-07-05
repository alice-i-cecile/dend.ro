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

truncate_rwl <- function(rwl, max_entries=10000)
{
  selected_entries <- 0
  selected_series <- rwl[, NULL]
  random_order <- sample(names(rwl))
  
  for (i in 1:ncol(rwl))
  {
    if (selected_entries > max_entries)
    {
      break
    }
    
    series_i <- rwl[random_order[i]]
    num_entries <- sum(!is.na(series_i])
    selected_entries <- selected_entries + num_entries
    selected_series <- cbind(selected_series, series_i)
  }

  return(series_i)

}

# Truncating chronologies ####

# Load the chronology

# Truncated the chronology

# Save the truncated chronologies

# Construct a list of the truncated chronology files

# Descripitive statistics ####

# Find the number of measurements (raw)
count_entries <- function(rwl)
{
  sum(!is.na(rwl))
}

# Find the number of measurements (truncated)

# Find the number of series (raw)
count_series <- function(rwl)
{
  ncol(rwl)
}


# Find the number of series (truncated)

# Find the number of chronologies

# Find the sample depth by time (raw)

# Find the sample depth by time (truncated)

# Plot the graph of sample depth vs. time