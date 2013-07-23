# Generate a list of file names ####

# Find the root folder with the data
data_folder <- "./Meta-analysis/Extracted"

# Find the folders for each region
region_folders <- list.dirs(data_folder)

# Make a list of each file, indexed by which region it's from
full_data_files <- melt(sapply(region_folders, list.files))
colnames(full_data_files) <- c("file", "folder")

# Reading and selecting data ####
read_data <- function(data_file)
{
  data_path <- paste(data_file$folder, data_file$file, sep="/")
  rwl <- read.rwl(data_path)
  return(rwl)
}

# Ensure that each file loads correctly ####
valid <- sapply(1:nrow(full_data_files), function(x) !inherits(try(read_data(full_data_files[x,])), "try-error"))

# Flag broken data sets
full_data_files$valid <- valid

# Save list of data sets with flags
write.csv(full_data_files, file="./Meta-analysis/Results/data_files.csv")

# Clean and pick an order ####
full_data_files <- read.csv("./Meta-analysis/Results/data_files.csv")

# Set RNG seed and number of chronoologies to analyse
# Ensure consistent RNG
set.seed(seed=42)

# Only select valid data sets
data_files <- full_data_files[full_data_files$valid==TRUE, ]
# Find a predetermined random order for the files
data_files <- data_files[sample(1:nrow(data_files)), ]
data_files$order <- 1:nrow(data_files)
write.csv(data_files, file="./Meta-analysis/Results/ordered_data_files.csv")

# Load from file
data_files <- read.csv(file="./Meta-analysis/Results/ordered_data_files.csv")

# Some data sets are too large for GAM (mgcv) to handle!
# Record and skip them
data_files$size_ok <- TRUE
data_files$size_ok[47] <- FALSE
data_files$size_ok[2085] <- FALSE
data_files$valid[2185] <- FALSE

write.csv(data_files, file="./Meta-analysis/Results/ordered_data_files.csv")
