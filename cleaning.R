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