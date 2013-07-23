# Function to perform heavy lifting of meta-analysis
meta_heavy <- function(data_files, chron_counter)
{
  for (i in 1:nrow(data_files))
  {
    # Read in a chronology
    data_file <- data_files[i,]
    print(data_file$order)
    rwl <- read_data(data_file)
    
    # Convert it to a tree ring array
    tra <- rwl.to.stra (rwl)
    
    # Descriptive ####
    n_entries <- count_entries(rwl)
    n_series <- ncol(rwl)
   
    # Save results
    path_entries <- "./Meta-analysis/Results/descriptive/entries"
    path_series <- "./Meta-analysis/Results/descriptive/series"
    
    id <- strsplit(as.character(data_file$file), split=".rwl")[[1]][1]
    file_name_entries <- paste(id, "ENTRIES.txt", sep="_")
    file_name_series <- paste(id, "SERIES.txt", sep="_")

    write(n_entries, file=paste(path_entries, file_name_entries, sep="/"))   
    write(n_series, file=paste(path_series, file_name_series, sep="/"))    

    # Sample depth ####
    sample_depth <- sample_depth_tra(tra, 2, sparse=TRUE)
    path <- "./Meta-analysis/Results/descriptive/sample_depth"
    id <- strsplit(as.character(data_file$file), split=".rwl")[[1]][1]
    file_name <- paste(id, "SAMPLE_DEPTH.csv", sep="_")
    write.csv(sample_depth, file=paste(path, file_name, sep="/"))
    
    # Model comparison ####
    # Try all model combinations
    analyses_list <- model_compare(rwl, method="sfs", sparse=TRUE)
    
    # Save the results
    path <- "./Meta-analysis/Results/all_models/analysis"
    id <- strsplit(as.character(data_file$file), split=".rwl")[[1]][1]
    file_name <- paste(id, "ALL_MODELS.RData", sep="_")
    save(analyses_list, file=paste(path, file_name, sep="/"))
    rm(analyses_list)
    
    # Modern sample bias ####
    gam_ita <- standardize_tra(tra, model=list(I=T, T=T, A=T), sparse=TRUE, method="gam")
    gam_ta <- standardize_tra(tra, model=list(I=F, T=T, A=T), sparse=TRUE, method="gam")
    
    analyses_list <- list(gam_ita=gam_ita, gam_ta=gam_ta)
    
    # Save the results
    path <- "./Meta-analysis/Results/gam/analysis"
    id <- strsplit(as.character(data_file$file), split=".rwl")[[1]][1]
    file_name <- paste(id, "GAM.RData", sep="_")
    save(analyses_list, file=paste(path, file_name, sep="/"))
    rm(analyses_list)
    
    # Update your place
    chron_counter <- chron_counter + 1
    write(chron_counter, file="./Meta-analysis/Results/chron_counter.txt")
  }
  return()
}

# Manual code to perform analysis in part ####

# Set size of next analysis
chron_size <- 500

# Initialize counter
# write(1, file="./Meta-analysis/Results/chron_counter.txt")

# Get current position
chron_counter <- scan(file="./Meta-analysis/Results/chron_counter.txt")

# Choose the data files to analyse
sel_data <- data_files[chron_counter:(chron_counter+chron_size-1), ]

# Do the heavy lifting for the next part
meta_heavy(sel_data, chron_counter)