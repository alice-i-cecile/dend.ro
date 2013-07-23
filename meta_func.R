# Libraries ####
library(ggplot2)
library(reshape2)
library(dplR)
library(plyr)
library(pastecs)

source("standardization.R")
source("coreTRA.R")

# Utility functions ####
meta_apply <- function(data_files, FUN, ...)
{
  out <- list()
  for (i in 1:nrow(data_files))
  {
    print(toString(data_files[i, "file"]))
    rwl <- read_data(data_files[i,])
    out_i <- FUN(rwl, ...)
    out[[i]] <- out_i
  }
  
  return(out)
}

# Reading data
read_data <- function(data_file)
{
  data_path <- paste(data_file$folder, data_file$file, sep="/")
  rwl <- read.rwl(data_path)
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

# Find the sample depth by time
sample_depth_rwl <- function(rwl, factor.dim=2, sparse=TRUE)
{
  if(sparse)
  {
    tra <- rwl.to.stra(rwl)
    
  } else
  {
    tra <- rwl.to.tra(rwl)
  }
  
  return(sample_depth_tra(tra, factor.dim, sparse))
}

# Comparing models ####
# Perform all 16 analyses on a single data set
# Exclude cases where form and error do not match
all_analyses <- function(tra, method="sfs", sparse=TRUE, ...)
{
  # Convert the tree-ring array to the appropriate form (sparse/full)
  if (sparse) 
  {
    if (!is.data.frame(tra))
    {
      tra <- sparse_tra(tra)
    }
  } else
  {
    if (is.data.frame(tra))
    {
      tra <- unsparse_tra(tra)
    }
  }
  
  # Building the set of choices for analysis
  forms <- c("additive", "multiplicative")
  errors <- c("norm", "lnorm")
  effect_I <- c(TRUE, FALSE)
  effect_T <- c(TRUE, FALSE)
  effect_A <- c(TRUE, FALSE)
  
  params <- expand.grid (forms, errors, effect_I, effect_T, effect_A)
  names(params) <- c("form", "error", "model.I", "model.T", "model.A")
  
  # Exclude models where form and error do not match
  add_norm <-  params$form == "additive" & params$error == "norm"
  mult_lnorm <-  params$form == "multiplicative" & params$error == "lnorm"
  params <- params [add_norm | mult_lnorm, ]
    
  # Wrapper for extracting parameters
  std_wrapper <- function(param, ...)
  {
    form <- param$form
    error <- param$error
    model <- list(I=param$model.I, T=param$model.T, A=param$model.A)
    
    out <- standardize_tra(tra, model, form, error, method, sparse, ...)
    return (out)
  }
  
  # Run each analysis on the data
  analyses_list <- lapply(1:nrow(params), 
                          function(x, ...)
                          {
                            std_wrapper(params[x,], ...)
                          },
                          ...
                          )
  
  return(analyses_list)
  
}

# Clean and perform all analyses
model_compare <- function (rwl, method="sfs", sparse=TRUE, ...)
{
  # Clean 0 and negative values out of the data
  rwl[rwl<=0] <- NA
  
  # Convert to tra
  if (sparse)
  {
    tra <- rwl.to.stra(rwl)
  } else {
    tra <- rwl.to.tra(rwl)
  }
  
  # Perform analyses
  analyses_list <- all_analyses(tra, method, sparse, ...)
  
  return(analyses_list)
}

# Extracting saved info
load_analysis <- function(data_file, path=".", extension=".RData")
{
  id <- strsplit(as.character(data_file$file), split=".rwl")[[1]][1]
  file_name <- paste(id, extension, sep="_")
  load(paste(path, file_name, sep="/"))
  out <- analyses_list
  return(out)
}

extract_fit_and_settings <- function(analyses)
{
  fit_df <- sapply(analyses, function(x) x$fit[3:13])
  settings_df <- sapply(analyses, function(x)
    {
     model_settings <- unlist(x$settings$model)
     x$settings$model <- NULL
     x$settings <- lapply(x$settings, as.character)
     x$settings <- c(x$settings[c("form", "error")], model_settings)
     out <- x$settings
     return(out)
    }
  )
  run_df <- t(rbind(fit_df, settings_df))
  
  return(run_df)
}

store_fit_and_settings <- function(data_file, path=".", extension=".RData")
{
  analyses <- load_analysis(data_file, path, extension)
  run_df <- extract_fit_and_settings(analyses)
  stored_df <- cbind(data_file[1:3], run_df)
  return(stored_df)
}

# Coding for graphs
generate_code <- function(model.I, model.T, model.A)
{
  if (sum(model.I, model.T, model.A)==0)
  {
    return ("Null")
  }
  
  code <- "" 
  if (model.I)
  {
   code <- paste(code, "I", sep=",")
  }
  if (model.T)
  {
    code <- paste(code, "T", sep=",")
  }
  if (model.A)
  {
    code <- paste(code, "A", sep=",")
  }
  
  code <- substr(code, 2, nchar(code))
  
  return(code)
}