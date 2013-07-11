rm(list=ls())
start_time <- Sys.time()

# Libraries ####
library(ggplot2)
library(reshape2)
library(dplR)
library(plyr)
library(pastecs)

source("standardization.R")
source("coreTRA.R")

# Saving parameters ####

# Set the height and width of graphics (cm)
save_width <- 17.2
save_height <- 17.2

# Set RNG seed and number of chronoologies to analyse ####
# Ensure consistent RNG
set.seed(seed=42)
N <- 10

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

# Truncate the data set to a resonable size ####
full_data_files <- read.csv("./Meta-analysis/Results/data_files.csv")

# Only select valid data sets
selected_data <- sample((1:nrow(full_data_files))[full_data_files$valid], size=N, replace=F)
data_files <- full_data_files[selected_data, ]

# Descripitive statistics ####
descriptive_stats <- data_files[1:2]

# Find the number of measurements
count_entries <- function(rwl)
{
  sum(!is.na(rwl))
}

n_entries <- unlist(meta_apply(data_files, count_entries))
descriptive_stats$entries <- n_entries 
print(paste(sum(n_entries), "entries analysed."))

# Find the number of series (raw)
count_series <- function(rwl)
{
  ncol(rwl)
}

n_series <- unlist(meta_apply(data_files, count_series))
descriptive_stats$series <- n_series
print(paste(sum(n_series), "series analysed."))

# Find the number of chronologies
n_chronologies <- nrow(data_files) 
print(paste(n_chronologies, "chronologies analysed."))

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

sample_depth_list <- meta_apply(data_files, sample_depth_rwl)

# Bind the results together into a matrix
first_years <- sapply(sample_depth_list, function(x){min(as.numeric(names(x)))})
last_years <- sapply(sample_depth_list, function(x){max(as.numeric(names(x)))})

first_year <- min(first_years)
last_year <- max(last_years)

sample_depth_matrix <- matrix(data=0, nrow=N, ncol=(last_year-first_year+1))
colnames(sample_depth_matrix) <- as.character(first_year:last_year)

for (i in 1:N)
{
  cols_i <- first_years[i] : last_years[i] - first_year + 1
  sample_depth_matrix[i, cols_i] <- sample_depth_list[[i]]
}

sample_depth_df <- cbind(data_files[1:2], sample_depth_matrix)

# Plot the graph of sample depth vs. time
summary_sample_depth <- data.frame(depth=colSums(sample_depth_matrix), Year=as.numeric(colnames(sample_depth_matrix)))

sample_depth_plot <- ggplot(summary_sample_depth, aes(y=depth, x=Year)) + geom_area() + theme_bw() + ylab("Number of series analyzed")

# Save results
write.csv(descriptive_stats, file="./Meta-analysis/Results/Descriptive/descriptive.csv")
write.csv(sample_depth_df, file="./Meta-analysis/Results/Descriptive/sample_depth.csv")
ggsave(filename="./Meta-analysis/Results/Descriptive/sample_depth.svg", plot=sample_depth_plot, width=save_width, height=save_height, units="cm")
ggsave(filename="./Meta-analysis/Results/Descriptive/sample_depth.pdf", plot=sample_depth_plot, width=save_width, height=save_height, units="cm")

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

# Loop for meta-analysis

# Analyses take a very long time. Save the results as you go
# Keep track of where you've stopped

# Initialize counter
write(1, file="./Meta-analysis/Results/all_models/chron_counter.txt")

# Pick up where you last left off
chron_counter <- scan(file="./Meta-analysis/Results/all_models/chron_counter.txt")

for (i in chron_counter:nrow(data_files))
{
  
  # Load and analyse the data
  rwl <- read_data(data_files[i,])
  analyses_list <- model_compare(rwl, method="sfs", sparse=TRUE)
  
  # Save the results
  path <- "./Meta-analysis/Results/all_models/analysis"
  id <- strsplit(as.character(data_files[i, "file"]), split=".rwl")[[1]][1]
  file_name <- paste(id, "ALL_MODELS.RData", sep="_")
  save(analyses_list, file=paste(path, file_name, sep="/"))
  rm(analyses_list)
  
  # Update your place
  chron_counter <- chron_counter + 1
  write(chron_counter, file="./Meta-analysis/Results/all_models/chron_counter.txt")
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

model_comp <- do.call(rbind, 
  lapply(1:nrow(data_files), 
    function(x, path, extension){
      store_fit_and_settings(data_files[x,], path, extension)
    }
    , path="./Meta-analysis/Results/all_models/analysis", extension="ALL_MODELS.RData"
  )
)

# Processing model comparison results
model_comp <- as.data.frame(model_comp)
model_comp[4:14] <- sapply (model_comp[4:14], as.numeric)
model_comp[15:16] <- sapply(model_comp[15:16], unlist)
model_comp[17:19] <- sapply(model_comp[17:19], unlist)

# Add delta *IC information
model_comp <- ddply(model_comp, "file", transform, 
                    dAIC=AIC-min(AIC),
                    dAICc=AICc-min(AICc),
                    dBIC=BIC-min(BIC)
                    )

# Add *IC model weights
model_comp <- ddply(model_comp, "file", transform, 
                    wAIC=exp(-0.5*dAIC) / sum(exp(-0.5*dAIC)),
                    wAICc=exp(-0.5*dAICc) / sum(exp(-0.5*dAICc)),
                    wBIC=exp(-0.5*dBIC) / sum(exp(-0.5*dBIC))
                    )

# Show which model is preferred by each metric
model_comp <- ddply(model_comp, "file", transform, 
                    sigma_choice=which.min(sigma),
                    Rsq_choice=which.max(Rsq),
                    adj.Rsq_choice=which.max(adj.Rsq),
                    llh_choice=which.max(llh),
                    AIC_choice=which.min(AIC),
                    AICc_choice=which.min(AICc),
                    BIC_choice=which.min(BIC)
)

# Compute best case values for each relevant metric
model_comp <-  ddply(model_comp, "file", transform, 
                     sigma_best=min(sigma),
                     Rsq_best=max(Rsq),
                     adj.Rsq_best=max(adj.Rsq),
                     wAIC_best=max(wAIC),
                     wAICc_best=max(wAICc),
                     wBIC_best=max(wBIC)
)

# Compute mean of best case values for each metric
mean_best_fit <- c(
                     sigma=mean(model_comp$sigma_best),
                     Rsq=mean(model_comp$Rsq_best),
                     adj.Rsq=mean(model_comp$adj.Rsq_best),
                     wAIC=mean(model_comp$wAIC_best),
                     wAICc=mean(model_comp$wAICc_best),
                     wBIC=mean(model_comp$wBIC_best)
                     )

# Compute sd of best case values for each metric
sd_best_fit <- c(
  sigma=sd(model_comp$sigma_best),
  Rsq=sd(model_comp$Rsq_best),
  adj.Rsq=sd(model_comp$adj.Rsq_best),
  wAIC=sd(model_comp$wAIC_best),
  wAICc=sd(model_comp$wAICc_best),
  wBIC=sd(model_comp$wBIC_best)
)

# Summarize typical values for best fitting models
best_fit_df <- data.frame(mean=mean_best_fit, sd=sd_best_fit)

# Summarize model fit stats for each class of of models
mean_class_fit <-   ddply(model_comp, c("form","error", "I", "T", "A"), summarise, 
                          sigma=mean(sigma),
                          Rsq=mean(Rsq),
                          adj.Rsq=mean(adj.Rsq),
                          dAIC=mean(dAIC),
                          dAICc=mean(dAICc),
                          dBIC=mean(dBIC),
                          wAIC=mean(wAIC),
                          wAICc=mean(wAICc),
                          wBIC=mean(wBIC)
                        )
  
  
sd_class_fit <-   ddply(model_comp, c("form","error", "I", "T", "A"), summarise, 
                          sigma=sd(sigma),
                          Rsq=sd(Rsq),
                          adj.Rsq=sd(adj.Rsq),
                          dAIC=sd(dAIC),
                          dAICc=sd(dAICc),
                          dBIC=sd(dBIC),
                          wAIC=sd(wAIC),
                          wAICc=sd(wAICc),
                          wBIC=sd(wBIC)
)

# Determine how many times each model was selected by each criteria
selection_occurences <- data.frame( 
                                    model=1:16,                                
                                    model_comp[1:16, 15:19],
                                    sigma_choice=NA,
                                    Rsq_choice=NA,
                                    adj.Rsq_choice=NA,
                                    llh_choice=NA,
                                    AIC_choice=NA,
                                    AICc_choice=NA,
                                    BIC_choice=NA
                                  )

for (i in 1:16)
{
  for(j in 7:13)
  {
    criteria <- names(selection_occurences)[j]
    selection_occurences[i, j] <- sum(model_comp[criteria]==i)
  }
}

selection_frequency <- selection_occurences
selection_frequency[7:13] <- selection_occurences[7:13]/colSums(selection_occurences[7:13])

# Graph the selection frequency
sel_freq_melt <- melt(selection_frequency, id.vars=c("model", "form", "error", "I", "T", "A"))
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
  
sel_freq_melt$code <-(mapply(generate_code, model.I=sel_freq_melt$I, model.T=sel_freq_melt$T, model.A=sel_freq_melt$A))

code_order <- rev(c("I,T,A", "I,T", "I,A", "T,A", "I", "T", "A", "Null"))


sel_freq_plot <- ggplot(sel_freq_melt, aes(x=variable, y=code, fill=value)) + geom_tile() + scale_fill_gradient(name="Selection frequency", low="white", high="black") + facet_grid(~form) + xlab("Selection criteria") + ylab("Model") + ylim(levels=code_order)

# Compare number of parameters in selected model to number of data points
selected_model <- ddply(model_comp, ~file, summarize, 
                        n=n,
                        sigma=sigma_choice,
                        Rsq=Rsq_choice,
                        adj.Rsq=adj.Rsq_choice,
                        llh=llh_choice,
                        AIC=AIC_choice,
                        AICc=AICc_choice,
                        BIC=BIC_choice)

num_param <- c(rep(3,2), rep(2,4), rep(1,2), rep(2,2), rep(1,4), rep(0,2))
selected_model[3:9] <- sapply(selected_model[3:9], function(x){num_param[x]})
param_melt <- melt(selected_model, id.vars=c("n", "file"))


n_vs_effects_plot <- ggplot(param_melt, aes(x=n, y=value)) + geom_point() + facet_wrap(~variable, ncol=2) + theme_bw() +ylab("Number of effects in selected model") + xlab("Number of data points")

# Save results
write.csv(model_comp, file="./Meta-analysis/Results/all_models/model_comparison.csv")

write.csv(best_fit_df, file="./Meta-analysis/Results/all_models/best_fit.csv")

write.csv(mean_class_fit, file="./Meta-analysis/Results/all_models/mean_class_fit.csv")
write.csv(sd_class_fit, file="./Meta-analysis/Results/all_models/mean_class_fit.csv")

write.csv(selection_occurences, file="./Meta-analysis/Results/all_models/selection_occurences.csv")
write.csv(selection_frequency, file="./Meta-analysis/Results/all_models/selection_frequency.csv")
ggsave(filename="./Meta-analysis/Results/all_models/selection_frequency.svg", plot=sel_freq_plot, width=save_width, height=save_height, units="cm")
ggsave(filename="./Meta-analysis/Results/all_models/selection_frequency.pdf", plot=sel_freq_plot, width=save_width, height=save_height, units="cm")

ggsave(filename="./Meta-analysis/Results/all_models/n_vs_effects.svg", plot=n_vs_effects_plot, width=save_width, height=save_height, units="cm")
ggsave(filename="./Meta-analysis/Results/all_models/n_vs_effects.pdf", plot=n_vs_effects_plot, width=save_width, height=save_height, units="cm")

# Testing for modern sample bias ####

# Perform GAM-(I)TA multiplicative, log-normal standardization on each data set

# Loop for meta-analysis

# Analyses take a very long time. Save the results as you go
# Keep track of where you've stopped

# Initialize counter
write(1, file="./Meta-analysis/Results/gam/chron_counter.txt")

# Pick up where you last left off
chron_counter <- scan(file="./Meta-analysis/Results/gam/chron_counter.txt")

for (i in chron_counter:nrow(data_files))
{
  
  # Load and analyse the data
  rwl <- read_data(data_files[i,])
  stra <- rwl.to.stra (rwl)
  gam_ita <- standardize_tra(stra, model=list(I=T, T=T, A=T), sparse=TRUE, method="gam")
  gam_ta <- standardize_tra(stra, model=list(I=F, T=T, A=T), sparse=TRUE, method="gam")
  
  analyses_list <- list(gam_ita=gam_ita, gam_ta=gam_ta)
  
  # Save the results
  path <- "./Meta-analysis/Results/gam/analysis"
  id <- strsplit(as.character(data_files[i, "file"]), split=".rwl")[[1]][1]
  file_name <- paste(id, "GAM.RData", sep="_")
  save(analyses_list, file=paste(path, file_name, sep="/"))
  rm(analyses_list)
  
  # Update your place
  chron_counter <- chron_counter + 1
  write(chron_counter, file="./Meta-analysis/Results/gam/chron_counter.txt")
}

# Load analyses and extract model fit statistics
model_comp <- do.call(rbind, 
                      lapply(1:nrow(data_files), 
                             function(x){
                               store_fit_and_settings(data_files[x,], "./Meta-analysis/Results/gam/analysis", "GAM.RData")
                             }
                      )
)

# Processing model comparison results
model_comp <- as.data.frame(model_comp)
model_comp[4:14] <- sapply (model_comp[4:14], as.numeric)
model_comp[15:16] <- sapply(model_comp[15:16], unlist)
model_comp[17:19] <- sapply(model_comp[17:19], unlist)

# Add delta *IC information
model_comp <- ddply(model_comp, "file", transform, 
                    dAIC=AIC-min(AIC),
                    dAICc=AICc-min(AICc),
                    dBIC=BIC-min(BIC)
)

# Add *IC model weights
model_comp <- ddply(model_comp, "file", transform, 
                    wAIC=exp(-0.5*dAIC) / sum(exp(-0.5*dAIC)),
                    wAICc=exp(-0.5*dAICc) / sum(exp(-0.5*dAICc)),
                    wBIC=exp(-0.5*dBIC) / sum(exp(-0.5*dBIC))
)

# Show which model is preferred by each metric
model_comp <- ddply(model_comp, "file", transform, 
                    sigma_choice=which.min(sigma),
                    Rsq_choice=which.max(Rsq),
                    adj.Rsq_choice=which.max(adj.Rsq),
                    llh_choice=which.max(llh),
                    AIC_choice=which.min(AIC),
                    AICc_choice=which.min(AICc),
                    BIC_choice=which.min(BIC)
)

# Report mean and sd of fit statistics for each model
mean_class_fit <-   ddply(model_comp, c("form","error", "I", "T", "A"), summarise, 
                          sigma=mean(sigma),
                          Rsq=mean(Rsq),
                          adj.Rsq=mean(adj.Rsq),
                          dAIC=mean(dAIC),
                          dAICc=mean(dAICc),
                          dBIC=mean(dBIC),
                          wAIC=mean(wAIC),
                          wAICc=mean(wAICc),
                          wBIC=mean(wBIC)
)

sd_class_fit <-   ddply(model_comp, c("form","error", "I", "T", "A"), summarise, 
                        sigma=sd(sigma),
                        Rsq=sd(Rsq),
                        adj.Rsq=sd(adj.Rsq),
                        dAIC=sd(dAIC),
                        dAICc=sd(dAICc),
                        dBIC=sd(dBIC),
                        wAIC=sd(wAIC),
                        wAICc=sd(wAICc),
                        wBIC=sd(wBIC)
)

# Determine how many times each model was selected by each criteria
selection_occurences <- data.frame( 
  model=1:2,                                
  model_comp[1:2, 15:19],
  sigma_choice=NA,
  Rsq_choice=NA,
  adj.Rsq_choice=NA,
  llh_choice=NA,
  AIC_choice=NA,
  AICc_choice=NA,
  BIC_choice=NA
)

for (i in 1:2)
{
  for(j in 7:13)
  {
    criteria <- names(selection_occurences)[j]
    selection_occurences[i, j] <- sum(model_comp[criteria]==i)
  }
}

selection_frequency <- selection_occurences
selection_frequency[7:13] <- selection_occurences[7:13]/colSums(selection_occurences[7:13])

# Graph the selection frequency
sel_freq_melt <- melt(selection_frequency, id.vars=c("model", "form", "error", "I", "T", "A"))
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

sel_freq_melt <- sel_freq_melt[!is.na(sel_freq_melt$form),]
sel_freq_melt$code <- mapply(generate_code, model.I=sel_freq_melt$I, model.T=sel_freq_melt$T, model.A=sel_freq_melt$A)

code_order <- rev(c("I,T,A", "T,A"))

sel_freq_plot <- ggplot(sel_freq_melt, aes(x=variable, y=code, fill=value)) + geom_tile() + scale_fill_gradient(name="Selection frequency", low="white", high="black") + ylab("Model") + ylim(levels=code_order)

# Extract the effect vectors of the gam models
extract_effects <- function(data_file,  path=".", extension=".RData")
{
  analyses <- load_analysis(data_file, path, extension)
  effects_list <- lapply(analyses, function(x) x$effects) 
  return(effects_list)
}

gam_effects <- lapply(1:nrow(data_files), function(x, path, extension)
  {
    extract_effects(data_files[x,], path, extension)
  }
  ,path="./Meta-analysis/Results/gam/analysis", extension="GAM.RData"
)

# Calculate the ratio (TA / ITA) of estimated time effects at each year
msb_ratios <- lapply(gam_effects, function(x)
  {
    T_ita <- x$gam_ita$T
    T_ta <- x$gam_ta$T
    return(T_ta/T_ita)
  }
)

# Check for trends in these ratios
msb_trends <- melt(sapply(msb_ratios, function(x) cor.test(x, 1:length(x), method="spearman")$estimate))

msb_trend_plot <- ggplot(msb_trends, aes(x=value)) + geom_density(fill="black") + theme_bw() + xlab("Spearman's rho test for trend in ratios of G=TA to G=ITA time effects") + ylab("Density")

# Graph the ratio of TA / ITA time effects by year
msb_points <- unlist(msb_ratios)
msb_df <- data.frame(ratio=msb_points, year=as.numeric(names(msb_points)))

msb_plot <- ggplot(msb_df, aes(x=year, y=ratio)) + geom_point(alpha=0.2) + geom_smooth(colour="black") + theme_bw() + xlab("Year") + ylab("Ratio of G=TA to G=ITA time effects") + geom_hline(y=1)

# Save results
write.csv(model_comp, file="./Meta-analysis/Results/gam/model_comparison.csv")

write.csv(mean_class_fit, file="./Meta-analysis/Results/gam/mean_class_fit.csv")
write.csv(sd_class_fit, file="./Meta-analysis/Results/gam/mean_class_fit.csv")

write.csv(selection_occurences, file="./Meta-analysis/Results/gam/selection_occurences.csv")
write.csv(selection_frequency, file="./Meta-analysis/Results/gam/selection_frequency.csv")
ggsave(filename="./Meta-analysis/Results/gam/selection_frequency.svg", plot=sel_freq_plot, width=save_width, height=save_height, units="cm")
ggsave(filename="./Meta-analysis/Results/gam/selection_frequency.pdf", plot=sel_freq_plot, width=save_width, height=save_height, units="cm")

save(msb_ratios, file="./Meta-analysis/Results/gam/msb_ratios.RData")
ggsave(filename="./Meta-analysis/Results/gam/msb_trends.svg", plot=msb_trend_plot, width=save_width, height=save_height, units="cm")
ggsave(filename="./Meta-analysis/Results/gam/msb_trends.pdf", plot=msb_trend_plot, width=save_width, height=save_height, units="cm")
ggsave(filename="./Meta-analysis/Results/gam/msb.svg", plot=msb_plot, width=save_width, height=save_height, units="cm")
ggsave(filename="./Meta-analysis/Results/gam/msb.pdf", plot=msb_plot, width=save_width, height=save_height, units="cm")

# Print time elapsed
finish_time <- Sys.time()
print(start_time)
print(finish_time)