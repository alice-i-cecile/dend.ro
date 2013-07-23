source("meta_func.R")

# Saving parameters ####

# Set the height and width of graphics (cm)
save_width <- 20
save_height <- 20

# Find which chronologies have been analysed ####

# Explore the directory to find the list of files
search_path <- "./Meta-analysis/Results/gam/analysis"
files_found <- list.files(search_path)

# Strip the suffix
completed <- sapply(files_found, function(x) strsplit(x, split="_")[[1]][1])

# Cross-reference against the list of files
clean_names <- sapply(data_files$file, function(x)
{
  x<-as.character(x) 
  out <- substr(x, 1, nchar(x)-4)
  return(out)
}
)

data_files$completed <- clean_names %in% completed
analysed <- data_files[data_files$completed, ]

# Record the list of analysed chronologies
write.csv(data_files, file="./Meta-analysis/Results/ordered_data_files.csv")

# Descriptive stats ####
descriptive_stats <- data.frame(id=clean_names[clean_names %in% completed])

# Number of entries
descriptive_stats$entries <- sapply(descriptive_stats$id, function(x)
  {
    suffix <- "ENTRIES.txt"
    dir <- "./Meta-analysis/Results/descriptive/entries"
    file <- paste(x, suffix, sep="_")
    path <- paste(dir, file, sep="/")
    
    entries <- scan(path)
    return(entries)
  }
)

print(paste(sum(descriptive_stats$entries), "entries analysed."))

# Number of series
descriptive_stats$series <- sapply(descriptive_stats$id, function(x)
  {
    suffix <- "SERIES.txt"
    dir <- "./Meta-analysis/Results/descriptive/series"
    file <- paste(x, suffix, sep="_")
    path <- paste(dir, file, sep="/")
    
    series <- scan(path)
    return(series)
  }
)

print(paste(sum(descriptive_stats$series), "series analysed."))

# Number of chronologies
n_chronologies <- nlevels(descriptive_stats$id)

print(paste(n_chronologies, "chronologies analysed."))

# Sample depth ####
sample_depth_list <- lapply(descriptive_stats$id, function(x)
  {
    suffix <- "SAMPLE_DEPTH.csv"
    dir <- "./Meta-analysis/Results/descriptive/sample_depth"
    file <- paste(x, suffix, sep="_")
    path <- paste(dir, file, sep="/")
    
    sample_depth <- read.csv(path)
    return(sample_depth)
  }
)

sample_depth_list <- lapply(sample_depth_list, function(x)
  {
    x <- as.data.frame(x)  
    names(x) <- c("year", "depth")
    return(x)
  }  
)

# Bind the results together into a matrix
first_years <- sapply(sample_depth_list, function(x){min(as.numeric(x$year))})
last_years <- sapply(sample_depth_list, function(x){max(as.numeric(x$year))})

first_year <- min(first_years)
last_year <- max(last_years)

sample_depth_matrix <- matrix(data=0, nrow=n_chronologies, ncol=(last_year-first_year+1))
colnames(sample_depth_matrix) <- as.character(first_year:last_year)

for (i in 1:n_chronologies)
{
  cols_i <- (first_years[i] : last_years[i]) - first_year + 1
  
  # Check that years match
  if(length(cols_i) != nrow(sample_depth_list[[i]]))
  {
    years <- first_years[i]:last_years[i]
    missing_years <- years[!(years %in% sample_depth_list[[i]]$year)]
    extra_rows <- data.frame(year=missing_years, depth=0)
    sample_depth_list[[i]] <- rbind(sample_depth_list[[i]], extra_rows)
    sample_depth_list[[i]] <- sample_depth_list[[i]][order(sample_depth_list[[i]]$year), ]
  }
  
  sample_depth_matrix[i, cols_i] <- sample_depth_list[[i]]$depth
}

sample_depth_df <- cbind(descriptive_stats$id, sample_depth_matrix)

# Plot the graph of sample depth vs. time
summary_sample_depth <- data.frame(depth=colSums(sample_depth_matrix), Year=as.numeric(colnames(sample_depth_matrix)))

sample_depth_plot <- ggplot(summary_sample_depth, aes(y=depth, x=Year)) + geom_area() + theme_bw() + ylab("Number of series analyzed")

# Save results
write.csv(descriptive_stats, file="./Meta-analysis/Results/descriptive/descriptive.csv")
write.csv(sample_depth_df, file="./Meta-analysis/Results/descriptive/sample_depth.csv")
ggsave(filename="./Meta-analysis/Results/descriptive/sample_depth.svg", plot=sample_depth_plot, width=save_width, height=save_height, units="cm")
ggsave(filename="./Meta-analysis/Results/descriptive/sample_depth.pdf", plot=sample_depth_plot, width=save_width, height=save_height, units="cm")

# Model comparison ####
model_comp <- do.call(rbind, 
                      lapply(1:nrow(analysed), 
                             function(x, path, extension){
                               store_fit_and_settings(analysed[x,], path, extension)
                             }
                             , path="./Meta-analysis/Results/all_models/analysis", extension="ALL_MODELS.RData"
                      )
)

# Processing model comparison results
model_comp <- as.data.frame(model_comp)
model_comp[4:14] <- sapply (model_comp[4:14], as.numeric)
model_comp[15:16] <- sapply(model_comp[15:16], unlist)
model_comp[17:19] <- sapply(model_comp[17:19], unlist)
model_comp$code <- 1:16

# Add delta *IC information. Compute relative to full multiplicative model
model_comp <- ddply(model_comp, "file", transform, 
                    dAIC=AIC-AIC[2],
                    dAICc=AICc-AICc[2],
                    dBIC=BIC-BIC[2]
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
mean_class_fit <-   ddply(model_comp, "code", summarise, 
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


sd_class_fit <-   ddply(model_comp, "code", summarise, 
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
write.csv(sd_class_fit, file="./Meta-analysis/Results/all_models/sd_class_fit.csv")

write.csv(selection_occurences, file="./Meta-analysis/Results/all_models/selection_occurences.csv")
write.csv(selection_frequency, file="./Meta-analysis/Results/all_models/selection_frequency.csv")
ggsave(filename="./Meta-analysis/Results/all_models/selection_frequency.svg", plot=sel_freq_plot, width=save_width, height=save_height, units="cm")
ggsave(filename="./Meta-analysis/Results/all_models/selection_frequency.pdf", plot=sel_freq_plot, width=save_width, height=save_height, units="cm")

ggsave(filename="./Meta-analysis/Results/all_models/n_vs_effects.svg", plot=n_vs_effects_plot, width=save_width, height=save_height, units="cm")
ggsave(filename="./Meta-analysis/Results/all_models/n_vs_effects.pdf", plot=n_vs_effects_plot, width=save_width, height=save_height, units="cm")

# Testing for modern sample bias ####
# Load analyses and extract model fit statistics
model_comp_gam <- do.call(rbind, 
                      lapply(1:nrow(analysed), 
                             function(x){
                               store_fit_and_settings(analysed[x,], "./Meta-analysis/Results/gam/analysis", "GAM.RData")
                             }
                      )
)

# Processing model comparison results
model_comp_gam <- as.data.frame(model_comp_gam)
model_comp_gam[4:14] <- sapply (model_comp_gam[4:14], as.numeric)
model_comp_gam[15:16] <- sapply(model_comp_gam[15:16], unlist)
model_comp_gam[17:19] <- sapply(model_comp_gam[17:19], unlist)
model_comp_gam$code <- 1:2


# Add delta *IC information
model_comp_gam <- ddply(model_comp_gam, "file", transform, 
                      dAIC=AIC-AIC[1],
                      dAICc=AICc-AICc[1],
                      dBIC=BIC-BIC[1]
)

# Add *IC model weights
model_comp_gam <- ddply(model_comp_gam, "file", transform, 
                    wAIC=exp(-0.5*dAIC) / sum(exp(-0.5*dAIC)),
                    wAICc=exp(-0.5*dAICc) / sum(exp(-0.5*dAICc)),
                    wBIC=exp(-0.5*dBIC) / sum(exp(-0.5*dBIC))
)

# Show which model is preferred by each metric
model_comp_gam <- ddply(model_comp_gam, "file", transform, 
                    sigma_choice=which.min(sigma),
                    Rsq_choice=which.max(Rsq),
                    adj.Rsq_choice=which.max(adj.Rsq),
                    llh_choice=which.max(llh),
                    AIC_choice=which.min(AIC),
                    AICc_choice=which.min(AICc),
                    BIC_choice=which.min(BIC)
)

# Report mean and sd of fit statistics for each model
mean_class_fit_gam <-   ddply(model_comp_gam, "code", summarise, 
                          sigma=mean(sigma, na.rm=T),
                          Rsq=mean(Rsq, na.rm=T),
                          adj.Rsq=mean(adj.Rsq, na.rm=T),
                          dAIC=mean(dAIC, na.rm=T),
                          dAICc=mean(dAICc, na.rm=T),
                          dBIC=mean(dBIC[is.finite(dBIC)], na.rm=T),
                          wAIC=mean(wAIC, na.rm=T),
                          wAICc=mean(wAICc, na.rm=T),
                          wBIC=mean(wBIC, na.rm=T)
)

sd_class_fit_gam <-   ddply(model_comp_gam, "code", summarise, 
                        sigma=sd(sigma, na.rm=T),
                        Rsq=sd(Rsq, na.rm=T),
                        adj.Rsq=sd(adj.Rsq, na.rm=T),
                        dAIC=sd(dAIC, na.rm=T),
                        dAICc=sd(dAICc, na.rm=T),
                        dBIC=sd(dBIC[is.finite(dBIC)], na.rm=T),
                        wAIC=sd(wAIC, na.rm=T),
                        wAICc=sd(wAICc, na.rm=T),
                        wBIC=sd(wBIC, na.rm=T)
)

# Determine how many times each model was selected by each criteria
selection_occurences_gam <- data.frame( 
  model=1:2,                                
  model_comp_gam[1:2, 15:19],
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
    criteria <- names(selection_occurences_gam)[j]
    selection_occurences_gam[i, j] <- sum(model_comp_gam[criteria]==i)
  }
}

selection_frequency_gam <- selection_occurences_gam
selection_frequency_gam[7:13] <- selection_occurences_gam[7:13]/colSums(selection_occurences_gam[7:13])

# Graph the selection frequency
sel_freq_melt <- melt(selection_frequency_gam, id.vars=c("model", "form", "error", "I", "T", "A"))
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

sel_freq_plot_gam <- ggplot(sel_freq_melt, aes(x=variable, y=code, fill=value)) + geom_tile() + scale_fill_gradient(name="Selection frequency", low="white", high="black") + ylab("Model") + xlab("Selection Criteria")+ylim(levels=code_order)

# Extract the effect vectors of the gam models
extract_effects <- function(data_file,  path=".", extension=".RData")
{
  
  return(effects_list)
}


# Extract effects
# Apply the post hoc correction
gam_effects <- lapply(1:nrow(analysed), function(x, path, extension)
{
  analyses <- load_analysis(analysed[x,], path, extension)
  
  #if(min(as.numeric(names(analyses$gam_ita$effects$T)))==1){
  #  names(analyses$gam_ita$effects$T) <- names(analyses$gam_ta$effects$T)
  #}
  
  analyses$gam_ita$effects <- post_hoc_intercession(analyses$gam_ita$effects, analyses$gam_ita$tra, sparse=TRUE)
  effects_list <- lapply(analyses, function(x) x$effects)   
  return(effects_list)
}
  ,path="./Meta-analysis/Results/gam/analysis", extension="GAM.RData"
)
  
# Calculate the ratio (TA / ITA) of estimated time effects at each year
msb_ratios <- lapply(gam_effects, function(x)
{
  
  T_ita <- x$gam_ita$T
  T_ta <- x$gam_ta$T
  overlap <- intersect(names(T_ita), names(T_ta))
  return(T_ta[overlap]/T_ita[overlap])
}
)

# Check for trends in these ratios
msb_trends <- melt(sapply(msb_ratios, function(x) cor.test(x, 1:length(x), method="spearman")$estimate))

msb_trend_plot <- ggplot(msb_trends, aes(x=value)) + geom_density(fill="black") + theme_bw() + xlab("Spearman's rho test for trend in ratios of G=TA to G=ITA time effects") + ylab("Density")

# Graph the ratio of TA / ITA time effects by year
msb_points <- unlist(msb_ratios)
msb_df <- data.frame(ratio=msb_points, year=as.numeric(names(msb_points)))

# Compute the quantiles at each year
msb_quantiles <- function (year, msb_df){
  slice <- msb_df[msb_df$year==year,]
  
  slice_quantiles <- quantile(slice$ratio, probs=seq(0.2, 0.8, 0.2))
  
  med <- median(slice_quantiles)
  
  out <- c(year=year, median=med, slice_quantiles)
  
  return(out)
  
}

msb_quant_df <- as.data.frame(t(sapply(unique(msb_df$year), msb_quantiles, msb_df=msb_df)))
names(msb_quant_df) <- c("year", "median", "q20", "q40", "q60", "q80")

melted_msb_df <- melt(msb_quant_df, id.vars="year")

msb_plot <- ggplot(msb_quant_df, aes(x=year, y=median)) + geom_ribbon(aes(ymin=q20, ymax=q80), fill="grey75") + geom_ribbon(aes(ymin=q40, ymax=q60), fill="grey50") + geom_line() + theme_bw() + xlab("Year") + ylab("Ratio of G=TA to G=ITA time effects") + geom_hline(y=1)

# Save results
write.csv(model_comp_gam, file="./Meta-analysis/Results/gam/model_comparison_gam.csv")

write.csv(mean_class_fit_gam, file="./Meta-analysis/Results/gam/mean_class_fit_gam.csv")
write.csv(sd_class_fit_gam, file="./Meta-analysis/Results/gam/sd_class_fit_gam.csv")

write.csv(selection_occurences_gam, file="./Meta-analysis/Results/gam/selection_occurences_gam.csv")
write.csv(selection_frequency_gam, file="./Meta-analysis/Results/gam/selection_frequency_gam.csv")
ggsave(filename="./Meta-analysis/Results/gam/selection_frequency_gam.svg", plot=sel_freq_plot_gam, width=save_width, height=save_height, units="cm")
ggsave(filename="./Meta-analysis/Results/gam/selection_frequency_gam.pdf", plot=sel_freq_plot_gam, width=save_width, height=save_height, units="cm")

save(gam_effects, file="./Meta-analysis/Results/gam/gam_effects.RData")
save(msb_ratios, file="./Meta-analysis/Results/gam/msb_ratios.RData")
ggsave(filename="./Meta-analysis/Results/gam/msb_trends.svg", plot=msb_trend_plot, width=save_width, height=save_height, units="cm")
ggsave(filename="./Meta-analysis/Results/gam/msb_trends.pdf", plot=msb_trend_plot, width=save_width, height=save_height, units="cm")
ggsave(filename="./Meta-analysis/Results/gam/msb.svg", plot=msb_plot, width=save_width, height=save_height, units="cm")
ggsave(filename="./Meta-analysis/Results/gam/msb.pdf", plot=msb_plot, width=save_width, height=save_height, units="cm")
ggsave(filename="./Meta-analysis/Results/gam/small_msb.svg", plot=msb_plot + ylim(c(0, 3)) + xlim(c(1000, 2000)), width=save_width, height=save_height, units="cm")
ggsave(filename="./Meta-analysis/Results/gam/small_msb.pdf", plot=msb_plot + ylim(c(0, 3)) + xlim(c(1000, 2000)), width=save_width, height=save_height, units="cm")
