# Set RNG seed
set.seed(42)

# Dependencies ####
library (dplR)
library(reshape)
library(ggplot2)
source("standardization.R")
source("coreTRA.R")
source("syntheticData.R")

# Generating time effect ####
low_freq_T <-  0.3*sin((2*pi/500)*1:2000)
high_freq_T <- 0.2*sin((2*pi/50)*1:2000)
noise_T <- rnorm(2000, sd=0.1)
synth_T <- exp(low_freq_T + high_freq_T + noise_T)

# Scaling to geometric mean of 1
synth_T <- synth_T / exp(mean(log(synth_T)))
names(synth_T) <- 1:2000
                                                  
# Generating age effect ####
synth_A <- constBAI_trend(x=1:300, k=50)
names(synth_A) <- 1:300
                         
plot(synth_A)
plot(cumsum(synth_A))

# Generating individual effect and DOB for each tree ####

# Ancient trees are easy                         
ancient_I <- rlnorm(200, sdlog=0.2)
names(ancient_I) <- paste0("A", 1:200)
                         
ancient_DOB <- sample(1:1699, 200, replace=TRUE)
ancient_DOD <- ancient_DOB + 300
                         
# Modern trees suffer from big tree selection bias
generate_modern <- function (cutoff, last_year, I_effect=rlnorm(1), T_effect, A_effect, DOB=sample(1700:1999, 1)){
  T_trunc <- T_effect[DOB:last_year]
  A_trunc <- A_effect[1:length(T_trunc)]
  
  rw <- I_effect*T_trunc*A_trunc
  
  dbh <- sum(rw, na.rm=TRUE)
  print(dbh)
  if(dbh>=cutoff){
    return(c(I=I_effect, DOB=DOB))
  } else {
    return()
  }
}

# Only include trees that pass the diameter cutoff
modern_summary <- data.frame(I=NA, DOB=NA)[0,]
                         
while(nrow(modern_summary)<200)
{
  test_tree <- generate_modern(50, 2000, rlnorm(1, sdlog=0.2), synth_T, synth_A)
  modern_summary <- rbind(modern_summary, test_tree)
}

names(modern_summary) <- c("I", "DOB")
modern_I <- modern_summary$I
names(modern_I) <- paste0("M", 1:200)

modern_DOB <- modern_summary$DOB
names(modern_DOB) <- paste0("M", 1:200)
                         
modern_DOD <- rep(2000, 200)
names(modern_DOD) <- paste0("M", 1:200)                         
                         
# Generate tree ring data ####
                         
ancient <- data.frame(G=NA, i=NA, t=NA, a=NA)[0,]
# Twenty at a time for memory
for (i in seq(from=1, to=200, by=20))
{
  print(i)
  trees_i <- base_synth_TRA(ancient_I[i:(i+19)], synth_T, synth_A, ancient_DOB[i:(i+19)], ancient_DOD[i:(i+19)], 0.3, TRUE)
  ancient <- rbind(ancient, sparse_tra(trees_i))
}
         
modern <- data.frame(G=NA, i=NA, t=NA, a=NA)[0,]
# Twenty at a time for memory
for (i in seq(from=1, to=200, by=20))
{
   print(i)
   trees_i <- base_synth_TRA(modern_I[i:(i+19)], synth_T, synth_A, modern_DOB[i:(i+19)], modern_DOD[i:(i+19)], 0.3, TRUE)
   modern <- rbind(modern, sparse_tra(trees_i))
}

demo <- rbind(ancient, modern)
                         
# Analyze the data ####
 
# Signal-free regional curve standardization
demo_sfs <- standardize_tra(demo, model=list(I=FALSE, T=TRUE, A=TRUE), form="multiplicative", error="lnorm", method="sfs", sparse=TRUE)
T_sfs <- demo_sfs$effects$T
ratio_sfs <- T_sfs / synth_T [(2001-length(T_sfs)):2000]

# FES with I
demo_fes <- standardize_tra(demo, model=list(I=TRUE, T=TRUE, A=TRUE), form="multiplicative", error="lnorm", method="sfs", sparse=TRUE)
T_fes <- demo_fes$effects$T
ratio_fes <- T_fes / synth_T [(2001-length(T_fes)):2000]

# Spline detrending
demo_rwl <- tra.to.rwl(unsparse_tra(demo))
demo_iss <- chron(detrend(demo_rwl, method="Spline"), biweight=F, prewhiten=F)
T_iss <- demo_iss$xxxstd
T_iss <- T_iss / exp(mean(log(T_iss)))
ratio_iss <- T_iss / synth_T [(2001-length(T_iss)):2000]
                         
# Plotting ####

# Overlaid true and reconstruction
T_df <- data.frame(T_cv=c(T_sfs, T_fes, T_iss),
t=rep(as.numeric(as.character(names(T_sfs))), 3),
model=c(rep("SFS-RCS", length(T_sfs)),
rep("FES", length(T_fes)),
rep("ISS-Spline", length(T_iss))                  
))
T_df$model <- factor(x=T_df$model, levels=c("ISS-Spline","SFS-RCS","FES"))
                         
true_T <- data.frame(T_cv=synth_T [(2001-length(T_fes)):2000], t=as.numeric(as.character(names(T_sfs))))

ggplot(T_df, aes(x=t, y=T_cv)) + geom_line(alpha=0.5, colour="red") + facet_grid(model~.) + ylim (c(0,2)) + geom_line(data=true_T, colour="black", alpha=0.5) + theme_bw() + xlab("Year") + ylab("Time effect") + geom_vline(x=1700) + geom_hline(y=1)
                         
# Ratio between true and reconstruction
                         
ratio_df <- data.frame(T_cv=c(ratio_sfs, ratio_fes, ratio_iss),
t=rep(as.numeric(as.character(names(T_sfs))), 3),
model=c(rep("SFS-RCS", length(T_sfs)),
rep("FES", length(T_fes)),
rep("ISS-Spline", length(T_iss))                  
))
ratio_df$model <- factor(x=ratio_df$model, levels=c("ISS-Spline","SFS-RCS","FES"))

                         
ggplot(ratio_df, aes(x=t, y=T_cv)) + geom_line() + facet_grid(model~.) + ylim (c(0,2)) + theme_bw() + xlab("Year") + ylab("Ratio between reconstructed and true time effect") + geom_vline(x=1700) + geom_hline(y=1)

# Saving ####
write.csv(demo, file="./Examples/synthetic/demo_chron.csv")
write.csv(T_df, file="./Examples/synthetic/time_effects.csv")
write.csv(ratio_df, file="./Examples/synthetic/ratios.csv")


