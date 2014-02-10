# Loading dependencies ####
library (dplR)
library(reshape)
library(ggplot2)
source("standardization.R")
source("coreTRA.R")

# Grabbing data ####
# Tornetrask tree ring data
# http://www.cru.uea.ac.uk/cru/papers/melvin2012holocene/
# Potential bias in 'updating' tree-ring chronologies using regional curve standardisation: 
# Re-processing 1500 years of Torneträsk density and ring-width data

# Thomas M Melvin, Håkan Grudd and Keith R Briffa 
# The Holocene (2013) 23, 364-373. 
# DOI: 10.1177/0959683612460791

torn_trw_path <- "./Examples/torn/alltrw.raw"
torn_po_path <- "./Examples/torn/alltrw.pth"

torn_trw <- read.rwl(torn_trw_path)
torn_po_raw <- read.csv(torn_po_path)

# Processing raw pith offset file (nonstandard format)
split_line <- function(line)
{
  verbose_split <- strsplit(as.character(line), split=" ")
  words <- verbose_split[[1]][nchar(verbose_split[[1]])>0]
  return(words)
}

torn_po <- as.data.frame(t(apply(torn_po_raw, MARGIN=1, split_line)))

# From readme_raw.prn
# First 3 columns: identifier, estimated year pith grew, estimated pith offset (cm)
names(torn_po) <- c("ID", "BirthYear", "PO")

torn_po$PO <- as.numeric(as.character(torn_po$PO))

# Pre-processing ####

# Convert to BAI data
bai_series <- function(series_name)
{
  if (series_name %in% torn_po$ID)
  {
    position <- which(torn_po$ID==series_name)
    po_df <- data.frame(series=series_name, d2pith=torn_po$PO[position])
    bai_series <- bai.in(torn_trw[series_name], po_df)
    return(bai_series)
  } 
  else 
  {
    return (NULL)
  }
}

torn_bai_raw <- sapply(colnames(torn_trw), bai_series)
torn_bai <- combine.rwl(torn_bai_raw[!sapply(torn_bai_raw, is.null)])

# Save BAI series
write.rwl(torn_bai, fname="./Examples/torn/torn_bai.rwl", format="tucson")

# Extract birth years
birth_years_series <- function(series_name)
{
  if (series_name %in% torn_po$ID)
  {
    position <- which(torn_po$ID==series_name)
    birth_year <- as.numeric(as.character(torn_po$BirthYear[position]))
    return(birth_year)
  } 
  else 
  {
    return ()
  }
}

torn_birth_years <- sapply(names(torn_bai), birth_years_series) 

# Convert to sparse tree-ring array
torn <- rwl.to.stra(torn_bai, birth_years=torn_birth_years)

# Standardization ####

# RCS
torn_rcs <- standardize_tra(torn, model=list(I=FALSE, T=TRUE, A=TRUE), form="multiplicative", error="lnorm", method="rcs", sparse=TRUE)
save(torn_rcs, file="./Examples/torn/trn_rcs.RData")
rm (torn_rcs)

# SFS TA
torn_sfs_ta <- standardize_tra(torn, model=list(I=FALSE, T=TRUE, A=TRUE), form="multiplicative", error="lnorm", method="sfs", sparse=TRUE)
save(torn_sfs_ta, file="./Examples/torn/trn_sfs_ta.RData")
rm(torn_sfs_ta)

# GAM TA
# Too large for memory
#torn_gam_ta <- standardize_tra(torn, model=list(I=FALSE, T=TRUE, A=TRUE), form="multiplicative", error="lnorm", method="gam", sparse=TRUE)
#save(torn_gam_ta, file="./Examples/torn/trn_gam_ta.RData")
#rm(torn_gam_ta)

# SFS ITA
torn_sfs_ita <- standardize_tra(torn, model=list(I=TRUE, T=TRUE, A=TRUE), form="multiplicative", error="lnorm", method="sfs", sparse=TRUE)
save(torn_sfs_ita, file="./Examples/torn/trn_sfs_ita.RData")
rm(torn_sfs_ita)

# GAM ITA
# Too large for memory
#torn_gam_ita <- standardize_tra(torn, model=list(I=TRUE, T=TRUE, A=TRUE), form="multiplicative", error="lnorm", method="rcs", sparse=TRUE) 
#save(torn_gam_ita, file="./Examples/torn/trn_gam_ita.RData")
#rm(torn_gam_ita)

# Analysis ####

# Reload standardizations
load("./Examples/torn/trn_rcs.RData")
load("./Examples/torn/trn_sfs_ta.RData")
load("./Examples/torn/trn_sfs_ita.RData")

# Combine fit data
torn_fit <- as.data.frame(t(data.frame(
  rcs=unlist(torn_rcs$fit[3:13]),
  sfs_ta=unlist(torn_sfs_ta$fit[3:13]),
  #    gam_ta=unlist(torn_gam_ta$fit[3:13]),
  sfs_ita=unlist(torn_sfs_ita$fit[3:13])#,
  #    gam_ita=unlist(torn_gam_ita$fit[3:13])
)
)
)

# Save and display fit data
print(torn_fit)

write.csv(torn_fit, file="./Examples/torn/trn_fit.csv")

# Combine effects
torn_effects_raw <- list(
  rcs=torn_rcs$effects,
  sfs_ta=torn_sfs_ta$effects,
  #  gam_ta=torn_gam_ta$effects,
  sfs_ita=torn_sfs_ita$effects#,
  #  gam_ita=torn_gam_ita$effects
)

dfify_effects <- function(effects_list)
{
  dfify_effect <- function(effect)
  {
    effect_df <- data.frame(id=names(unlist(effect)), effect=unlist(effect))
    return(effect_df)
  }
  
  effects <- lapply(effects_list, dfify_effect)
  
  for (i in 1:length(names(effects)))
  {
    effects[[i]]$type <- names(effects)[i]
  }
  
  effects <- merge_all(effects)
  
  return (effects)
}

torn_effects_list <- list()
for (i in 1:length(torn_effects_raw))
{
  torn_effects_list[[i]] <- dfify_effects(torn_effects_raw[[i]])
}
names(torn_effects_list) <- names(torn_effects_raw)

for (i in 1:length(torn_effects_list))
{
  torn_effects_list[[i]]$model <- names(torn_effects_list)[i] 
}

torn_effects <- head(torn_effects_list[[1]], 0)
for (i in 1:length(torn_effects_list))
{
  effects_i <- torn_effects_list[[i]]
  torn_effects <- rbind(torn_effects, effects_i)
}

# Save effects
write.csv(torn_effects, file="./Examples/torn/trn_effects.csv")

# Plotting ####

# Plot sample depth
torn_sample_depth <- sample_depth_tra(torn, factor.dim=2, sparse=TRUE)

torn_sd_plot <- ggplot(data.frame(sd=torn_sample_depth, year=as.numeric(names(torn_sample_depth))), aes(x=year, y=sd)) + geom_area() + theme_bw()
print(torn_sd_plot)

# Plot I
I_effects <- torn_effects[torn_effects$type=="I"&torn_effects$model=="sfs_ita",]

torn_I_hist <- ggplot(I_effects, aes(x=effect, colour=model, fill=model)) + geom_density(alpha=0.7) + theme_bw()
print (torn_I_hist)

# Plot I vs. birth year
I_year_df <- I_effects

match_birth_year <- function(r){
  id <- I_year_df[r, "id"]
  birth_year <- torn_birth_years[names(torn_birth_years)==id]
  return(birth_year)
}

I_year_df$birth_year <- sapply(1:nrow(I_year_df), match_birth_year)

torn_I_birth_plot <- ggplot(I_year_df, aes(x=birth_year, y=effect, colour=model)) + geom_point() + geom_smooth()
print(torn_I_birth_plot)

# Correlation test for MSB
torn_I_birth_trend <- cor.test(x=I_year_df$birth_year, y=I_year_df$effect, method="spearman")
print(torn_I_birth_trend)
save(torn_I_birth_trend, file="./Examples/torn/trn_I_birth_trend.RData")

# Plot T
T_effects <- torn_effects[torn_effects$type=="T",]
T_effects$id <- as.numeric(as.character(T_effects$id))

# Truncate around Year 0
T_effects <- T_effects[T_effects$id > 0,]

torn_T_plot <- ggplot(T_effects, aes(y=effect, x=id, colour=model)) + geom_line() + facet_grid(model~.) + theme_bw()
print(torn_T_plot)

# Plot MSB distortion and SFS effects
torn_rcs <- T_effects[T_effects$model=="rcs", ]
torn_rcs <- torn_rcs[order(torn_rcs$id),]

torn_sfs_ta <- T_effects[T_effects$model=="sfs_ta", ]
torn_sfs_ta <- torn_sfs_ta[order(torn_sfs_ta$id),]

torn_sfs_ita <- T_effects[T_effects$model=="sfs_ita", ]
torn_sfs_ita <- torn_sfs_ita[order(torn_sfs_ita$id),]

torn_sfs_effect <- torn_rcs
torn_sfs_effect$effect <- torn_sfs_ta$effect / torn_rcs$effect

torn_msb_effect <- torn_rcs
torn_msb_effect$effect <- torn_sfs_ita$effect / torn_sfs_ta$effect

write.csv(torn_sfs_effect, file="./Examples/torn/torn_sfs_effect.csv")
write.csv(torn_msb_effect, file="./Examples/torn/torn_msb_effect.csv")


torn_sfs_effect_plot <- ggplot(torn_sfs_effect, aes(y=effect, x=id)) + geom_line() + geom_hline(y=1)
print(torn_sfs_effect_plot)

torn_msb_effect_plot <- ggplot(torn_msb_effect, aes(y=effect, x=id)) + geom_line() + geom_hline(y=1)
print(torn_msb_effect_plot)

# Plot A
A_effects <- torn_effects[torn_effects$type=="A",]
A_effects$id <- as.numeric(as.character(A_effects$id))

torn_A_plot <- ggplot(A_effects, aes(y=effect, x=id, colour=model)) + geom_line() + facet_grid(model~.) + theme_bw()
print(torn_A_plot)