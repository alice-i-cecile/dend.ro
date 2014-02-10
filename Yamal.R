# Loading dependencies ####
library (dplR)
library(reshape)
library(ggplot2)
source("standardization.R")
source("coreTRA.R")

# Grabbing data ####
# Yamal tree ring data from:
# http://www.cru.uea.ac.uk/cru/papers/briffa2013qsr/
# raw.zip
# Reassessing the evidence for tree-growth and inferred temperature change during the Common Era in Yamalia, northwest Siberia
# Keith R. Briffa, Thomas M. Melvin, Timothy J. Osborn, Rashit M. Hantemirov, Alexander V. Kirdyanov, Valeriy Mazepa, Stepan G. Shiyatov and Jan Esper 
# Quaternary Science Reviews (2013) 72, 83-107. doi: 10.1016/j.quascirev.2013.04.008

yamal_trw_path <- "./Examples/yam/yml-all.raw"
yamal_po_path <- "./Examples/yam/yml-all.pth"

yamal_trw <- read.rwl(yamal_trw_path)
yamal_po_raw <- read.csv(yamal_po_path)

# Processing raw pith offset file (nonstandard format)
split_line <- function(line)
{
    verbose_split <- strsplit(as.character(line), split=" ")
    words <- verbose_split[[1]][nchar(verbose_split[[1]])>0]
    return(words)
}

yamal_po_words <- apply(yamal_po_raw, MARGIN=1, split_line)

fill_short_lines <- function(line)
{
  line_length <- length(line)
  
  if (line_length < 7)
  {
    line <- c(line, rep(NA, 7-line_length))
  }
  
  return(line)
}

yamal_po <- as.data.frame(t(sapply(yamal_po_words, fill_short_lines)))
# From readme_raw.prn
# First 3 columns: identifier, estimated year pith grew, estimated pith offset (cm)
names(yamal_po) <- c("ID", "BirthYear", "PO", "Unknown1", "Unknown2", "Unknown3", "Unknown4")

yamal_po$PO <- as.numeric(as.character(yamal_po$PO))

# Pre-processing ####

# Convert to BAI data
bai_series <- function(series_name)
{
  if (series_name %in% yamal_po$ID)
  {
    position <- which(yamal_po$ID==series_name)
    po_df <- data.frame(series=series_name, d2pith=yamal_po$PO[position])
    bai_series <- bai.in(yamal_trw[series_name], po_df)
    return(bai_series)
  } 
  else 
  {
    return (NULL)
  }
}

yamal_bai_raw <- sapply(colnames(yamal_trw), bai_series)
yamal_bai <- combine.rwl(yamal_bai_raw[!sapply(yamal_bai_raw, is.null)])

# Save BAI series
write.rwl(yamal_bai, fname="./Examples/yam/yamal_bai.rwl", format="tucson")

# Extract birth years
birth_years_series <- function(series_name)
{
  if (series_name %in% yamal_po$ID)
  {
    position <- which(yamal_po$ID==series_name)
    birth_year <- as.numeric(as.character(yamal_po$BirthYear[position]))
    return(birth_year)
  } 
  else 
  {
    return ()
  }
}

yamal_birth_years <- sapply(names(yamal_bai), birth_years_series) 

# Convert to sparse tree-ring array
yamal <- rwl.to.stra(yamal_bai, birth_years=yamal_birth_years)

# Standardization ####

# RCS
yamal_rcs <- standardize_tra(yamal, model=list(I=FALSE, T=TRUE, A=TRUE), form="multiplicative", error="lnorm", method="rcs", sparse=TRUE)
save(yamal_rcs, file="./Examples/yam/yml_rcs.RData")
rm (yamal_rcs)

# SFS TA
yamal_sfs_ta <- standardize_tra(yamal, model=list(I=FALSE, T=TRUE, A=TRUE), form="multiplicative", error="lnorm", method="sfs", sparse=TRUE)
save(yamal_sfs_ta, file="./Examples/yam/yml_sfs_ta.RData")
rm(yamal_sfs_ta)

# GAM TA
# Too large for memory
#yamal_gam_ta <- standardize_tra(yamal, model=list(I=FALSE, T=TRUE, A=TRUE), form="multiplicative", error="lnorm", method="gam", sparse=TRUE)
#save(yamal_gam_ta, file="./Examples/yam/yml_gam_ta.RData")
#rm(yamal_gam_ta)

# SFS ITA
yamal_sfs_ita <- standardize_tra(yamal, model=list(I=TRUE, T=TRUE, A=TRUE), form="multiplicative", error="lnorm", method="sfs", sparse=TRUE)
save(yamal_sfs_ita, file="./Examples/yam/yml_sfs_ita.RData")
rm(yamal_sfs_ita)

# GAM ITA
# Too large for memory
#yamal_gam_ita <- standardize_tra(yamal, model=list(I=TRUE, T=TRUE, A=TRUE), form="multiplicative", error="lnorm", method="rcs", sparse=TRUE) 
#save(yamal_gam_ita, file="./Examples/yam/yml_gam_ita.RData")
#rm(yamal_gam_ita)

# Analysis ####

# Reload standardizations
load("./Examples/yam/yml_rcs.RData")
load("./Examples/yam/yml_sfs_ta.RData")
load("./Examples/yam/yml_sfs_ita.RData")

# Combine fit data
yamal_fit <- as.data.frame(t(data.frame(
    rcs=unlist(yamal_rcs$fit[3:13]),
    sfs_ta=unlist(yamal_sfs_ta$fit[3:13]),
#    gam_ta=unlist(yamal_gam_ta$fit[3:13]),
    sfs_ita=unlist(yamal_sfs_ita$fit[3:13])#,
#    gam_ita=unlist(yamal_gam_ita$fit[3:13])
    )
  )
)
  
# Save and display fit data
print(yamal_fit)

write.csv(yamal_fit, file="./Examples/yam/yml_fit.csv")

# Combine effects
yamal_effects_raw <- list(
  rcs=yamal_rcs$effects,
  sfs_ta=yamal_sfs_ta$effects,
#  gam_ta=yamal_gam_ta$effects,
  sfs_ita=yamal_sfs_ita$effects#,
#  gam_ita=yamal_gam_ita$effects
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

yamal_effects_list <- list()
for (i in 1:length(yamal_effects_raw))
{
  yamal_effects_list[[i]] <- dfify_effects(yamal_effects_raw[[i]])
}
names(yamal_effects_list) <- names(yamal_effects_raw)

for (i in 1:length(yamal_effects_list))
{
    yamal_effects_list[[i]]$model <- names(yamal_effects_list)[i] 
}

yamal_effects <- head(yamal_effects_list[[1]], 0)
for (i in 1:length(yamal_effects_list))
{
  effects_i <- yamal_effects_list[[i]]
  yamal_effects <- rbind(yamal_effects, effects_i)
}

# Save effects
write.csv(yamal_effects, file="./Examples/yam/yml_effects.csv")

# Plotting ####

# Plot sample depth
yamal_sample_depth <- sample_depth_tra(yamal, factor.dim=2, sparse=TRUE)

yamal_sd_plot <- ggplot(data.frame(sd=yamal_sample_depth, year=as.numeric(names(yamal_sample_depth))), aes(x=year, y=sd)) + geom_area() + theme_bw()
print(yamal_sd_plot)

# Plot I
I_effects <- yamal_effects[yamal_effects$type=="I"&yamal_effects$model=="sfs_ita",]

yamal_I_hist <- ggplot(I_effects, aes(x=effect, colour=model, fill=model)) + geom_density(alpha=0.7) + theme_bw()
print (yamal_I_hist)

# Plot I vs. birth year
I_year_df <- I_effects

match_birth_year <- function(r){
  id <- I_year_df[r, "id"]
  birth_year <- yamal_birth_years[names(yamal_birth_years)==id]
  return(birth_year)
}

I_year_df$birth_year <- sapply(1:nrow(I_year_df), match_birth_year)

yamal_I_birth_plot <- ggplot(I_year_df, aes(x=birth_year, y=effect, colour=model)) + geom_point() + geom_smooth()
print(yamal_I_birth_plot)

# Correlation test for MSB
yamal_I_birth_trend <- cor.test(x=I_year_df$birth_year, y=I_year_df$effect, method="spearman")
print(yamal_I_birth_trend)
save(yamal_I_birth_trend, file="./Examples/yam/yml_I_birth_trend.RData")

# Plot T
T_effects <- yamal_effects[yamal_effects$type=="T",]
T_effects$id <- as.numeric(as.character(T_effects$id))

yamal_T_plot <- ggplot(T_effects, aes(y=effect, x=id, colour=model)) + geom_line() + facet_grid(model~.) + theme_bw()
print(yamal_T_plot)

# Plot MSB distortion and SFS effects
yamal_rcs <- T_effects[T_effects$model=="rcs", ]
yamal_rcs <- yamal_rcs[order(yamal_rcs$id),]

yamal_sfs_ta <- T_effects[T_effects$model=="sfs_ta", ]
yamal_sfs_ta <- yamal_sfs_ta[order(yamal_sfs_ta$id),]

yamal_sfs_ita <- T_effects[T_effects$model=="sfs_ita", ]
yamal_sfs_ita <- yamal_sfs_ita[order(yamal_sfs_ita$id),]

yamal_sfs_effect <- yamal_rcs
yamal_sfs_effect$effect <- yamal_sfs_ta$effect / yamal_rcs$effect

yamal_msb_effect <- yamal_rcs
yamal_msb_effect$effect <- yamal_sfs_ita$effect / yamal_sfs_ta$effect

write.csv(yamal_sfs_effect, file="./Examples/yam/yamal_sfs_effect.csv")
write.csv(yamal_msb_effect, file="./Examples/yam/yamal_msb_effect.csv")


yamal_sfs_effect_plot <- ggplot(yamal_sfs_effect, aes(y=effect, x=id)) + geom_line() + geom_hline(y=1)
print(yamal_sfs_effect_plot)

yamal_msb_effect_plot <- ggplot(yamal_msb_effect, aes(y=effect, x=id)) + geom_line() + geom_hline(y=1)
print(yamal_msb_effect_plot)

# Plot A
A_effects <- yamal_effects[yamal_effects$type=="A",]
A_effects$id <- as.numeric(as.character(A_effects$id))

yamal_A_plot <- ggplot(A_effects, aes(y=effect, x=id, colour=model)) + geom_line() + facet_grid(model~.) + theme_bw()
print(yamal_A_plot)