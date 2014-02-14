# Libraries ####
library(dplR)
library(ggplot2)

# Load data from each example ####
yamal_effects <- read.csv(file="./Examples/yam/yml_effects.csv", header=TRUE)
torn_effects <- read.csv(file="./Examples/torn/trn_effects.csv", header=TRUE)

yamal_trw_path <- "./Examples/yam/yml-all.raw"
torn_trw_path <- "./Examples/torn/alltrw.raw"
yamal_trw <- read.rwl(yamal_trw_path)
torn_trw <- read.rwl(torn_trw_path)

yamal_po <- read.csv("./Examples/yam/yml_po.csv")
torn_po <- read.csv("./Examples/torn/torn_po.csv")

# Splitting nicely
yamal_I <- yamal_effects[yamal_effects$type=="I" & (yamal_effects$model=="ITA" | yamal_effects$model=="TA"),]
yamal_T <-  yamal_effects[yamal_effects$type=="T" & (yamal_effects$model=="ITA" | yamal_effects$model=="TA"),]
yamal_A <-  yamal_effects[yamal_effects$type=="A" & (yamal_effects$model=="ITA" | yamal_effects$model=="TA"),]

torn_I <- torn_effects[torn_effects$type=="I" & (torn_effects$model=="ITA" | torn_effects$model=="TA"),]
torn_T <-  torn_effects[torn_effects$type=="T" & (torn_effects$model=="ITA" | torn_effects$model=="TA"),]
torn_A <-  torn_effects[torn_effects$type=="A" & (torn_effects$model=="ITA" | torn_effects$model=="TA"),]

# Plot I vs. age ####

# Extract date of birth
yamal_I_year_df <- yamal_I
torn_I_year_df <- torn_I

birth_years_series_YAMAL <- function(series_name)
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

birth_years_series_TORN <- function(series_name)
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

yamal_birth_years <- unlist(sapply(names(yamal_trw), birth_years_series_YAMAL))
yamal_birth_years <- data.frame(id=names(yamal_birth_years), birth_year=yamal_birth_years)
torn_birth_years <- unlist(sapply(names(torn_trw), birth_years_series_TORN))
torn_birth_years <- data.frame(id=names(torn_birth_years), birth_year=torn_birth_years)

yamal_I_year_df <- merge(yamal_I_year_df, yamal_birth_years, by="id")
torn_I_year_df <- merge(torn_I_year_df, torn_birth_years, by="id")

# Combining data
yamal_I_year_df$Chronology <- "Yamal"
torn_I_year_df$Chronology <- "Tornetrask"
I_year_df <- rbind(yamal_I_year_df, torn_I_year_df)

I_birth_plot <- ggplot(I_year_df, aes(x=birth_year, y=effect)) + geom_point() + geom_smooth(colour="red") + ylab("Individual effect") + xlab("Year of birth") + theme_bw() + facet_grid(Chronology~.)
print(yamal_I_birth_plot)

# Plot time signals ####
yamal_T$id <- as.numeric(as.character(yamal_T$id))
torn_T$id <- as.numeric(as.character(torn_T$id))

yamal_T$Chronology <- "Yamal"
torn_T$Chronology <- "Tornetrask"

joint_T <- rbind(yamal_T, torn_T)
T_plot <- ggplot(joint_T, aes(y=effect, x=id, colour=model)) + geom_line(alpha=0.5) + theme_bw() + geom_hline(y=1) + ylab("Time effect") + xlab("Year") + theme_bw() + scale_colour_manual(values=c("black", "red"), guide="none") + facet_grid(Chronology~.) + ylim (c(0, 3))
print(T_plot) 

# Plot MSB effect ####
joint_T_TA <- joint_T[joint_T$model=="TA",]
joint_T_ITA <- joint_T[joint_T$model=="ITA",]

msb_effect <- joint_T_TA
msb_effect$effect <- joint_T_TA$effect / joint_T_ITA$effect

msb_effect_plot <- ggplot(msb_effect, aes(y=effect, x=id)) + geom_line() + geom_hline(y=1) + ylab("Ratio between uncorrected and corrected time effect") + xlab("Year") + theme_bw() + facet_grid(Chronology~.)
print(msb_effect_plot)

