# Load data from each example ####
yamal_effects <- read.csv(file="./Examples/yam/yml_effects.csv", header=TRUE)
torn_effects <- read.csv(file="./Examples/torn/trn_effects.csv", header=TRUE)

# Don't load in rownames as a new colummn
yamal_effects <- yamal_effects[,2:ncol(yamal_effects)]
torn_effects <- torn_effects[,2:ncol(torn_effects)]

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
  if (series_name %in% yamal_po$ID)
  {
    position <- which(yamal_po$ID==series_name)
    birth_year <- as.numeric(as.character(torn_po$BirthYear[position]))
    return(birth_year)
  } 
  else 
  {
    return ()
  }
}

yamal_birth_years <- sapply(names(yamal_trw), birth_years_series_YAMAL)
torn_birth_years <- sapply(names(torn_trw),birth_years_series_TORN)

match_birth_year <- function(r, birth_years, df){
  id <- df[r, "id"]
  birth_year <- birth_years[names(birth_years)==id]
  return(birth_year)
}

yamal_I_year_df$birth_year <- sapply(1:nrow(yamal_I_year_df), match_birth_year, birth_years=yamal_birth_years, df=yamal_I_year_df)
torn_I_year_df$birth_year <- sapply(1:nrow(torn_I_year_df), match_birth_year, birth_years=torn_birth_years, df=torn_I_year_df)

yamal_I_year_df$Chronology <- "Yamal"
torn_I_year_df$Chronology <- "Tornetrask"
I_year_df <- rbind(yamal_I_year_df, torn_I_year_df)


yamal_I_birth_plot <- ggplot(I_year_df, aes(x=birth_year, y=effect)) + geom_point() + geom_smooth(colour="red") + ylab("Individual effect") + xlab("Year of birth") + theme_bw()
print(yamal_I_birth_plot)


# Plot time signals ####

# Plot MSB effect ####