# Load data and libraries ####

library(ggplot2)

load("./Meta-analysis/Meta/gam/msb_ratios.RData")

# Extract ages and trends ####

# First year
chron_start <- sapply(msb_ratios, function(x){
  min(as.numeric(names(x)))
})


# Trend
chron_trend <- sapply(msb_ratios, function(x){
  cor.test(x, as.numeric(names(x)), method="spearman")$estimate
})



# Combining findings
trend_df <- data.frame(start=chron_start, trend=chron_trend)

trend_df$group <- sapply(chron_start, function(x){
  if (x < 1700){
    return("Ancient")
  } else {
    return("Modern")
  }
})


# Average bias by year for each group
ratio_group <- lapply(msb_ratios, function(x){
  if(min(as.numeric(names(x))) < 1700)
  {
    df <- data.frame(ratio=x, year=as.numeric(names(x)), group="Ancient")
  } else {
    df <- data.frame(ratio=x, year=as.numeric(names(x)), group="Modern")
  }
    return (df)
})

ratio_df <- Reduce(rbind, ratio_group)

# Plotting ####

ggplot(trend_df, aes(x=start, y=trend)) + geom_point() + geom_smooth()

ggplot(trend_df, aes(x=trend, colour=group, fill=group)) + geom_density(alpha=0.5) + theme_bw()

ggplot(ratio_df, aes(x=year, y=ratio, colour=group)) + geom_density2d(size=1) + xlab("Year") + ylab("Ratio between uncorrected and corrected chronologies") + theme_bw() + scale_colour_manual(values=c("dodgerblue3", "orange"), guide=FALSE) + ylim(c(0,2)) + xlim(c(1500, 2015)) + theme(panel.background = element_rect(fill="grey85")) + geom_hline(y=1)