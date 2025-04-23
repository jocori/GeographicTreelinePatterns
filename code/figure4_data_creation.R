setwd("~/Desktop/KU/Projects/GeographicTreelinePatterns")
reg<-read.csv("data/Regressions_28Oct.csv")
loc<-reg[,c("Peak","Lat","Long")]
locs<-unique(loc[,"Peak"])
library(dplyr)

# Use dplyr to get one latitude and longitude for each unique Peak
unique_peaks <- reg %>%
  group_by(Peak) %>%
  summarize(Lat = first(Lat), Long = first(Long), Slope = first(Slope))

# View the result
print(unique_peaks)
#Remove mt. washington as it's too spatially different from others
unique_peaks<-unique_peaks[unique_peaks$Peak != "Mt. Washington",]

# add column for sign of slope
#if positive -> 1
#if negative -> 0

reg$sign <- ifelse(reg$Slope > 0, 1,0)


write.csv(unique_peaks, "data/figure4_peaks.csv")
