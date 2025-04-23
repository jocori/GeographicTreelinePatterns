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
#for all unqiue peaks, 
#get average slope, latitude, and longitude for each mountain, 
#then put these averages as new columns in unique peaks 
# compute the averages from the full `reg` data
avg_stats <- reg %>%
  filter(Peak != "Mt. Washington") %>%
  group_by(Peak) %>%
  summarise(
    slope     = mean(Slope,     na.rm = TRUE),
    latitude  = mean(Lat,  na.rm = TRUE),
    longitude = mean(Long, na.rm = TRUE)
  )

# join them back onto your unique_peaks
unique_peaks <- unique_peaks %>%
  left_join(avg_stats, by = "Peak")
# add column for sign of slope
#if positive -> 1
#if negative -> 0

unique_peaks$sign <- ifelse(unique_peaks$Slope > 0, 1,0)


write.csv(unique_peaks, "data/figure4_peaks.csv")
