setwd("~/Desktop/KU/Projects/GeographicTreelinePatterns")
reg<-read.csv("data/Regressions_28Oct.csv")
loc<-reg[,c("Peak","Lat","Long")]
locs<-unique(loc[,"Peak"])
library(dplyr)

# Use dplyr to get one latitude and longitude for each unique Peak
#for all unqiue peaks, 
#get median slope, averaage latitude, and mean longitude for each mountain, 
#then put these summaries as new columns in unique peaks 
unique_peaks <- reg %>%
  filter(Peak != "Mt. Washington") %>% #Remove mt. washington as it's too spatially different from others
  group_by(Peak) %>%
  summarize(Lat = mean(Lat), Long = mean(Long), Slope = median(Slope))

# View the result
print(unique_peaks)





# add column for sign of slope
#if positive -> 1
#if negative -> 0

unique_peaks$sign <- ifelse(unique_peaks$slope > 0, 1,0)
mean(reg$Slope)
mean(unique_peaks$Slope)
write.csv(unique_peaks, "data/figure4_peaks.csv")
