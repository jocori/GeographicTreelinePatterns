setwd("~/Desktop/KU/Projects/GeographicTreelinePatterns")
reg<-read.csv("data/Regressions_3Jan.csv")
loc<-reg[,c("Peak","Lat","Long")]
locs<-unique(loc[,"Peak"])
# Install and load dplyr if not already installed
# install.packages("dplyr")
library(dplyr)

# Use dplyr to get one latitude and longitude for each unique Peak
unique_peaks <- reg %>%
  group_by(Peak) %>%
  summarize(Lat = first(Lat), Long = first(Long))

# View the result
print(unique_peaks)

write.csv(unique_peaks, "data/figure1_peaks.csv")

#Remove mt. washington as it's too spatially different from others
unique_peaks<-unique_peaks[unique_peaks$Peak != "Mt. Washington",]
write.csv(unique_peaks, "data/figure1_peaks.csv")
