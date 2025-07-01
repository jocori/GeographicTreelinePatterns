setwd("/Users/j364g325/Documents/mountain peaks")
peaks<- read.csv("Peaks.csv", header = TRUE)
head(peaks)
peaks<-na.omit(peaks)
write.csv(peaks, "peaks_clean.csv")
##convert to spatial data
install.packages("sf")

# Load required library
library(sf)


# Convert to a spatial data frame

spatial_peaks <- st_as_sf(peaks, coords = c("long", "lat"), crs = 4326)

# Check the structure
print(spatial_peaks)
st_write(spatial_peaks, "spatial_peaks.shp")
