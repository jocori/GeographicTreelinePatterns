setwd("C:/Users/j364g325/Documents/mountain peaks")
regs <- read.csv("data/Regressions_3Jan.csv")
install.packages("dplyr")
library(dplyr)
install.packages("tidyr")
library(tidyr)
##make database for town's map figure
averaged_data <- regs %>%
  group_by(Peak, Year) %>%
  summarize(
    Significance = mean(Significance, na.rm = TRUE),
    Intercept = mean(Intercept, na.rm = TRUE),
    Slope = mean(Slope, na.rm = TRUE),
    AIC = mean(AIC, na.rm = TRUE),
    Stations_After_Treeline = mean(Stations_After_Treeline, na.rm = TRUE),
    max_NDVI = mean(max_NDVI, na.rm = TRUE),
    min_NDVI = mean(min_NDVI, na.rm = TRUE),
    Peak_ID = first(Peak_ID),  # Assuming Peak_ID is constant for each Peak across Directions
    Long = first(Long),        # Assuming Long is constant for each Peak across Directions
    Lat = first(Lat),          # Assuming Lat is constant for each Peak across Directions
    Country = first(Country),  # Assuming Country is constant for each Peak across Directions
    Regression_Type = first(Regression_Type),  # Assuming you might still want this info, taking the first one as representative
    Best_Model = first(Best_Model),  # Assuming you might still want this info, taking the first one as representative
    .groups = "drop"  # Removes grouping
  )
## find difference in max NDVI
# Step 1: Filter data for the two specific years
filtered_data <- regs %>%
  filter(Year %in% c("Avg1317", "Avg8488")) %>%
  select(Peak, Year, max_NDVI, Lat, Long)

# Step 2: Compute the average max_NDVI for each Peak within the filtered years
avg_ndvi_by_peak <- filtered_data %>%
  group_by(Peak, Year) %>%
  summarize(avg_max_NDVI = mean(max_NDVI, na.rm = TRUE), .groups = "drop")

# Step 3: Spread the data so each Peak has Avg_1317 and Avg_8488 in separate columns
avg_ndvi_wide <- avg_ndvi_by_peak %>%
  pivot_wider(names_from = Year, values_from = avg_max_NDVI)

# Step 4: Calculate the change in max_NDVI between Avg_1317 and Avg_8892 for each Peak
avg_ndvi_wide <- avg_ndvi_wide %>%
  mutate(Change_in_max_NDVI = `Avg1317` - `Avg8488`)

# View the results
print(avg_ndvi_wide)

## find difference in min NDVI
# Step 1: Filter data for the two specific years
filtered_data_min <- regs %>%
  filter(Year %in% c("Avg1317", "Avg8488")) %>%
  select(Peak, Year, min_NDVI)

# Step 2: Compute the average max_NDVI for each Peak within the filtered years
avg_min_ndvi_by_peak <- filtered_data_min %>%
  group_by(Peak, Year) %>%
  summarize(avg_min_NDVI = mean(min_NDVI, na.rm = TRUE), .groups = "drop")

# Step 3: Spread the data so each Peak has Avg_1317 and Avg_8488 in separate columns
avg_min_ndvi_wide <- avg_min_ndvi_by_peak %>%
  pivot_wider(names_from = Year, values_from = avg_min_NDVI)

# Step 4: Calculate the change in max_NDVI between Avg_1317 and Avg_8892 for each Peak
avg_min_ndvi_wide <- avg_min_ndvi_wide %>%
  mutate(Change_in_min_NDVI = `Avg1317` - `Avg8488`)

# View the results
print(avg_min_ndvi_wide)


# Assuming avg_ndvi_wide and averaged_data are your data frames
# and both have a 'Peak' column to join on

# Add Lat and Long from averaged_data to avg_ndvi_wide
avg_ndvi_wide_with_coords <- avg_ndvi_wide %>%
  left_join(averaged_data %>% select(Peak, Long, Lat) %>% distinct(Peak, .keep_all = TRUE), by = "Peak")

# Check the first few rows to confirm Lat and Long are added
head(avg_ndvi_wide_with_coords)
#add color for map plot
avg_ndvi_wide_with_coords$color_indicator <- with(avg_ndvi_wide_with_coords, ifelse(Change_in_max_NDVI > 0, "Positive", ifelse(Change_in_max_NDVI < 0, "Negative", "Neutral")))
#import map .shp
install.packages("ggplot2")
library(ggplot2)
install.packages("maps")
library(maps)
install.packages("crayon")
library(crayon)
# Get North America map data
north_america_map <- map_data("world", region = c("USA", "Canada", "Mexico"))

# Plot the map and the NDVI differences
ggplot() +
  geom_polygon(data = north_america_map, aes(x = long, y = lat, group = group), fill = "lightgray", color = "white") +
  geom_point(data = avg_ndvi_wide_with_coords, aes(x = Long, y = Lat, size = abs(Change_in_max_NDVI), color = color_indicator), alpha = 0.6) +
  scale_size(range = c(1,2)) +
  scale_color_manual(values = c("Positive" = "green", "Negative" = "red", "Neutral" = "grey")) +
  labs(title = "Change in max_NDVI Across Peaks", x = "Longitude", y = "Latitude", size = "NDVI Difference", color = "Change") +
  theme_minimal() +
  coord_fixed(1.3)


#make map graph
#divide data by latitude