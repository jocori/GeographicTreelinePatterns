regs_plot %>%
  ggplot(aes(x = Lat, y = change_in_treeline_elevation, color = Long)) + 
  geom_point(size = 3, alpha = 0.7) + 
  scale_color_gradient(name = "Longitude", low = "blue", high = "red") + 
  xlab("Latitude") + 
  ylab("Change in Treeline Elevation (2017 - 1984)") + theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = 
          element_line(colour = "black", linewidth = 0.25), 
        axis.text = element_text(family = "sans"), 
        axis.ticks = element_line(),
        axis.title = element_text(size = 13, family = "sans"))

#popocatepetl
dat<- read.csv("data/20231208_Averages.csv")
#popo <- dat[dat$Name == "Volcan Popocatepetl",]

#popo_data<-popo %>%
 # summarise(Period = c("Avg8488", "Avg1317"),
  #          Mean_NDVI = c(mean(Avg8488, na.rm = TRUE), mean(Avg1317, na.rm = TRUE)))


###next attempt
library(tidyverse)
# Reshape the data for raw NDVI points
popo_raw_data <- dat %>%
  filter(Name == "Volcan Popocatepetl") %>%
  select(Avg8488, Avg1317) %>%
  pivot_longer(cols = c(Avg8488, Avg1317), names_to = "Period", values_to = "NDVI")

# Prepare the summary data for the line plot
popo_summary <- popo_raw_data %>%
  group_by(Period) %>%
  summarize(Mean_NDVI = mean(NDVI, na.rm = TRUE), .groups = "drop")

# Plot with raw data points and summary line
ggplot() +
  # Raw data points (transparent)
  geom_point(data = popo_raw_data, aes(x = Period, y = NDVI), alpha = 0.3, color = "gray", size = 2) +
  # Line and point plot for mean NDVI
  geom_line(data = popo_summary, aes(x = Period, y = Mean_NDVI, group = 1), size = 1, color = "maroon") +
  geom_point(data = popo_summary, aes(x = Period, y = Mean_NDVI), size = 3, color = "darkred") +
  ylab("Mean NDVI") +
  xlab("Period") +
  ggtitle("Change in Greenness (NDVI) for Volcán Popocatepetl") +
  scale_x_discrete(
    limits = c("Avg8488", "Avg1317"),  # Specifies the order
    labels = c("Avg8488" = "1984-1988", "Avg1317" = "2013-2017")  # Custom labels
  ) + theme_minimal() + 
  theme(panel.grid = element_blank(), axis.line = 
              element_line(colour = "black", linewidth = 0.25), 
            axis.text = element_text(family = "sans"), 
            axis.ticks = element_line(),
            axis.title = element_text(size = 13, family = "sans"))

