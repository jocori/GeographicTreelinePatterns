#set working directory
setwd("~/Desktop/KU/Projects/GeographicTreelinePatterns")
regs<- read.csv("data/regs_final22jan25.csv")

#load packages
library(dplyr)
library(ggplot2)
library(tidyverse)
library(patchwork)
# figure 5

#get the slope for lines on the graph
m1<-lm(regs$change_in_treeline_elevation~regs$Lat)
summary(m1)
m2<-lm(regs$change_in_treeline_elevation~regs$Long)
summary(m2)
#########################
library(ggplot2)
library(patchwork)

#Figure 4
# Define a function to format coordinates with degree symbols
format_lat <- function(lat) {
  ifelse(lat >= 0, paste0(lat, "° N"), paste0(abs(lat), "° S"))
}

format_long <- function(long) {
  ifelse(long >= 0, paste0(long, "° E"), paste0(abs(long), "° W"))
}

#define panel labels based on journal specifications
tags<- c("(a)","(b)")

# create latitude plot
png("figures/figure4.png", width = 7.5, height = 5, units = "in", res = 300)
lat_plot <- regs %>%
  ggplot(aes(x = Lat, y = change_in_treeline_elevation)) +
  geom_point(colour = "#a845b9") +
  geom_abline(slope = -14.63, intercept = 611.69, colour = "black", lwd = 0.75) +
  labs(x = "Latitude", y = "Change in Treeline Elevation (m)") +
  scale_x_continuous(limits = c(0, 55),
                     breaks = seq(0, 55, by = 10),
                     labels = format_lat) +
  annotate("text", 
           x = 0,  # Shift slightly inside
           y = 4000,  
           label = "(a)", 
           hjust = 0, 
           vjust = 1, 
           size = 5)+

  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.25),
    axis.text = element_text(family = "sans"),
    axis.ticks = element_line(),
    plot.title = element_blank(),
    axis.title = element_text(size = 13, family = "sans"),
    margin(0,0,l = 5,0)
  )


# Create the Longitude plot
long_plot <- regs %>%
  ggplot(aes(x = Long, y = change_in_treeline_elevation)) +
  geom_point(colour = "#5b6c65") +
  geom_abline(slope = 14.313, intercept = 1636.478, colour = "black", lwd = 0.75) +
  labs(x = "Longitude", y = element_blank()) +
  scale_x_continuous(
    limits = c(-130, -80),
    breaks = seq(-130, -80, by = 10),
    labels = function(x) paste0(abs(x), "°W")  # Convert negatives to absolute values with "°W"
  ) +  # Apply longitude formatting
  annotate("text", x = min(regs$Long), y = 4000,
           label = "(b)", hjust = -0.1, vjust = 0.2, size = 5) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.25),
    axis.text = element_text(family = "sans"),
    axis.ticks = element_line(),
    plot.title = element_blank(),
    axis.title = element_text(size = 13, family = "sans")
  )

# Combine the plots
lat_plot + long_plot # Combines the two plots side by side
dev.off()




#malinche (figure 2)
dat<- read.csv("data/20231208_Averages.csv")

regs<-regs[,-c(20:24)]
second_highest <- regs %>%
  arrange(desc(change_in_treeline_elevation)) %>%  # Sort by change_in_treeline in descending order
  slice(3)  # Select the second row

# Display the corresponding Peak
second_highest$Peak
regs[regs$Peak == "Malinche",]
# Filter data for Cerro Chirippo
malin_data <- dat %>%
  filter(Name == "Malinche") %>%
  select(elevation, Avg8488, Avg1317) %>%
  pivot_longer(cols = c(Avg8488, Avg1317), names_to = "Period", values_to = "NDVI")
png(filename = "figures/figure2b.png", width = 17, units = "cm", height = 11,res = 300)
# Plot with regression lines
ggplot(malin_data, aes(x = elevation, y = NDVI, color = Period)) +
  geom_point(alpha = 0.4) +  # Plot raw data points
  geom_smooth(method = "lm", se = FALSE) +  # Add regression lines
  scale_color_manual(values = c("Avg8488" = "darkred", "Avg1317" = "deepskyblue2"), 
                     labels = c("Avg8488" = "1984-1988", "Avg1317" = "2013-2017")) +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = 
          element_line(colour = "black", linewidth = 0.25), 
        axis.text = element_text(family = "sans", color = "black"), 
        axis.ticks = element_line(color = "black"),
        title = element_blank())
dev.off()


#chir<-dat[dat$Name == "Cerro Chirippo",]
#chir<-chir[,c(3,11,14:16)]
#write.csv(chir,"data/chirippo.csv")

###figure 5
#get the slope for lines on the graph
m1<-lm(regs$change_in_treeline_NDVI~regs$Lat)
summary(m1)
m2<-lm(regs$change_in_treeline_NDVI~regs$Long)
summary(m2)
#########################
library(ggplot2)
library(patchwork)

#Figure 5
# Define a function to format coordinates with degree symbols
format_lat <- function(lat) {
  ifelse(lat >= 0, paste0(lat, "° N"), paste0(abs(lat), "° S"))
}

format_long <- function(long) {
  ifelse(long >= 0, paste0(long, "° E"), paste0(abs(long), "° W"))
}


# create latitude plot
png("figures/figure05.png", width = 7.5, height = 5, units = "in", res = 300)
regs %>%
  ggplot(aes(x = Lat, y = change_in_treeline_NDVI)) +
  geom_point(colour = "#a845b9") +
  geom_abline(slope = -0.0035301, intercept = 0.1308954, colour = "black", lwd = 0.75) +
  labs(x = "Latitude", y = "Change in Treeline NDVI") +
  scale_x_continuous(limits = c(0, 55),
                     breaks = seq(0, 55, by = 10),
                     labels = format_lat) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.25),
    axis.text = element_text(family = "sans"),
    axis.ticks = element_line(),
    plot.title = element_blank(),
    axis.title = element_text(size = 13, family = "sans")
  )
dev.off()
# Create the Longitude plot
#long_plot <- regs %>%
#  ggplot(aes(x = Long, y = change_in_treeline_NDVI)) +
#  geom_point(colour = "#5b6c65") +
#  geom_abline(slope = 0.0039667, intercept = 0.4376782, colour = "black", lwd = 0.75) +
#  labs(x = "Longitude", y = element_blank()) +
#  scale_x_continuous(labels = format_long) +  # Apply longitude formatting
#  annotate("text", x = min(regs$Long), y = max(regs$change_in_treeline_NDVI),
#           label = "(b)", hjust = -0.1, vjust = 0.2, size = 5) +
#  theme_minimal() +
#  theme(
#    panel.grid = element_blank(),
#    axis.line = element_line(colour = "black", linewidth = 0.25),
#    axis.text = element_text(family = "sans"),
#    axis.ticks = element_line(),
#    plot.title = element_blank(),
#    axis.title = element_text(size = 13, family = "sans")
#  )

# Combine the plots
#lat_plot + long_plot # Combines the two plots side by side


##figure 6.....didn't end up using these but scared to delete the code just in case
#read in full dataset 
regs_fig6<-read.csv("data/full_regs_for_figure6.csv")
#fix dataframe
regs_fig6<-regs_fig6 %>%
  mutate(
    treeline_elevation = case_when(
      Year == "Avg1317" ~ present_elevation,      # Keep values from column1 for year 2017
      Year == "Avg8488" ~ past_elevation,      # Use values from column2 for year 1984
      TRUE ~ NA_real_              # Handle any other cases (if applicable)
    )
  )

regs_fig6 %>%
  ggplot(aes(x = Year, y = treeline_elevation)) +
  geom_boxplot()+
  labs(x = "Year", y = "Treeline Elevation (m)") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.25),
    axis.text = element_text(family = "sans"),
    axis.ticks = element_line(),
    plot.title = element_blank(),
    axis.title = element_text(size = 13, family = "sans")
  )
regs_fig6 %>%
  ggplot(aes(x = Year, y = treeline_ndvi)) +
  geom_boxplot()+
  labs(x = "Year", y = "NDVI at 2017 Treeline Elevation") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.25),
    axis.text = element_text(family = "sans"),
    axis.ticks = element_line(),
    plot.title = element_blank(),
    axis.title = element_text(size = 13, family = "sans")
  )
