#set working directory
setwd("~/Desktop/KU/Projects/GeographicTreelinePatterns")
regs<- read.csv("data/Regressions12Nov24.csv")

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

# create latitude plot
png("figures/figure5.png")
lat_plot <- regs %>%
  ggplot(aes(x = Lat, y = change_in_treeline_elevation)) +
  geom_point(colour = "#a845b9") +
  geom_abline(slope = -14.63, intercept = 611.69, colour = "black", lwd = 0.75)+
  labs(x = "Latitude", y = "Change in Treeline Elevation", title = "Latitude vs Treeline") +
  theme_minimal() + 
  theme(panel.grid = element_blank(), axis.line = 
                           element_line(colour = "black", linewidth = 0.25), 
                         axis.text = element_text(family = "sans"), 
                         axis.ticks = element_line(),
                         plot.title = element_blank(),
                         axis.title = element_text(size = 13, family = "sans"))
#create longitude plot
long_plot <- regs %>%
  ggplot(aes(x = Long, y = change_in_treeline_elevation)) +
  geom_point(colour = "#5b6c65") +
  labs(x = "Longitude", y = element_blank(), title = "Longitude vs Treeline") +
  geom_abline(slope = 14.313, intercept = 1636.478, colour = "black", lwd = 0.75, )+
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = 
          element_line(colour = "black", linewidth = 0.25), 
        axis.text = element_text(family = "sans"), 
        axis.ticks = element_line(),
        plot.title = element_blank(),
        axis.title = element_text(size = 13, family = "sans"))

#combine the plots
lat_plot + long_plot  # Combines the two plots side by side

dev.off()

#distance to coast vs longitude
regs %>%
  ggplot(aes(x = dist_coast, y = Long)) + 
  geom_point(size = 3, alpha = 0.7) + 
  #geom_abline(slope = -14.63, intercept = 611.69, colour = "orange", lwd = 1.5) +
  #geom_abline(slope = 14.313, intercept = 1636.478, colour = "black", lwd = 1.5) +
  xlab("Distance to Coast (m)") + 
  ylab("Longitude") + theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = 
          element_line(colour = "black", linewidth = 0.25), 
        axis.text = element_text(family = "sans"), 
        axis.ticks = element_line(),
        axis.title = element_text(size = 13, family = "sans"))




#malinche
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
