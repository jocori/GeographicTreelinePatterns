#set working directory
setwd("~/Desktop/KU/Projects/GeographicTreelinePatterns")
# Create results directory if it doesn't exist
if (!dir.exists("results")) {
  dir.create("results")
}
#load in packages
library(dplyr)
library(lme4)
library(AICcmodavg)
library(car)
library(corrplot)
library(spdep)
library(sp)
library(spaMM) # spatial mixed models
library(RSpectra) # for fast calculation of extreme eigenvector in spamm models
library(xtable) #for creating latex table
library(broom.mixed) #for adding random effects to latex table
library(texreg)
#read in data
regs<- read.csv("data/regsJan2025.csv")
dist_coast<- read.csv("data/Regressions12Nov24.csv")

# Add only the 'dist_coast' column to 'regs'
regs <- regs %>%
  left_join(dist_coast %>% select(Peak, Direction, dist_coast), by = c("Peak", "Direction"))

#Add columns for scaled longitude and latitude
regs <- regs %>%
  mutate(Lat_scaled = scale(Lat),
         Long_scaled = scale(Long))
#save final dataset
write.csv(regs, "data/regs_final22jan25.csv")

#calculate k nearest neighbors for moran's test (test for spatial autocorrelation)
coords <- as.matrix(cbind(regs$Long, regs$Lat))
neighbors <- dnearneigh(coords, d1 = 0, d2 = 1500, longlat = TRUE)  # Large threshold for connectivity
listw <- nb2listw(neighbors, style = "W")
moran.test(regs$change_in_treeline_NDVI, listw = listw)
summary(regs$change_in_treeline_NDVI)
summary(regs$change_in_treeline_elevation)


#linear mixed models
m1<-lmer(change_in_treeline_NDVI~Lat +Long + (1|Peak_ID), data = regs)
summary(m1)
confint(m1)
AIC(m1)
vif(m1)
#moran.test(resid(m1),listw =listw)
m2<-lmer(change_in_treeline_NDVI~Lat *Long + (1|Peak_ID), data = regs)
summary(m2)
AIC(m2)
vif(m2)
#moran.test(resid(m2),listw =listw)
m3<-lmer(change_in_treeline_NDVI~Lat +Long + Stations_After_Treeline +(1|Peak_ID), data = regs)
summary(m3)
AIC(m3)
vif(m3)
m4<-lmer(change_in_treeline_NDVI~  Stations_After_Treeline +Lat*Long+(1|Peak_ID), data = regs)
summary(m4)
AIC(m4)
vif(m4)
#moran.test(resid(m4),listw =listw)
m5<-lmer(change_in_treeline_NDVI~Lat +Long + Direction +(1|Peak_ID), data = regs)
summary(m5)
AIC(m5)
vif(m5)
#moran.test(resid(m5),listw =listw)
m6<-lmer(change_in_treeline_NDVI~ Direction +Lat*Long+(1|Peak_ID), data = regs)
summary(m6)
AIC(m6)
vif(m6)
#moran.test(resid(m6),listw =listw)
m7<-lmer(change_in_treeline_NDVI~Lat +Long + Stations_After_Treeline +Direction+(1|Peak_ID), data = regs)
summary(m7)
AIC(m7)
vif(m7)
#moran.test(resid(m7),listw =listw)
m8<-lmer(change_in_treeline_NDVI~Stations_After_Treeline +Direction+Lat*Long+(1|Peak_ID), data = regs)
summary(m8)
confint(m8)
AIC(m8)
vif(m8)
#moran.test(resid(m8),listw =listw)
m9<-lmer(change_in_treeline_NDVI~Lat+Long +dist_coast+(1|Peak_ID), data = regs)
summary(m9)
AIC(m9)
vif(m9)
#moran.test(resid(m9),listw =listw)
m10<-lmer(change_in_treeline_NDVI~dist_coast+Lat*Long +(1|Peak_ID), data = regs)
summary(m10)
AIC(m10)
vif(m10)
#moran.test(resid(m10),listw =listw)
m11<-lmer(change_in_treeline_NDVI~dist_coast+Stations_After_Treeline+Lat+Long +(1|Peak_ID), data = regs)
summary(m11)
AIC(m11)
vif(m11)
#moran.test(resid(m11),listw =listw)
m12<-lmer(change_in_treeline_NDVI~dist_coast+Stations_After_Treeline+Lat*Long +(1|Peak_ID), data = regs)
summary(m12)
AIC(m12)
vif(m12)
#moran.test(resid(m12),listw =listw)
m13<-lmer(change_in_treeline_NDVI~dist_coast+Direction+Lat+Long +(1|Peak_ID), data = regs)
summary(m13)
AIC(m13)
vif(m13)
#moran.test(resid(m13),listw =listw)
m14<-lmer(change_in_treeline_NDVI~dist_coast+Direction+Lat*Long +(1|Peak_ID), data = regs)
summary(m14)
AIC(m14)
vif(m14)
#moran.test(resid(m14),listw =listw)
m15<-lmer(change_in_treeline_NDVI~dist_coast+Direction+Stations_After_Treeline+Lat+Long +(1|Peak_ID), data = regs)
summary(m15)
confint(m15)
AIC(m15)
vif(m15)
#moran.test(resid(m15),listw =listw)
m16<-lmer(change_in_treeline_NDVI~dist_coast+Direction+
            Stations_After_Treeline+Lat*Long +(1|Peak_ID), data = regs)
summary(m16)
confint(m16)
AIC(m16)
vif(m16)
#moran.test(resid(m16),listw =listw)

# best linear mixed model is m1
mods<-list(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16)
aictab(cand.set = mods)
confint(m8) #significant interaction between long and lat

# Helper function to format numbers
format_numeric <- function(x) {
  format(x, digits = 4, scientific = TRUE)
}

# Define model terms for lmer models
lmer_terms <- c(
  "Latitude + Longitude",
  "Latitude x Longitude",
  "Latitude + Longitude + # Stations After Treeline",
  "# Stations After Treeline + Latitude x Longitude",
  "Latitude + Longitude + Direction",
  "Direction + Latitude x Longitude",
  "Latitude + Longitude + # Stations After Treeline  + Direction",
  "# Stations After Treeline  + Direction + Latitude x Longitude",
  "Latitude + Longitude + Distance to the Coast (m)",
  "Distance to the Coast (m) + Latitude x Longitude",
  "Distance to the Coast (m) + # Stations After Treeline  + Latitude + Longitude",
  "Distance to the Coast (m) + # Stations After Treeline  + Latitude x Longitude",
  "Distance to the Coast (m) + Direction + Latitude + Longitude",
  "Distance to the Coast (m) + Direction + Latitude x Longitude",
  "Distance to the Coast (m) + Direction + # Stations After Treeline  + Latitude + Longitude",
  "Distance to the Coast (m) + Direction + # Stations After Treeline  + Latitude x Longitude"
)

# Calculate delta AIC and weights
lmer_aic_values <- sapply(mods[1:16], AIC)
#lmer_edf_values <- sapply(mods[1:16], extractAIC)[1,]
delta_aic_lmer <- lmer_aic_values - min(lmer_aic_values)
aic_weights_lmer <- exp(-0.5 * delta_aic_lmer) / sum(exp(-0.5 * delta_aic_lmer))

# Create the AIC table
lmer_aic_table <- data.frame(
  Terms = lmer_terms,
  AIC = format_numeric(lmer_aic_values),
  Delta_AIC = format_numeric(delta_aic_lmer),
  Weight = format_numeric(aic_weights_lmer)
)

# Sort by AIC
lmer_aic_table <- lmer_aic_table[order(as.numeric(lmer_aic_table$AIC)), ]
colnames(lmer_aic_table)[colnames(lmer_aic_table) == "Delta_AIC"] <- "Delta AIC"

# Save to CSV
write.csv(
  lmer_aic_table,
  file = "results/lmer_AIC_ndvi.csv",
  row.names = FALSE
)

## m1 is best model

## spatial mixed models 
#Fit a spatially correlated random effect to check consistency
m1_spamm <- fitme(change_in_treeline_NDVI~Lat +Long+ 
                    Matern(1|Long +Lat), data = regs, method = "REML")
summary(m1_spamm)
AIC(m1_spamm)
moran.test(resid(m1_spamm),listw =listw)

m2_spamm<-fitme(change_in_treeline_NDVI~Lat *Long +
                  Matern(1|Long +Lat), data = regs, method = "REML")
summary(m2_spamm)
AIC(m2_spamm)
moran.test(resid(m2_spamm),listw =listw)

m3_spamm<-fitme(change_in_treeline_NDVI~Lat +Long + 
                  Stations_After_Treeline +Matern(1|Long +Lat),
                data = regs, method = "REML")
summary(m3_spamm)
AIC(m3_spamm)
moran.test(resid(m3_spamm),listw =listw)

m4_spamm<-fitme(change_in_treeline_NDVI~  Stations_After_Treeline +
                  Lat*Long+Matern(1|Long +Lat), 
                data = regs, method = "REML")
summary(m4_spamm)
AIC(m4_spamm)
moran.test(resid(m4_spamm),listw =listw)

m5_spamm<-fitme(change_in_treeline_NDVI~Lat +
                  Long + Direction +Matern(1|Long +Lat), 
                data = regs, method = "REML")
summary(m5_spamm)
AIC(m5_spamm)
moran.test(resid(m5_spamm),listw =listw)

m6_spamm<-fitme(change_in_treeline_NDVI~ Direction +
                  Lat*Long+Matern(1|Long +Lat), 
                data = regs, method = "REML")
summary(m6_spamm)
AIC(m6_spamm)
moran.test(resid(m6_spamm),listw =listw)

m7_spamm<-fitme(change_in_treeline_NDVI~Lat +
                  Long + Stations_After_Treeline +
                  Direction+Matern(1|Long +Lat), data = regs, 
                method = "REML")
summary(m7_spamm)
AIC(m7_spamm)
moran.test(resid(m7_spamm),listw =listw)

m8_spamm<-fitme(change_in_treeline_NDVI~Stations_After_Treeline +
                  Direction+Lat*Long+Matern(1|Long +Lat), 
                data = regs, method = "REML")
summary(m8_spamm)
AIC(m8_spamm)
moran.test(resid(m8_spamm),listw =listw)

m9_spamm<-fitme(change_in_treeline_NDVI~Lat+
                  Long +dist_coast+Matern(1|Long +Lat), 
                data = regs, method = "REML")
summary(m9_spamm)
AIC(m9_spamm)
moran.test(resid(m9_spamm),listw =listw)

m10_spamm<-fitme(change_in_treeline_NDVI~dist_coast+
                   Lat*Long +Matern(1|Long +Lat), 
                 data = regs, method = "REML")
summary(m10_spamm)
AIC(m10_spamm)
moran.test(resid(m10_spamm),listw =listw)

m11_spamm<-fitme(change_in_treeline_NDVI~dist_coast+
                   Stations_After_Treeline+Lat+
                   Long +Matern(1|Long +Lat), 
                 data = regs, method = "REML")
summary(m11_spamm)
AIC(m11_spamm)
moran.test(resid(m11_spamm),listw =listw)

m12_spamm<-fitme(change_in_treeline_NDVI~dist_coast+
                   Stations_After_Treeline+Lat*Long +
                   Matern(1|Long +Lat), data = regs, method = "REML")
summary(m12_spamm)
AIC(m12_spamm)
moran.test(resid(m12_spamm),listw =listw)

m13_spamm<-fitme(change_in_treeline_NDVI~dist_coast+
                   Direction+Lat+Long +
                   Matern(1|Long +Lat), data = regs, method = "REML")
summary(m13_spamm)
AIC(m13_spamm)
moran.test(resid(m13_spamm),listw =listw)

m14_spamm<-fitme(change_in_treeline_NDVI~dist_coast+
                   Direction+Lat*Long +
                   Matern(1|Long +Lat), data = regs, method = "REML")
summary(m14_spamm)
AIC(m14_spamm)
moran.test(resid(m14_spamm),listw =listw)

m15_spamm<-fitme(change_in_treeline_NDVI~dist_coast+
                   Direction+Stations_After_Treeline+
                   Lat+Long +Matern(1|Long +Lat), 
                 data = regs, method = "REML")
summary(m15_spamm)

AIC(m15_spamm)
moran.test(resid(m15_spamm),listw =listw)

m16_spamm<-fitme(change_in_treeline_NDVI~dist_coast+
                   Direction+Stations_After_Treeline+
                   Lat*Long +Matern(1|Long +Lat), 
                 data = regs, method = "REML")
summary(m16_spamm)
AIC(m16_spamm)
moran.test(resid(m16_spamm),listw =listw)
#best spatial mixed model
mods_spamm<-list(m1_spamm,m2_spamm,m3_spamm,m4_spamm,m5_spamm,m6_spamm,m7_spamm,
                 m8_spamm,m9_spamm,m10_spamm,m11_spamm,m12_spamm,m13_spamm,m14_spamm,
                 m15_spamm,m16_spamm)
# Extract AIC for each model
aic_values <- sapply(mods_spamm, extractAIC)

# Extract edf and AIC values from the aic_values matrix
#edf_values <- aic_values["edf", ]
aic_values_only <- aic_values["AIC", ]
spamm_terms <- lmer_terms
delta_AIC_spamm <-aic_values_only - min(aic_values_only) # calculate delta AIC
AIC_weight_spamm <-exp(-0.5 * delta_AIC_spamm) / 
  sum(exp(-0.5 * delta_AIC_spamm)) #calculate AIC weights
# Create a data frame
model_names <- paste0("m", 1:16, "_spamm")
comparison_table_spamm <- data.frame(
  Terms = spamm_terms,
  AIC = format_numeric(aic_values_only),
  Delta_AIC = format_numeric(delta_AIC_spamm),
  Weight = format_numeric(AIC_weight_spamm)
)


# Print the resulting table
print(comparison_table_spamm)

# Sort by AIC
comparison_table_spamm <- comparison_table_spamm[order(as.numeric(comparison_table_spamm$AIC)), ]
colnames(comparison_table_spamm)[colnames(comparison_table_spamm) == "Delta_AIC"] <- "Delta AIC"
colnames(comparison_table_spamm)[colnames(comparison_table_spamm) == "Type"] <- "Model Type"
# Save to CSV
write.csv(
  comparison_table_spamm,
  file = "results/spamm_AIC_ndvi.csv",
  row.names = FALSE
)
#best model is m1_spamm
confint(m2_spamm, parm = names(fixef(m2_spamm))) #significant interaction
confint(m1_spamm, parm = names(fixef(m1_spamm)))


#calculate k nearest neighbors for moran's test (test for spatial autocorrelation)
coords <- as.matrix(cbind(regs$Long, regs$Lat))
neighbors <- dnearneigh(coords, d1 = 0, d2 = 1500, longlat = TRUE)  # Large threshold for connectivity
listw <- nb2listw(neighbors, style = "W")
moran.test(regs$change_in_treeline_NDVI, listw = listw)


#####################################
#######best model overall out of linear mixed models, spatial mixed models, and pcnm
all_mods <- c(mods, mods_spamm)
# Extract AIC for each model
aic_values <- sapply(all_mods, extractAIC)

# Extract edf and AIC values from the aic_values matrix
edf_values <- aic_values[1, ]
aic_values_only <- aic_values[2, ]

# Create a data frame
# Generate model names for the first group
model_names_group1 <- rep("Linear",16)
# Generate model names for the second group
model_names_group2 <- rep("Spatial",16)

# Combine all the model names into one vector
model_names <- c(model_names_group1, model_names_group2)


delta_AIC <-aic_values_only - min(aic_values_only) # calculate delta AIC
AIC_weight <-exp(-0.5 * delta_AIC) / 
  sum(exp(-0.5 * delta_AIC)) #calculate AIC weights
terms<- c(lmer_terms, spamm_terms)
# Create a data frame

comparison_table<- data.frame(
  Type = model_names,
  Terms = terms,
  AIC = format_numeric(aic_values_only),
  DeltaAIC = format_numeric(delta_AIC),
  Weight = format_numeric(AIC_weight),check.names = FALSE
)


# Sort by AIC for better comparison
comparison_table_all <- comparison_table[order(comparison_table$DeltaAIC), ]

# Rename the DeltaAIC column to Delta AIC
colnames(comparison_table_all)[colnames(comparison_table_all) == "DeltaAIC"] <- "Delta AIC"
colnames(comparison_table_all)[colnames(comparison_table_all) == "Type"] <- "Model Type"
# Print the resulting table
print(comparison_table_all)

# Save to CSV
write.csv(
  comparison_table_all,
  file = "results/total_AIC_ndvi.csv",
  row.names = FALSE
)
#best overall is m1_spamm

######best model out of top models in each category
mods_best<-list(m1, m1_spamm)
model_names <- c("m1","m1_spamm")
# Extract AIC for each model
aic_values <- sapply(mods_best, extractAIC)

# Extract edf and AIC values from the aic_values matrix
#edf_values <- aic_values[1, ]
aic_values_only <- aic_values[2, ]
delta_AIC <-aic_values_only - min(aic_values_only) # calculate delta AIC
AIC_weight <-exp(-0.5 * delta_AIC) / 
  sum(exp(-0.5 * delta_AIC)) #calculate AIC weights
m1_terms <- c("Latitude + Longitude")
m1_spamm_terms <- c("Latitude + Longitude ")
model_names_best <- c("Linear", "Spatial")
#m9_pcnm_terms <- c("PCNM1 + PCNM2 + dist_coast")
terms<- c(m1_terms,m1_spamm_terms) # terms for table calrity
comparison_table <- data.frame(
  Type = model_names_best,
  Terms = terms,
  AIC = format_numeric(aic_values_only),
  Delta_AIC = format_numeric(delta_AIC),
  Weight = format_numeric(AIC_weight),
  check.names = FALSE
)
comparison_table_best <- comparison_table[order(comparison_table$Delta_AIC), ]
colnames(comparison_table_best)[colnames(comparison_table_best) == "Delta_AIC"] <- "Delta AIC"
colnames(comparison_table_best)[colnames(comparison_table_best) == "Type"] <- "Model Type"

comparison_table_best
#plot(residuals(m2_spamm)) #residuals are fine

# Save to CSV
write.csv(
  comparison_table_best,
  file = "results/best_AIC_ndvi.csv",
  row.names = FALSE
)

##### make tables for the 3 best models (one of each analysis type: liner, spatial, and pcnm)

# Extract and save outputs for the three best models
#m1_output <- extract_model_output(m1, "m1")

#fixed_effects <- fixef(m1_spamm)
#best_fit_params <- confint(m2_spamm,parm = names(fixed_effects))$confint_best_fit

# Refit the model with these parameters
#refitted_model <- update(m2_spamm, init = best_fit_params)

# Use the refitted model for analysis
#summary(refitted_model)
#confint(refitted_model, parm = names(fixed_effects))
###HAD TO MANUALLY TYPE TABLE FOR THE SPAMM MODEL

#for manually extracting random effects and adding them to the table
#summary(m8)
#summary(m9_pcnm)


##get LaTeX code for all tables

#LMER
# Generate the LaTeX table
latex_table_lmer <- xtable(lmer_aic_table, 
                      label = "tab:lmer_aic",
                      digits = 4, align = c(rep("l",5)), display = c("s","G","G","G","G"))

# Print the LaTeX code
print(latex_table_lmer, 
      include.rownames = FALSE, # Do not include row names
      floating = TRUE,          # Add LaTeX float environment
      table.placement = "t")    # Placement option [t] for top
# Save the LaTeX code to a file
sink("results/lmer_aic_ndvi.tex")
print(latex_table_lmer, 
      include.rownames = FALSE, 
      floating = TRUE, 
      table.placement = "t")
sink()

#SPAMM
# Generate the LaTeX table
latex_table <- xtable(comparison_table_spamm, 
                      label = "tab:spamm_aic",
                      digits = 4, align = c(rep("l",5)), display = c("s","G","G","G","G"))

# Print the LaTeX code
print(latex_table, 
      include.rownames = FALSE, # Do not include row names
      floating = TRUE,          # Add LaTeX float environment
      table.placement = "t")    # Placement option [t] for top
# Save the LaTeX code to a file
sink("results/spamm_aic_ndvi.tex")
print(latex_table, 
      include.rownames = FALSE, 
      floating = TRUE, 
      table.placement = "t")
sink()

# all models
# Generate the LaTeX table
latex_table <- xtable(comparison_table_all, 
                      label = "tab:total_aic",
                      digits = 4, align = c(rep("l",6)), display = c("s","G","G","G","G","G"))

# Print the LaTeX code
print(latex_table, 
      include.rownames = FALSE, # Do not include row names
      floating = TRUE,          # Add LaTeX float environment
      table.placement = "t")    # Placement option [t] for top
# Save the LaTeX code to a file
sink("results/total_aic_ndvi.tex")
print(latex_table, 
      include.rownames = FALSE, # Do not include row names
      floating = TRUE,          # Add LaTeX float environment
      table.placement = "t") 
sink()
#best models only
# Generate the LaTeX table
latex_table <- xtable(comparison_table_best, 
                      label = "tab:best_aic",
                      digits = 4, align = c(rep("l",6)), display = c("s","s","G","G","G","G"))

# Print the LaTeX code
print(latex_table, 
      include.rownames = FALSE, # Do not include row names
      floating = TRUE,          # Add LaTeX float environment
      table.placement = "t")    # Placement option [t] for top
# Save the LaTeX code to a file
sink("results/best_aic_ndvi.tex")
print(latex_table, 
      include.rownames = FALSE, 
      floating = TRUE, 
      table.placement = "t")
sink()


## now I need model output tables for m1, m1_spamm, and m8

###m1

# Extract fixed effects summary from the model
fixed_effects_df <- as.data.frame(coef(summary(m1)))

# Add confidence intervals if needed
conf_intervals <- confint(m1)
conf_intervals<-conf_intervals[3:5,]

# Add confidence intervals to the fixed effects table
fixed_effects_df$CI.Lower <- conf_intervals[, 1]
fixed_effects_df$CI.Upper <- conf_intervals[, 2]
# Rename columns for clarity
colnames(fixed_effects_df) <- c("Estimate", "SE", "T-Value", "Lower 95% CI", "Upper 95% CI")
rownames(fixed_effects_df) <- c("Intercept","Latitude","Longitude")

#extract random effects
random_effects <- VarCorr(m1)
random_effects_summary <- as.data.frame(random_effects)

# Extract random intercept variance and standard deviation
random_variance <- random_effects_summary$vcov[1]
random_sd <- sqrt(random_variance)

# Create a data frame for random effects
random_effects_df <- data.frame(
  Term = c("Random intercept (variance)", "Random intercept (std. dev.)"),
  Estimate = c(random_variance, random_sd),
  SE = NA,  # Random effects typically don't have standard errors
  'TValue' = NA,    # No t-value for random effects
  'LowerCI' = NA,   # No confidence intervals for random effects
  'UpperCI' = NA
)
colnames(random_effects_df) <- c("Term","Estimate", "SE", "T-Value", "Lower 95% CI", "Upper 95% CI")
#combine fixed and random effects into single table
# Add term names to the fixed effects
fixed_effects_df$Term <- rownames(fixed_effects_df)

# Combine fixed and random effects
combined_effects <- rbind(
  fixed_effects_df[, c("Term", "Estimate", "SE", "T-Value", "Lower 95% CI", "Upper 95% CI")],
  random_effects_df
)

# Convert to LaTeX table
latex_table <- xtable(
  combined_effects,
  label = "tab:model_summary",
  digits = 4, align = c(rep("l",7)), display = c("s","G","G","G","G","G","G"))

# Save as a .tex file
sink("results/model_summary_table_NDVI_m1.tex")
print(latex_table, include.rownames = FALSE)
sink()

#### m1_spamm
# Extract fixed effects summary from the model
fixed_effects_df <- as.data.frame(summary(m1_spamm, details = FALSE, verbose = FALSE)$beta_table)

# Add confidence intervals if needed
conf_intervals <- confint(m1_spamm, parm = names(fixef(m1_spamm)))
conf_intervals<-attr(conf_intervals,"table")

# Add confidence intervals to the fixed effects table
fixed_effects_df$CI.Lower <- conf_intervals[, 1]
fixed_effects_df$CI.Upper <- conf_intervals[, 2]

# Rename columns for clarity
colnames(fixed_effects_df) <- c("Estimate", "SE", "T-Value", "Lower 95% CI", "Upper 95% CI")
rownames(fixed_effects_df) <- c("Intercept","Latitude","Longitude")
#extract random effects
random_effects <- VarCorr(m1_spamm)
random_effects_summary <- as.data.frame(random_effects[1,])

# Create a data frame for random effects
random_effects_df <- data.frame(
  Term = c("Random intercept (variance)", "Random intercept (std. dev.)"),
  Estimate = c(random_effects_summary[,3],random_effects_summary[,4]),
  Std.Error = NA,  # Random effects typically don't have standard errors
  t.value = NA,    # No t-value for random effects
  CI.Lower = NA,   # No confidence intervals for random effects
  CI.Upper = NA
)
colnames(random_effects_df) <- c("Term","Estimate", "SE", "T-Value", "Lower 95% CI", "Upper 95% CI")
#combine fixed and random effects into single table
# Add term names to the fixed effects
fixed_effects_df$Term <- rownames(fixed_effects_df)

# Combine fixed and random effects
combined_effects <- rbind(
  fixed_effects_df[, c("Term", "Estimate", "SE", "T-Value", "Lower 95% CI", "Upper 95% CI")],
  random_effects_df
)

# Convert to LaTeX table
latex_table <- xtable(
  combined_effects,
  label = "tab:model_summary_m1spamm",
  digits = 4, align = c(rep("l",7)), display = c("s","G","G","G","G","G","G")
)

# Save as a .tex file
sink("results/model_summary_table_NDVI_m1_spamm.tex")
print(latex_table, include.rownames = FALSE)
sink()


