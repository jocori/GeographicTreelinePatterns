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
library(vegan) # PCNM
library(RSpectra) # for fast calculation of extreme eigenvector in spamm models
library(MuMIn) #for calculating R^2 values
library(xtable) #for generating LaTeX tables
#read in data
regs<- read.csv("data/regs_final22jan25.csv")
#check for correlation among predictor variables
numeric_data <- regs[, sapply(regs, is.numeric)]
correlation_matrix <- cor(numeric_data)

#Add columns for scaled longitude and latitude
regs <- regs %>%
  mutate(Lat_scaled = scale(Lat),
         Long_scaled = scale(Long))

#calculate k nearest neighbors for moran's test (test for spatial autocorrelation)
coords <- as.matrix(cbind(regs$Long, regs$Lat))
neighbors <- dnearneigh(coords, d1 = 0, d2 = 1500, longlat = TRUE)  # Large threshold for connectivity
listw <- nb2listw(neighbors, style = "W")
moran.test(regs$change_in_treeline_elevation, listw = listw)

##Plot the correlation matrix as a heatmap
corrplot(correlation_matrix, method = "color", 
         tl.col = "black", tl.srt = 45,
         type = "upper") # autocorrelation between lat and long is concerning

#linear mixed models
m1<-lmer(change_in_treeline_elevation~Lat +Long + (1|Peak_ID), data = regs)
summary(m1)
AIC(m1)
vif(m1)
moran.test(resid(m1),listw =listw)
m2<-lmer(change_in_treeline_elevation~Lat *Long + (1|Peak_ID), data = regs)
summary(m2)
AIC(m2)
vif(m2)
moran.test(resid(m2),listw =listw)
m3<-lmer(change_in_treeline_elevation~Lat +Long + Stations_After_Treeline +(1|Peak_ID), data = regs)
summary(m3)
AIC(m3)
vif(m3)
moran.test(resid(m3),listw =listw)
m4<-lmer(change_in_treeline_elevation~  Stations_After_Treeline +Lat*Long+(1|Peak_ID), data = regs)
summary(m4)
AIC(m4)
vif(m4)
moran.test(resid(m4),listw =listw)
m5<-lmer(change_in_treeline_elevation~Lat +Long + Direction +(1|Peak_ID), data = regs)
summary(m5)
AIC(m5)
vif(m5)
moran.test(resid(m5),listw =listw)
m6<-lmer(change_in_treeline_elevation~ Direction +Lat*Long+(1|Peak_ID), data = regs)
summary(m6)
AIC(m6)
vif(m6)
moran.test(resid(m6),listw =listw)
m7<-lmer(change_in_treeline_elevation~Lat +Long + Stations_After_Treeline +Direction+(1|Peak_ID), data = regs)
summary(m7)
AIC(m7)
vif(m7)
moran.test(resid(m7),listw =listw)
m8<-lmer(change_in_treeline_elevation~Stations_After_Treeline +Direction+Lat*Long+(1|Peak_ID), data = regs)
summary(m8)
AIC(m8)
vif(m8)
moran.test(resid(m8),listw =listw)
m9<-lmer(change_in_treeline_elevation~Lat+Long +dist_coast+(1|Peak_ID), data = regs)
summary(m9)
AIC(m9)
vif(m9)
moran.test(resid(m9),listw =listw)
m10<-lmer(change_in_treeline_elevation~dist_coast+Lat*Long +(1|Peak_ID), data = regs)
summary(m10)
AIC(m10)
vif(m10)
moran.test(resid(m10),listw =listw)
m11<-lmer(change_in_treeline_elevation~dist_coast+Stations_After_Treeline+Lat+Long +(1|Peak_ID), data = regs)
summary(m11)
AIC(m11)
vif(m11)
moran.test(resid(m11),listw =listw)
m12<-lmer(change_in_treeline_elevation~dist_coast+Stations_After_Treeline+Lat*Long +(1|Peak_ID), data = regs)
summary(m12)
AIC(m12)
vif(m12)
moran.test(resid(m12),listw =listw)
m13<-lmer(change_in_treeline_elevation~dist_coast+Direction+Lat+Long +(1|Peak_ID), data = regs)
summary(m13)
AIC(m13)
vif(m13)
moran.test(resid(m13),listw =listw)
m14<-lmer(change_in_treeline_elevation~dist_coast+Direction+Lat*Long +(1|Peak_ID), data = regs)
summary(m14)
AIC(m14)
vif(m14)
moran.test(resid(m14),listw =listw)
m15<-lmer(change_in_treeline_elevation~dist_coast+Direction+Stations_After_Treeline+Lat+Long +(1|Peak_ID), data = regs)
summary(m15)
confint(m15)
AIC(m15)
vif(m15)
moran.test(resid(m15),listw =listw)
m16<-lmer(change_in_treeline_elevation~dist_coast+Direction+
            Stations_After_Treeline+Lat*Long +(1|Peak_ID), data = regs)
summary(m16)
confint(m16)
AIC(m16)
vif(m16)
moran.test(resid(m16),listw =listw)
m17 <- lmer(change_in_treeline_elevation~Lat + (1|Peak_ID), data = regs)
m18 <- lmer(change_in_treeline_elevation~ Long + (1|Peak_ID), data = regs)

# best linear mixed model is m8
mods<-list(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18)
aictab(cand.set = mods, modnames = names(mods))
confint(m8) #significant interaction between long and lat

#calculate marginal and conditional r-squared values
r.squaredGLMM(m8)

# Helper function to format numbers
format_numeric <- function(x) {
  format(x, digits = 4, scientific = TRUE)
}

# Define model terms for lmer models
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
  "Distance to the Coast (m) + Direction + # Stations After Treeline  + Latitude x Longitude",
  "Latitude",
  "Longitude"
)

# Calculate delta AIC and weights
lmer_aic_values <- sapply(mods[1:18], AIC)
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
lmer_aic_table <- lmer_aic_table[order(as.numeric(lmer_aic_table$Delta_AIC)), ]
colnames(lmer_aic_table)[colnames(lmer_aic_table) == "Delta_AIC"] <- "Delta AIC"
colnames(lmer_aic_table)[colnames(lmer_aic_table) == "Terms"] <- "Model Terms"
# Save to CSV
write.csv(
  lmer_aic_table,
  file = "results/lmer_AIC_elevation.csv",
  row.names = FALSE
)


## spatial mixed models 
#Fit a spatially correlated random effect to check consistency
m1_spamm <- fitme(change_in_treeline_elevation~Lat +Long+ 
                   Matern(1|Long +Lat), data = regs, method = "REML")
summary(m1_spamm)
AIC(m1_spamm)
moran.test(resid(m1_spamm),listw =listw)

m2_spamm<-fitme(change_in_treeline_elevation~Lat *Long +
                  Matern(1|Long +Lat), data = regs, method = "REML")
summary(m2_spamm)
AIC(m2_spamm)
moran.test(resid(m2_spamm),listw =listw)

m3_spamm<-fitme(change_in_treeline_elevation~Lat +Long + 
                  Stations_After_Treeline +Matern(1|Long +Lat),
                data = regs, method = "REML")
summary(m3_spamm)
AIC(m3_spamm)
moran.test(resid(m3_spamm),listw =listw)

m4_spamm<-fitme(change_in_treeline_elevation~  Stations_After_Treeline +
                  Lat*Long+Matern(1|Long +Lat), 
                data = regs, method = "REML")
summary(m4_spamm)
AIC(m4_spamm)
moran.test(resid(m4_spamm),listw =listw)

m5_spamm<-fitme(change_in_treeline_elevation~Lat +
                  Long + Direction +Matern(1|Long +Lat), 
                data = regs, method = "REML")
summary(m5_spamm)
AIC(m5_spamm)
moran.test(resid(m5_spamm),listw =listw)

m6_spamm<-fitme(change_in_treeline_elevation~ Direction +
                  Lat*Long+Matern(1|Long +Lat), 
                data = regs, method = "REML")
summary(m6_spamm)
AIC(m6_spamm)
moran.test(resid(m6_spamm),listw =listw)

m7_spamm<-fitme(change_in_treeline_elevation~Lat +
                  Long + Stations_After_Treeline +
                  Direction+Matern(1|Long +Lat), data = regs, 
                method = "REML")
summary(m7_spamm)
AIC(m7_spamm)
moran.test(resid(m7_spamm),listw =listw)

m8_spamm<-fitme(change_in_treeline_elevation~Stations_After_Treeline +
                  Direction+Lat*Long+Matern(1|Long +Lat), 
                data = regs, method = "REML")
summary(m8_spamm)
AIC(m8_spamm)
moran.test(resid(m8_spamm),listw =listw)

m9_spamm<-fitme(change_in_treeline_elevation~Lat+
                  Long +dist_coast+Matern(1|Long +Lat), 
                data = regs, method = "REML")
summary(m9_spamm)
AIC(m9_spamm)
moran.test(resid(m9_spamm),listw =listw)

m10_spamm<-fitme(change_in_treeline_elevation~dist_coast+
                   Lat*Long +Matern(1|Long +Lat), 
                 data = regs, method = "REML")
summary(m10_spamm)
AIC(m10_spamm)
moran.test(resid(m10_spamm),listw =listw)

m11_spamm<-fitme(change_in_treeline_elevation~dist_coast+
                   Stations_After_Treeline+Lat+
                   Long +Matern(1|Long +Lat), 
                 data = regs, method = "REML")
summary(m11_spamm)
AIC(m11_spamm)
moran.test(resid(m11_spamm),listw =listw)

m12_spamm<-fitme(change_in_treeline_elevation~dist_coast+
                   Stations_After_Treeline+Lat*Long +
                   Matern(1|Long +Lat), data = regs, method = "REML")
summary(m12_spamm)
AIC(m12_spamm)
moran.test(resid(m12_spamm),listw =listw)

m13_spamm<-fitme(change_in_treeline_elevation~dist_coast+
                   Direction+Lat+Long +
                   Matern(1|Long +Lat), data = regs, method = "REML")
summary(m13_spamm)
AIC(m13_spamm)
moran.test(resid(m13_spamm),listw =listw)

m14_spamm<-fitme(change_in_treeline_elevation~dist_coast+
                   Direction+Lat*Long +
                   Matern(1|Long +Lat), data = regs, method = "REML")
summary(m14_spamm)
AIC(m14_spamm)
moran.test(resid(m14_spamm),listw =listw)

m15_spamm<-fitme(change_in_treeline_elevation~dist_coast+
                   Direction+Stations_After_Treeline+
                   Lat+Long +Matern(1|Long +Lat), 
                 data = regs, method = "REML")
summary(m15_spamm)
AIC(m15_spamm)
moran.test(resid(m15_spamm),listw =listw)

m16_spamm<-fitme(change_in_treeline_elevation~dist_coast+
                   Direction+Stations_After_Treeline+
                   Lat*Long +Matern(1|Long +Lat), 
                 data = regs, method = "REML")
summary(m16_spamm)
AIC(m16_spamm)
moran.test(resid(m16_spamm),listw =listw)
m17_spamm<-fitme(change_in_treeline_elevation~Lat+Matern(1|Long +Lat), 
                 data = regs, method = "REML")
summary(m17_spamm)

#calculate pseudo R2 
pseudoR2(m17_spamm, nullform = ~1)

m18_spamm<-fitme(change_in_treeline_elevation~Long +Matern(1|Long +Lat), 
                 data = regs, method = "REML")
#best spatial mixed model
mods_spamm<-list(m1_spamm,m2_spamm,m3_spamm,m4_spamm,m5_spamm,m6_spamm,m7_spamm,
           m8_spamm,m9_spamm,m10_spamm,m11_spamm,m12_spamm,m13_spamm,m14_spamm,
           m15_spamm,m16_spamm,m17_spamm,m18_spamm)
# Extract AIC for each model
aic_values_spamm <- sapply(mods_spamm, extractAIC)

# Extract  and AIC values from the aic_values matrix

aic_values_only_spamm <- aic_values_spamm["AIC", ]
spamm_terms <- lmer_terms
delta_AIC_spamm <-aic_values_only_spamm - min(aic_values_only_spamm) # calculate delta AIC
AIC_weight_spamm <-exp(-0.5 * delta_AIC_spamm) / 
  sum(exp(-0.5 * delta_AIC_spamm)) #calculate AIC weights
# Create a data frame

comparison_table_spamm <- data.frame(
  Terms = spamm_terms,
  AIC = format_numeric(aic_values_only_spamm),
  Delta_AIC = format_numeric(delta_AIC_spamm),
  Weight = format_numeric(AIC_weight_spamm)
)


# Print the resulting table
print(comparison_table_spamm)

# Sort by AIC
comparison_table_spamm <- comparison_table_spamm[order(as.numeric(comparison_table_spamm$Delta_AIC)), ]
colnames(comparison_table_spamm)[colnames(comparison_table_spamm) == "Delta_AIC"] <- "Delta AIC"
colnames(comparison_table_spamm)[colnames(comparison_table_spamm) == "Type"] <- "Model Type"
# Save to CSV
write.csv(
  comparison_table_spamm,
  file = "results/spamm_AIC.csv",
  row.names = FALSE
)
#best model is m17_spamm
#confint(m17_spamm, parm = names(fixef(m17_spamm))) #significant interaction

latex_table_lmer <- xtable(lmer_aic_table, 
                           label = "tab:spamm_aic",
                           digits = 4, align = c(rep("l",5)), display = c("s","G","G","G","G"))

# Print the LaTeX code
print(latex_table_lmer, 
      include.rownames = FALSE, # Do not include row names
      floating = TRUE,          # Add LaTeX float environment
      table.placement = "t")    # Placement option [t] for top
# Save the LaTeX code to a file
sink("results/lmer_aic_elevation.tex")
print(latex_table_lmer, 
      include.rownames = FALSE, 
      floating = TRUE, 
      table.placement = "t")
sink()

#SPAMM aic
# Generate the LaTeX table
latex_table_spamm <- xtable(comparison_table_spamm, 
                            label = "tab:spamm_aic",
                            digits = 4, align = c(rep("l",5)), display = c("s","G","G","G","G"))

# Print the LaTeX code
print(latex_table_spamm, 
      include.rownames = FALSE, # Do not include row names
      floating = TRUE,          # Add LaTeX float environment
      table.placement = "t")    # Placement option [t] for top
# Save the LaTeX code to a file
sink("results/spamm_aic_elevation.tex")
print(latex_table_spamm, 
      include.rownames = FALSE, 
      floating = TRUE, 
      table.placement = "t")
sink()

#######best model overall out of linear mixed models, spatial mixed models
all_mods <- c(mods, mods_spamm)
# Extract AIC for each model
aic_values_all <- sapply(all_mods, extractAIC)

# Extract edf and AIC values from the aic_values matrix
edf_values <- aic_values_all[1, ]
aic_values_only_all <- aic_values_all[2, ]

# Create a data frame
# Generate model names for the first group
model_names_group1 <- rep("Linear",18)
# Generate model names for the second group
model_names_group2 <- rep("Spatial",18)

# Combine all the model names into one vector
model_names_all <- c(model_names_group1, model_names_group2)


delta_AIC_all <-aic_values_only_all - min(aic_values_only_all) # calculate delta AIC
AIC_weight_all <-exp(-0.5 * delta_AIC_all) / 
  sum(exp(-0.5 * delta_AIC_all)) #calculate AIC weights
terms_all<- c(lmer_terms, spamm_terms)
# Create a data frame

comparison_table_all<- data.frame(
  Type = model_names_all,
  Terms = terms_all,
  AIC = format_numeric(aic_values_only_all),
  DeltaAIC = format_numeric(delta_AIC_all),
  Weight = format_numeric(AIC_weight_all),check.names = FALSE
)


# Sort by AIC for better comparison
comparison_table_all <- comparison_table_all[order(comparison_table_all$DeltaAIC), ]

# Rename the DeltaAIC column to Delta AIC
colnames(comparison_table_all)[colnames(comparison_table_all) == "DeltaAIC"] <- "Delta AIC"
colnames(comparison_table_all)[colnames(comparison_table_all) == "Type"] <- "Model Type"
# Print the resulting table
print(comparison_table_all)


# Save to CSV
write.csv(
  comparison_table_all,
  file = "results/total_AIC.csv",
  row.names = FALSE
)
latex_table_all <- xtable(comparison_table_all, 
                          label = "tab:total_aic_elevation",
                          digits = 4, align = c(rep("l",6)), display = c("s","G","G","G","G","G"))

# Print the LaTeX code
print(latex_table_all, 
      include.rownames = FALSE, # Do not include row names
      floating = TRUE,          # Add LaTeX float environment
      table.placement = "t")    # Placement option [t] for top
# Save the LaTeX code to a file
sink("results/total_aic_elevation.tex")
print(latex_table_all, 
      include.rownames = FALSE, # Do not include row names
      floating = TRUE,          # Add LaTeX float environment
      table.placement = "t") 
sink()
#best is m17_spamm
######best model out of top models in each category
mods_best<-list(m8, m17_spamm)
model_names_best <- c("Linear","Spatial")
# Extract AIC for each model
aic_values_best <- sapply(mods_best, extractAIC)

# Extract edf and AIC values from the aic_values matrix
edf_values <- aic_values_best[1, ]
aic_values_only_best <- aic_values_best[2, ]
delta_AIC_best <-aic_values_only_best - min(aic_values_only_best) # calculate delta AIC
AIC_weight_best <-exp(-0.5 * delta_AIC_best) / 
  sum(exp(-0.5 * delta_AIC_best)) #calculate AIC weights
m8_terms <- c("# Stations After Treeline + Direction + Latitude * Longitude")
m17_spamm_terms <- c("Latitude")
#m9_pcnm_terms <- c("PCNM1 + PCNM2 + dist_coast")
terms_best<- c(m8_terms,m17_spamm_terms) # terms for table calrity

comparison_table_best <- data.frame(
  Type = model_names_best,
  Terms = terms_best,
  AIC = format_numeric(aic_values_only_best),
  DeltaAIC = format_numeric(delta_AIC_best),
  Weight = format_numeric(AIC_weight_best)
)
comparison_table_best <- comparison_table_best[order(comparison_table_best$AIC), ]
colnames(comparison_table_best)[colnames(comparison_table_best) == "DeltaAIC"] <- "Delta AIC"
colnames(comparison_table_best)[colnames(comparison_table_best) == "Type"] <- "Model Type"
comparison_table_best


# Save to CSV
write.csv(
  comparison_table_best,
  file = "results/best_AIC_elevation.csv",
  row.names = FALSE
)
latex_table_best <- xtable(comparison_table_best, 
                           label = "tab:best_aic",
                           digits = 4, align = c(rep("l",6)), display = c("s","s","G","G","G","G"))

# Print the LaTeX code
print(latex_table_best, 
      include.rownames = FALSE, # Do not include row names
      floating = TRUE,          # Add LaTeX float environment
      table.placement = "t")    # Placement option [t] for top
# Save the LaTeX code to a file
sink("results/best_aic_elevation.tex")
print(latex_table_best, 
      include.rownames = FALSE, 
      floating = TRUE, 
      table.placement = "t")
sink()

### Model Summary Table for m8

# Extract fixed effects
m8.fixed_effects_df <- as.data.frame(coef(summary(m8)))

# Add confidence intervals
m8.conf_intervals <- confint(m8)
m8.conf_intervals <- m8.conf_intervals[3:14, ]  # Remove intercept + random effects CIs

# Add CIs to fixed effects
m8.fixed_effects_df$`Lower 95% CI` <- m8.conf_intervals[, 1]
m8.fixed_effects_df$`Upper 95% CI` <- m8.conf_intervals[, 2]

# Rename columns
colnames(m8.fixed_effects_df)[1:3] <- c("Estimate", "SE", "T-Value")

# Add term names (this order should match model output)
m8.fixed_effects_df$Term <- c(
  "Intercept",
  "# Stations After Treeline",
  "Direction (North)",
  "Direction (Northeast)",
  "Direction (Northwest)",
  "Direction (South)",
  "Direction (Southeast)",
  "Direction (Southwest)",
  "Direction (West)",
  "Latitude",
  "Longitude",
  "Latitude × Longitude"
)

# Reorder columns
m8.fixed_effects_df <- m8.fixed_effects_df[, c("Term", "Estimate", "SE", "T-Value", "Lower 95% CI", "Upper 95% CI")]

# --- RANDOM EFFECTS ---

# Extract random intercept variance
m8.random_effects <- VarCorr(m8)
m8.random_intercept_var <- attr(m8.random_effects$Peak_ID, "stddev")^2
m8.random_intercept_sd  <- sqrt(m8.random_intercept_var)

m8.random_effects_df <- data.frame(
  Term = c("Random intercept (variance)", "Random intercept (std. dev.)"),
  Estimate = c(m8.random_intercept_var, m8.random_intercept_sd),
  SE = NA,
  `T-Value` = NA,
  `Lower 95% CI` = NA,
  `Upper 95% CI` = NA
)

# --- RESIDUAL VARIANCE ---

m8.residual_var <- attr(m8.random_effects, "sc")^2
m8.residual_sd  <- sqrt(m8.residual_var)

m8.residual_df <- data.frame(
  Term = c("Residual (variance)", "Residual (std. dev.)"),
  Estimate = c(m8.residual_var, m8.residual_sd),
  SE = NA,
  `T-Value` = NA,
  `Lower 95% CI` = NA,
  `Upper 95% CI` = NA
)

# --- COMBINE ALL TOGETHER ---

# ensure column names match
colnames(m8.random_effects_df) <- colnames(m8.fixed_effects_df)
colnames(m8.residual_df) <- colnames(m8.fixed_effects_df)

m8.combined_effects <- rbind(
  m8.fixed_effects_df,
  m8.random_effects_df,
  m8.residual_df
)

# Create LaTeX table
m8.latex_table <- xtable(
  m8.combined_effects,
  label = "tab:model_summary_m8",
  digits = 4,
  align = c("l", rep("l", 6)),
  display = c("s", "s", "g", "g", "g", "g", "g")
)

# Export to .tex file
sink("results/model_summary_table_elevation_m8.tex")
print(m8.latex_table, include.rownames = FALSE)
sink()

# ---- Model Summary Table for m17_spamm ----

# Extract fixed effects summary
m17.fixed_df <- as.data.frame(summary(m17_spamm, details = FALSE)$beta_table)

# Add confidence intervals
m17.confint <- confint(m17_spamm, parm = names(fixef(m17_spamm)))
m17.confint <- attr(m17.confint, "table")

# Add confidence intervals to fixed effects table
m17.fixed_df$`Lower 95% CI` <- m17.confint[, 1]
m17.fixed_df$`Upper 95% CI` <- m17.confint[, 2]

# Rename columns
colnames(m17.fixed_df)[1:3] <- c("Estimate", "SE", "T-Value")

# Add term names
m17.fixed_df$Term <- c("Intercept","Latitude")

# Reorder columns
m17.fixed_df <- m17.fixed_df[, c("Term", "Estimate", "SE", "T-Value", "Lower 95% CI", "Upper 95% CI")]

# ---- Extract Variance Components ----

# Random effect variance (lambda) and SD
m17.re.var <- m17_spamm$lambda[[1]]
m17.re.sd  <- sqrt(m17.re.var)

# Residual variance (phi) and SD
m17.resid.var <- m17_spamm$phi
m17.resid.sd  <- sqrt(m17.resid.var)

# ---- Create Random and Residual Effects Data Frames ----

# Random intercept effects
m17.re_df <- data.frame(
  Term = c("Random intercept (variance)", "Random intercept (std. dev.)"),
  Estimate = c(m17.re.var, m17.re.sd),
  SE = NA,
  `T-Value` = NA,
  `Lower 95% CI` = NA,
  `Upper 95% CI` = NA
)

# Residual variance
m17.resid_df <- data.frame(
  Term = c("Residual (variance)", "Residual (std. dev.)"),
  Estimate = c(m17.resid.var, m17.resid.sd),
  SE = NA,
  `T-Value` = NA,
  `Lower 95% CI` = NA,
  `Upper 95% CI` = NA
)

# ---- Standardize column names for binding ----
colnames(m17.re_df) <- colnames(m17.fixed_df)
colnames(m17.resid_df) <- colnames(m17.fixed_df)

# ---- Combine all into final summary table ----
m17.combined <- rbind(
  m17.fixed_df,
  m17.re_df,
  m17.resid_df
)



# ---- Convert to LaTeX table ----

m17.latex_table <- xtable(
  m17.combined,
  label = "tab:model_summary_m17spamm",
  digits = 4,
  align = c("l", rep("l", 6)),
  display = c("s", "s", "g", "g", "g", "g", "g")
)

# ---- Save as .tex file ----
sink("results/model_summary_table_elevation_m17_spamm.tex")
print(m17.latex_table, include.rownames = FALSE)
sink()

#calculate k nearest neighbors for moran's test (test for spatial autocorrelation)
coords <- as.matrix(cbind(regs$Long, regs$Lat))
neighbors <- dnearneigh(coords, d1 = 0, d2 = 1500, longlat = TRUE)  # Large threshold for connectivity
listw <- nb2listw(neighbors, style = "W")
moran.test(regs$change_in_treeline_elevation, listw = listw)

## Eigenvector Filtering or PCNM
# explore Principal Coordinates of Neighbor Matrices (PCNM) 
# to capture latent spatial patterns
#pcnm_res <- pcnm(dist(coords))
#significant_axes <- pcnm_res$vectors[, 1:5] 
#regs <- cbind(regs, significant_axes)

#m1_pcnm<-lmer(change_in_treeline_elevation~PCNM1 + PCNM2+ (1|Peak_ID),
#              data = regs)
#summary(m1_pcnm)
#AIC(m1_pcnm)
#vif(m1_pcnm)
#moran.test(resid(m1),listw =listw)

#m2_pcnm<-lmer(change_in_treeline_elevation~PCNM1*PCNM2+ (1|Peak_ID), data = regs)
#summary(m2_pcnm)
#AIC(m2_pcnm)
#vif(m2_pcnm)
#moran.test(resid(m2_pcnm),listw =listw)

#m3_pcnm<-lmer(change_in_treeline_elevation~Stations_After_Treeline + PCNM1 + PCNM2+(1|Peak_ID), data = regs)
#summary(m3_pcnm)
#AIC(m3_pcnm)
#vif(m3_pcnm)
#moran.test(resid(m3_pcnm),listw =listw)

#m4_pcnm<-lmer(change_in_treeline_elevation~  Stations_After_Treeline + PCNM1*PCNM2+(1|Peak_ID), data = regs)
#summary(m4_pcnm)
#AIC(m4_pcnm)
#vif(m4_pcnm)
#moran.test(resid(m4_pcnm),listw =listw)

#m5_pcnm<-lmer(change_in_treeline_elevation~ PCNM1 + PCNM2+ Direction +(1|Peak_ID), data = regs)
#summary(m5_pcnm)
#AIC(m5_pcnm)
#vif(m5_pcnm)
#moran.test(resid(m5_pcnm),listw =listw)

#m6_pcnm<-lmer(change_in_treeline_elevation~ Direction + PCNM1*PCNM2+(1|Peak_ID), data = regs)
#summary(m6_pcnm)
#AIC(m6_pcnm)
#vif(m6_pcnm)
#moran.test(resid(m6_pcnm),listw =listw)

#m7_pcnm<-lmer(change_in_treeline_elevation~  PCNM1 + PCNM2 + Stations_After_Treeline +Direction+(1|Peak_ID), data = regs)
#summary(m7_pcnm)
#AIC(m7_pcnm)
#vif(m7_pcnm)
#moran.test(resid(m7_pcnm),listw =listw)

#m8_pcnm<-lmer(change_in_treeline_elevation~Stations_After_Treeline +Direction+ PCNM1*PCNM2+(1|Peak_ID), data = regs)
#summary(m8_pcnm)
#AIC(m8_pcnm)
#vif(m8_pcnm)
#moran.test(resid(m8_pcnm),listw =listw)

#m9_pcnm<-lmer(change_in_treeline_elevation~PCNM1 + PCNM2+dist_coast+(1|Peak_ID), data = regs)
#summary(m9_pcnm)
#AIC(m9_pcnm)
#vif(m9_pcnm)
#moran.test(resid(m9_pcnm),listw =listw)

#m10_pcnm<-lmer(change_in_treeline_elevation~dist_coast+ PCNM1*PCNM2 +(1|Peak_ID), data = regs)
#summary(m10_pcnm)
#AIC(m10_pcnm)
#vif(m10_pcnm)
#moran.test(resid(m10_pcnm),listw =listw)

#m11_pcnm<-lmer(change_in_treeline_elevation~dist_coast+ PCNM1 + PCNM2+Stations_After_Treeline+(1|Peak_ID), data = regs)
#summary(m11_pcnm)
#AIC(m11_pcnm)
#vif(m11_pcnm)
#moran.test(resid(m11_pcnm),listw =listw)

#m12_pcnm<-lmer(change_in_treeline_elevation~dist_coast+Stations_After_Treeline+ PCNM1*PCNM2+(1|Peak_ID), data = regs)
#summary(m12_pcnm)
#AIC(m12_pcnm)
#vif(m12_pcnm)
#moran.test(resid(m12_pcnm),listw =listw)

#m13_pcnm<-lmer(change_in_treeline_elevation~dist_coast+Direction+ PCNM1 + PCNM2+ +(1|Peak_ID), data = regs)
#summary(m13_pcnm)
#AIC(m13_pcnm)
#vif(m13_pcnm)
#moran.test(resid(m13_pcnm),listw =listw)

#m14_pcnm<-lmer(change_in_treeline_elevation~dist_coast+Direction+ PCNM1*PCNM2 +(1|Peak_ID), data = regs)
#summary(m14_pcnm)
#AIC(m14_pcnm)
#vif(m14_pcnm)
#moran.test(resid(m14_pcnm),listw =listw)

#m15_pcnm<-lmer(change_in_treeline_elevation~dist_coast+Direction+ PCNM1 + PCNM2+Stations_After_Treeline+(1|Peak_ID), data = regs)
#summary(m15_pcnm)
#AIC(m15_pcnm)
#vif(m15_pcnm)
#moran.test(resid(m15_pcnm),listw =listw)

#m16_pcnm<-lmer(change_in_treeline_elevation~dist_coast+Direction+Stations_After_Treeline+ PCNM1*PCNM2 +(1|Peak_ID), data = regs)
#summary(m16_pcnm)
#AIC(m16_pcnm)
#vif(m16_pcnm)
#moran.test(resid(m16_pcnm),listw =listw)

#best PCNM model
#mods_pcnm<-list(m1_pcnm,m2_pcnm,m3_pcnm,m4_pcnm,m5_pcnm,
#           m6_pcnm,m7_pcnm,m8_pcnm,m9_pcnm,m10_pcnm,
#           m11_pcnm,m12_pcnm,m13_pcnm,m14_pcnm,m15_pcnm,m16_pcnm)

#aictab(cand.set = mods_pcnm)
#model 8 is best
#confint(m8_pcnm) #only PCNM 1 is significant
# Define terms for pcnm models
#pcnm_terms <- c(
#  "PCNM1 + PCNM2",
#  "PCNM1 * PCNM2",
#  "Stations_After_Treeline + PCNM1 + PCNM2",
#  "Stations_After_Treeline + PCNM1 * PCNM2",
#  "PCNM1 + PCNM2 + Direction",
#  "Direction + PCNM1 * PCNM2",
#  "PCNM1 + PCNM2 + Stations_After_Treeline + Direction",
#  "Stations_After_Treeline + Direction + PCNM1 * PCNM2",
#  "PCNM1 + PCNM2 + dist_coast",
#  "dist_coast + PCNM1 * PCNM2",
#  "dist_coast + PCNM1 + PCNM2 + Stations_After_Treeline",
#  "dist_coast + Stations_After_Treeline + PCNM1 * PCNM2",
#  "dist_coast + Direction + PCNM1 + PCNM2",
#  "dist_coast + Direction + PCNM1 * PCNM2",
#  "dist_coast + Direction + Stations_After_Treeline + PCNM1 + PCNM2",
#  "dist_coast + Direction + Stations_After_Treeline + PCNM1 * PCNM2"
#)
# Extract AIC values from the aic_values matrix
#aic_values_pcnm <- sapply(mods_pcnm, extractAIC)
#aic_values_pcnm <-aic_values_pcnm[2,]
#delta_AIC_pcnm <-aic_values_pcnm - min(aic_values_pcnm) # calculate delta AIC
#AIC_weight_pcnm <-exp(-0.5 * delta_AIC_pcnm) / 
#  sum(exp(-0.5 * delta_AIC_pcnm)) #calculate AIC weights
# Create a data frame
#model_names <- paste0("m", 1:16, "_pcnm")
#comparison_table_pcnm <- data.frame(
#  Model = model_names,
#  Terms = pcnm_terms,
# AIC = format_numeric(aic_values_pcnm),
# Delta_AIC = format_numeric(delta_AIC_pcnm),
# Weight = format_numeric(AIC_weight_pcnm)
#)


# Sort by AIC
#comparison_table_pcnm <- comparison_table_pcnm[order(as.numeric(comparison_table_pcnm$AIC)), ]

# Print the resulting table
#print(comparison_table_pcnm)

# Save to CSV
#write.csv(
#  comparison_table_pcnm,
#  file = "results/pcnm_AIC.csv",
#  row.names = FALSE
#)


#####################################
#LMER
# Generate the LaTeX table
