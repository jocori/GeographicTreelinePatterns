### changing data to include at least the last 10 stations of each transect, and treeline
#set working directory
setwd("C:/Users/j364g325/Documents/mountain peaks")
# Read in your data
peaks <- read.csv("Reg_Master.csv", header = TRUE)
#in case there are "X" columns at the end
peaks<-peaks[,-c(47:53)]
#125 are missing directions
missing.dirs<-peaks[peaks$dir == "",]
#write.csv(missing.dirs,"missing.dirs.csv")

##get rid of them for now...figure out working code first!
peaks<-peaks[peaks$dir!="",]
#get rid of NA values
peaks<-na.omit(peaks)
# Extract years columns
year_columns <- grep("^Y[0-9]{4}$", names(peaks), value = TRUE)

# Function to create nested list for a single peak
create_peak_list <- function(peak_data) {
  # Split by direction
  transects <- split(peak_data, peak_data$dir)
  
  # For each transect, split by year and keep only relevant columns
  lapply(transects, function(transect) {
    lapply(year_columns, function(yr) {
      transect_data <- transect[c("wkt_geom", "Peak_ID", "Country", "dir", "dis_m_", "elevation", "long", "lat", "TreesBegin", yr)]
      names(transect_data)[ncol(transect_data)] <- "NDVI"
      return(transect_data)
    })
  })
}

# Create nested list structure for all peaks
nested_peaks <- lapply(split(peaks, peaks$Name), create_peak_list)

#function to select relevant stations with sorting
select_stations <- function(df) {
  # Sort the dataframe by distance from the peak (assuming 'dis_m_' is the column name)
  df <- df[order(df$dis_m_), ]
  
  # Find the index where TreesBegin is 1
  treeline_index <- which(df$TreesBegin == 1)
  
  if (length(treeline_index) > 0) {
    # Calculate stations after treeline
    stations_after_treeline <- nrow(df) - treeline_index[1] + 1
    stations_to_use <- max(10, 2 * stations_after_treeline + 1)
    
    # Select the relevant stations
    if (stations_to_use >= nrow(df)) {
      selected_df <- df
    } else {
      selected_df <- df[(nrow(df) - stations_to_use + 1):nrow(df), ]
    }
  } else {
    # If no treeline, select last 10 stations
    selected_df <- tail(df, 10)
  }
  
  return(selected_df)
}

# Apply the updated selection function to each year in each direction in each peak
selected_nested_peaks <- lapply(nested_peaks, function(peak) {
  lapply(peak, function(direction) {
    lapply(direction, select_stations)
  })
})

###now it's regression time!


#####keep all three models while assigning the best
perform_regressions <- function(df, peak_name, direction, year) {
  if (!is.data.frame(df) || nrow(df) <= 1) {
    return(data.frame(
      Peak = rep(peak_name, 3),
      Direction = rep(direction, 3),
      Year = rep(year, 3),
      Regression_Type = c("Linear", "Reciprocal_Linear", "Reciprocal_Quadratic"),
      Significance = rep(NA, 3),
      Intercept = rep(NA, 3),
      Slope = rep(NA, 3),
      AIC = rep(NA, 3),
      Best_Model = rep(NA, 3),
      Stations_After_Treeline = rep(NA, 3)
    ))
  }
  
  # Exclude rows where NDVI is NA or zero
  df <- df[!is.na(df$NDVI) & df$NDVI != 0,]
  
  # If there are not enough valid rows after filtering, return NA for all values
  if (nrow(df) <= 1) {
    return(data.frame(
      Peak = rep(peak_name, 3),
      Direction = rep(direction, 3),
      Year = rep(year, 3),
      Regression_Type = c("Linear", "Reciprocal_Linear", "Reciprocal_Quadratic"),
      Significance = rep(NA, 3),
      Intercept = rep(NA, 3),
      Slope = rep(NA, 3),
      AIC = rep(NA, 3),
      Best_Model = rep(NA, 3),
      Stations_After_Treeline = rep(NA, 3)
    ))
  }
  



  # Fit the models
  model_linear <- lm(NDVI ~ elevation, data = df)
  model_reciprocal_linear <- lm(1/NDVI ~ elevation, data = df)
  model_reciprocal_quadratic <- lm(1/NDVI ~ I(1/elevation^2) + elevation, data = df)
  
  # Calculate AICs
  aic_linear <- AIC(model_linear)
  aic_reciprocal_linear <- AIC(model_reciprocal_linear)
  aic_reciprocal_quadratic <- AIC(model_reciprocal_quadratic)
  
  # Calculate coefficients and significance for each model
  coef_linear <- coef(model_linear)
  coef_reciprocal_linear <- coef(model_reciprocal_linear)
  coef_reciprocal_quadratic <- coef(model_reciprocal_quadratic)
  
  significance_linear <- summary(model_linear)$coefficients[2, 'Pr(>|t|)']
  significance_reciprocal_linear <- summary(model_reciprocal_linear)$coefficients[2, 'Pr(>|t|)']
  significance_reciprocal_quadratic <- summary(model_reciprocal_quadratic)$coefficients[2, 'Pr(>|t|)']
  
  # Determine the best model by AIC
  aic_values <- c(Linear = aic_linear, Reciprocal_Linear = aic_reciprocal_linear, Reciprocal_Quadratic = aic_reciprocal_quadratic)
  best_model_type <- names(which.min(aic_values))
  
  # Construct and return a dataframe with results for all model types
  results_df <- data.frame(
    Peak = rep(peak_name, 3),
    Direction = rep(direction, 3),
    Year = rep(year, 3),
    Regression_Type = c("Linear", "Reciprocal_Linear", "Reciprocal_Quadratic"),
    Significance = c(significance_linear, significance_reciprocal_linear, significance_reciprocal_quadratic),
    Intercept = c(coef_linear[1], coef_reciprocal_linear[1], coef_reciprocal_quadratic[1]),
    Slope = c(coef_linear[2], coef_reciprocal_linear[2], coef_reciprocal_quadratic[2]),
    AIC = c(aic_linear, aic_reciprocal_linear, aic_reciprocal_quadratic),
    Best_Model = rep(best_model_type, 3),
    Stations_After_Treeline = rep(nrow(df), 3)
  )
  
  return(results_df)
}

all_results <- list()

# to assign correct years
yearsss <- 1982:2017

for (peak_name in names(selected_nested_peaks)) {
  peak_data <- selected_nested_peaks[[peak_name]]
  for (direction in names(peak_data)) {
    direction_data <- peak_data[[direction]]
    for (year in 1:length(direction_data)) {
      year_data <- direction_data[[year]]
      # Perform regressions and store results
      regression_results <- perform_regressions(year_data, peak_name, direction, yearsss[year])
      all_results <- c(all_results, list(regression_results))
    }
  }
}

# Combine all results into a single dataframe
final_results_df <- do.call(rbind, all_results) #this works but only preserves rows for the best model


write.csv(final_results_df,"peak_regressions_atLeastLast10Stations.csv")







######most of this works, but the output isn't correct........I am too terrified to delete it though

library(dplyr)

# Function to perform regression and return results
perform_regressions <- function(df, peak_name, direction, year, type) {
  if (!is.data.frame(df) || nrow(df) <= 1) {
    return(data.frame(
      Peak = peak_name,
      Direction = direction,
      Year = year,
      Regression_Type = type,
      Significance = NA,
      Intercept = NA,
      Slope = NA,
      AIC = NA,
      Stations_After_Treeline = NA
    ))
  }
  
  # Count the number of stations after treeline
  treeline_index <- which.max(df$TreesBegin)
  stations_after_treeline <- ifelse(is.na(treeline_index), 0, nrow(df) - treeline_index)
  
  # Exclude rows with NDVI NA or zero
  df <- df[!is.na(df$NDVI) & df$NDVI != 0,]
  
  # Check if there are enough rows for regression
  if (nrow(df) > 1) {
    # Fit the model
    model <- switch(type,
                    "Linear" = lm(NDVI ~ elevation, data = df),
                    "Reciprocal_Linear" = lm(1/NDVI ~ elevation, data = df),
                    "Reciprocal_Quadratic" = lm(1/NDVI ~ I(1/elevation^2) + elevation, data = df)
    )
    
    # Extract model details
    aic_model <- AIC(model)
    coef_model <- coef(model)
    significance <- summary(model)$coefficients[2, 'Pr(>|t|)']
    
    # Return model details
    return(data.frame(
      Peak = peak_name,
      Direction = direction,
      Year = year,
      Regression_Type = type,
      Significance = significance,
      Intercept = coef_model[1],
      Slope = coef_model[2],
      AIC = aic_model,
      Best_by_AIC = 
        Stations_After_Treeline = stations_after_treeline
    ))
  } else {
    return(data.frame(
      Peak = peak_name,
      Direction = direction,
      Year = year,
      Regression_Type = type,
      Significance = NA,
      Intercept = NA,
      Slope = NA,
      AIC = NA,
      Stations_After_Treeline = stations_after_treeline
    ))
  }
}

# Example usage
# Assuming 'selected_nested_peaks' is your nested list of dataframes
all_results <- list()

# to assign correct years
yearsss <- 1982:2017

for (peak_name in names(selected_nested_peaks)) {
  peak_data <- selected_nested_peaks[[peak_name]]
  for (direction in names(peak_data)) {
    direction_data <- peak_data[[direction]]
    for (year in 1:length(direction_data)) {
      year_data <- direction_data[[year]]
      # Perform regressions for each type and store results
      for (type in c("Linear", "Reciprocal_Linear", "Reciprocal_Quadratic")) {
        regression_results <- perform_regressions(year_data, peak_name, direction, year, type)
        regression_results$Year <- yearsss[year]
        all_results <- c(all_results, list(regression_results))
      }
    }
  }
}

# Combine all results into a single dataframe
final_results_df <- do.call(rbind, all_results)



##############adding best model
perform_regressions <- function(df, peak_name, direction, year) {
  if (!is.data.frame(df) || nrow(df) <= 1) {
    return(data.frame(
      Peak = peak_name,
      Direction = direction,
      Year = year,
      Regression_Type = NA,
      Significance = NA,
      Intercept = NA,
      Slope = NA,
      AIC = NA,
      Best_Model = NA,
      Stations_After_Treeline = NA
    ))
  }
  
  # Count the number of stations after treeline
  treeline_index <- which.max(df$TreesBegin)
  stations_after_treeline <- ifelse(is.na(treeline_index), 0, nrow(df) - treeline_index)
  
  # Exclude rows with NDVI NA or zero
  df <- df[!is.na(df$NDVI) & df$NDVI != 0,]
  
  # Check if there are enough rows for regression
  if (nrow(df) > 1) {
    # Fit the models
    model_linear <- lm(NDVI ~ elevation, data = df)
    model_reciprocal_linear <- lm(1/NDVI ~ elevation, data = df)
    model_reciprocal_quadratic <- lm(1/NDVI ~ I(1/elevation^2) + elevation, data = df)
    
    # Calculate AICs
    aic_linear <- AIC(model_linear)
    aic_reciprocal_linear <- AIC(model_reciprocal_linear)
    aic_reciprocal_quadratic <- AIC(model_reciprocal_quadratic)
    
    # Determine the best model by AIC
    aic_values <- c(Linear = aic_linear, Reciprocal_Linear = aic_reciprocal_linear, Reciprocal_Quadratic = aic_reciprocal_quadratic)
    best_model_type <- names(which.min(aic_values))
    
    # Select the best model based on AIC
    best_model <- switch(best_model_type,
                         "Linear" = model_linear,
                         "Reciprocal_Linear" = model_reciprocal_linear,
                         "Reciprocal_Quadratic" = model_reciprocal_quadratic)
    
    # Extract model details
    coef_model <- coef(best_model)
    significance <- summary(best_model)$coefficients[2, 'Pr(>|t|)']
    
    # Return model details
    return(data.frame(
      Peak = peak_name,
      Direction = direction,
      Year = year,
      Regression_Type = best_model_type,
      Significance = significance,
      Intercept = coef_model[1],
      Slope = coef_model[2],
      AIC = aic_values[best_model_type],
      Best_Model = best_model_type,
      Stations_After_Treeline = stations_after_treeline
    ))
  } else {
    return(data.frame(
      Peak = peak_name,
      Direction = direction,
      Year = year,
      Regression_Type = NA,
      Significance = NA,
      Intercept = NA,
      Slope = NA,
      AIC = NA,
      Best_Model = NA,
      Stations_After_Treeline = stations_after_treeline
    ))
  }
}
all_results <- list()

# to assign correct years
yearsss <- 1982:2017

for (peak_name in names(selected_nested_peaks)) {
  peak_data <- selected_nested_peaks[[peak_name]]
  for (direction in names(peak_data)) {
    direction_data <- peak_data[[direction]]
    for (year in 1:length(direction_data)) {
      year_data <- direction_data[[year]]
      # Perform regressions and store results
      regression_results <- perform_regressions(year_data, peak_name, direction, yearsss[year])
      all_results <- c(all_results, list(regression_results))
    }
  }
}

# Combine all results into a single dataframe
final_results_df <- do.call(rbind, all_results) #this works but only preserves rows for the best model