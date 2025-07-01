################################################################################
# Project: NDVI of North American Mountain Peaks through time
# Authors: Joanna Corimanya, Daniel JimEnez-Garcia, A. Townsend Peterson
# Date: 09/11/2023 (dd/mm/yyyy)
# Funding:
################################################################################
#set working directory
setwd("C:/Users/j364g325/Documents/mountain peaks")
#read in data
peaks<- read.csv("Reg_Master.csv", header = TRUE)
head(peaks)
str(peaks)
#in case there are "X" columns at the end
#peaks<-peaks[,-c(47:53)]
allp <- unique(peaks$Name)
alld <- na.omit(unique(peaks$dir))
aaly <- colnames(peaks)[5:40]
#el<-na.omit(peaks$elevation)
#sev <- na.omit(peaks$Y2017)
# Create a directory to save the graphs (change the path as needed)
output_directory <- "C:/Users/j364g325/Documents/mountain peaks/peak_graphs_png"
dir.create(output_directory, showWarnings = FALSE)

# Loop through each combination of peak, direction, and year
#######CHANGE TO ONLY 2017#########
#######change to whole dataset########
for (Name in allp) {
  for (dir in alld) {
      # Create a unique filename for each combination
      filename <- paste(Name, dir, ".png", sep = "_")
      
      # Subset the data for the current combination
      subsetex <- na.omit(peaks[peaks$Name == Name & peaks$dir == dir, c("elevation", "Y2017")])
      
      # Check if there is any data to create a plot
      if (nrow(subsetex) > 0) {
      
      # Save the plot to a file in the output directory
      png(file.path(output_directory, filename))
        # Create a plot for the current combination
      plot(subsetex, main = paste("Peak:", Name, "Direction:", dir), ylim = c(-0.3,0.8))

      dev.off()
    }
  }
}

########DISTANCE########
output_directory <- "C:/Peaks/peak_graphs_dist"
dir.create(output_directory, showWarnings = FALSE)

# Loop through each combination of peak, direction, and year
#######CHANGE TO ONLY 2017#########
#######change to whole dataset########
for (Name in allp) {
  for (dir in alld) {
    # Create a unique filename for each combination
    filename <- paste(Name, dir, ".png", sep = "_")
    
    # Subset the data for the current combination
    subsetex <- na.omit(peaks[peaks$Name == Name & peaks$dir == dir, c("dis_m_", "Y2017")])
    
    # Check if there is any data to create a plot
    if (nrow(subsetex) > 0) {
      
      # Save the plot to a file in the output directory
      png(file.path(output_directory, filename))
      # Create a plot for the current combination
      plot(subsetex, main = paste("Peak:", Name, "Direction:", dir), ylim = c(-0.3,0.8) )
      
      dev.off()
    }
  }
}










#multiple years, for regression
allp <- unique(peaks$Name)
alld <- na.omit(unique(peaks$dir))
aaly <- colnames(peaks)[-c(1:8)]
#list(testpeaks, peaks)
for (Name in allp) {
  for (dir in alld) {
    for (yr in aaly) {
      # Create a unique filename for each combination
      #filename <- paste(Name, dir, yr, ".pdf", sep = "_")
      
      # Subset the data for the current combination
      subsetex <- na.omit(peaks[peaks$Name == Name & peaks$dir == dir, c("elevation", yr)])
      
      # Create a plot for the current combination
      
      # create empty dataframe
      peak_models <- data.frame()
      #pdf(file.path(output_directory, filename))
      lm(subsetex ~ peaks$elevation, main = paste("Peak:", Name, "Direction:", dir, "Year:", yr))
      peak_models$
      dev.off()
    }
  }
}
# Initialize an empty data frame to store the results
results <- data.frame(
  Peak = character(0),
  Direction = character(0),
  Regression_Type = character(0),
  Significance = numeric(0),
  Parameters = character(0),
  Best_via_AIC = character(0)
)

# Loop through all combinations of Name, dir, and yr
for (Name in allp) {
  for (dir in alld) {
    for (yr in aaly) {
      # Subset the data for the current combination
      subsetex <- na.omit(peaks[peaks$Name == Name & peaks$dir == dir, c("elevation", yr)])
      
      if (nrow(subsetex) > 1) {  # Check if there's enough data to fit a model
        # Fit linear regression
        linear_model <- lm(subsetex[[2]] ~ subsetex[[1]])
        
        # Fit reciprocal-linear regression
        reciprocal_linear_model <- lm(1/subsetex[[2]] ~ subsetex[[1]])
        
        # Fit reciprocal-quadratic regression
        reciprocal_quadratic_model <- lm(1/subsetex[[2]] ~ I(1/subsetex[[2]]^2) + subsetex[[1]])
        
        # Calculate AIC values for each model
        aic_linear <- AIC(linear_model)
        aic_reciprocal_linear <- AIC(reciprocal_linear_model)
        aic_reciprocal_quadratic <- AIC(reciprocal_quadratic_model)
        
        # Determine the best model via AIC
        best_model <- which.min(c(aic_linear, aic_reciprocal_linear, aic_reciprocal_quadratic))
        
        # Store the results in the data frame
        results <- rbind(results, data.frame(
          Peak = Name,
          Direction = dir,
          Regression_Type = c("Linear", "Reciprocal-Linear", "Reciprocal-Quadratic")[best_model],
          Significance = summary(linear_model)$coefficients[2, 'Pr(>|t|)'],
          Parameters = paste(coef(linear_model), collapse = ", "),
          Best_via_AIC = c("Linear", "Reciprocal-Linear", "Reciprocal-Quadratic")[best_model]
        ))
      }
    }
  }
}

# Print or save the results data frame
print(results)

# Initialize an empty data frame to store the results
results <- data.frame(
  Peak = character(0),
  Direction = character(0),
  Regression_Type = character(0),
  Significance = numeric(0),
  Parameters = character(0),
  Best_via_AIC = character(0)
)

# Loop through all combinations of Name, dir, and yr
for (Name in allp) {
  for (dir in alld) {
    for (yr in aaly) {
      # Subset the data for the current combination
      subsetex <- na.omit(peaks[peaks$Name == Name & peaks$dir == dir, c("elevation", yr)])
      
      # Check if there are any missing or infinite values in the response variable
      if (sum(!complete.cases(subsetex)) == 0 && nrow(subsetex) > 1) {
        tryCatch({
          # Fit linear regression
          linear_model <- lm(subsetex[[2]] ~ subsetex[[1]])
          
          # Fit reciprocal-linear regression
          reciprocal_linear_model <- lm(1/subsetex[[2]] ~ subsetex[[1]])
          
          # Fit reciprocal-quadratic regression
          reciprocal_quadratic_model <- lm(1/subsetex[[2]] ~ I(1/subsetex[[2]]^2) + subsetex[[1]])
          
          # Calculate AIC values for each model
          aic_linear <- AIC(linear_model)
          aic_reciprocal_linear <- AIC(reciprocal_linear_model)
          aic_reciprocal_quadratic <- AIC(reciprocal_quadratic_model)
          
          # Determine the best model via AIC
          best_model <- which.min(c(aic_linear, aic_reciprocal_linear, aic_reciprocal_quadratic))
          
          # Store the results in the data frame
          results <- rbind(results, data.frame(
            Peak = Name,
            Direction = dir,
            Regression_Type = c("Linear", "Reciprocal-Linear", "Reciprocal-Quadratic")[best_model],
            Significance = summary(linear_model)$coefficients[2, 'Pr(>|t|)'],
            Parameters = paste(coef(linear_model), collapse = ", "),
            Best_via_AIC = c("Linear", "Reciprocal-Linear", "Reciprocal-Quadratic")[best_model]
          ))
        }, error = function(e) {
          # Handle errors (e.g., NA/NaN/Inf in 'y') by continuing to the next iteration
          cat("Error fitting a model for Peak:", Name, "Direction:", dir, "Year:", yr, "\n")
          next
        })
      }
    }
  }
}

# Print or save the results data frame
print(results)

# Initialize an empty data frame to store the results
results <- data.frame(
  Peak = character(0),
  Year = numeric(0),
  Direction = character(0),
  Regression_Type = character(0),
  Significance = numeric(0),
  Parameters = character(0),
  Best_via_AIC = character(0)
)

# Loop through all combinations of Name, dir, and yr
for (Name in allp) {
  for (dir in alld) {
    for (yr in aaly) {
      # Extract the year from 'yr' (assuming 'yr' is in the format "Yyyyy")
      year <- as.integer(substr(yr, 2, 5))
      # Subset the data for the current combination
      subsetex <- na.omit(peaks[peaks$Name == Name & peaks$dir == dir, c("elevation", yr)])
      
      # Check if there are any missing or infinite values in the response variable
      if (sum(!complete.cases(subsetex)) == 0 && nrow(subsetex) > 1) {
        # Attempt to fit the models
        try({
          # Fit linear regression
          linear_model <- lm(subsetex[[2]] ~ subsetex[[1]])
          
          # Fit reciprocal-linear regression
          reciprocal_linear_model <- lm(1/subsetex[[2]] ~ subsetex[[1]])
          
          # Fit reciprocal-quadratic regression
          reciprocal_quadratic_model <- lm(1/subsetex[[2]] ~ I(1/subsetex[[2]]^2) + subsetex[[1]])
          
          # Calculate AIC values for each model
          aic_linear <- AIC(linear_model)
          aic_reciprocal_linear <- AIC(reciprocal_linear_model)
          aic_reciprocal_quadratic <- AIC(reciprocal_quadratic_model)
          
          # Determine the best model via AIC
          best_model <- which.min(c(aic_linear, aic_reciprocal_linear, aic_reciprocal_quadratic))
          
          # Store the results in the data frame
          results <- rbind(results, data.frame(
            Peak = Name,
            Year = year,
            Direction = dir,
            Regression_Type = c("Linear", "Reciprocal-Linear", "Reciprocal-Quadratic")[best_model],
            Significance = summary(linear_model)$coefficients[2, 'Pr(>|t|)'],
            Parameters = paste(coef(linear_model), collapse = ", "),
            Best_via_AIC = c("Linear", "Reciprocal-Linear", "Reciprocal-Quadratic")[best_model]
          ))
        }, silent = TRUE)  # Use silent = TRUE to suppress error messages
      }
    }
  }
}

# Print or save the results data frame
print(results)


# Initialize an empty data frame to store the results
results <- data.frame(
  Peak = character(0),
  Direction = character(0),
  Year = character(0),  # Add a Year column
  Regression_Type = character(0),
  Significance = numeric(0),
  Parameters = character(0),
  AIC = numeric(0)  # Add AIC column
)

# Initialize an empty data frame to store the results
results <- data.frame(
  Peak = character(0),
  Direction = character(0),
  Year = character(0),
  Regression_Type = character(0),
  Significance = numeric(0),
  Parameters = character(0),
  AIC = numeric(0)
)

# Initialize an empty data frame to store the results
results <- data.frame(
  Peak = character(0),
  Direction = character(0),
  Year = character(0),
  Regression_Type = character(0),
  Significance = numeric(0),
  Parameters = character(0),
  AIC = numeric(0)
)

# Loop through all combinations of Name, dir, and yr
for (Name in allp) {
  for (dir in alld) {
    for (yr in aaly) {
      # Extract the year from 'yr' (assuming 'yr' is in the format "Yyyyy")
      year <- as.integer(substr(yr, 2, 5))
      
      # Subset the data for the current combination
      subsetex <- na.omit(peaks[peaks$Name == Name & peaks$dir == dir, c("elevation", yr)])
      
      # Check if there are any missing or problematic values in the response variable
      if (sum(is.na(subsetex[[2]]) | is.nan(subsetex[[2]]) | is.infinite(subsetex[[2]])) == 0 &&
          nrow(subsetex) > 1) {
        # Fit linear regression
        linear_model <- tryCatch(
          lm(subsetex[[2]] ~ subsetex[[1]]),
          error = function(e) NULL
        )
        if (!is.null(linear_model)) {
          aic_linear <- AIC(linear_model)
          # Store the results for linear regression
          results <- rbind(results, data.frame(
            Peak = Name,
            Direction = dir,
            Year = year,  # Add Year information
            Regression_Type = "Linear",
            Significance = summary(linear_model)$coefficients[2, 'Pr(>|t|)'],
            Parameters = paste(coef(linear_model), collapse = ", "),
            AIC = aic_linear
          ))
        }
        
        # Fit reciprocal-linear regression
        reciprocal_linear_model <- tryCatch(
          lm(1/subsetex[[2]] ~ subsetex[[1]]),
          error = function(e) NULL
        )
        if (!is.null(reciprocal_linear_model)) {
          aic_reciprocal_linear <- AIC(reciprocal_linear_model)
          # Store the results for reciprocal-linear regression
          results <- rbind(results, data.frame(
            Peak = Name,
            Direction = dir,
            Year = year,  # Add Year information
            Regression_Type = "Reciprocal-Linear",
            Significance = summary(reciprocal_linear_model)$coefficients[2, 'Pr(>|t|)'],
            Parameters = paste(coef(reciprocal_linear_model), collapse = ", "),
            AIC = aic_reciprocal_linear
          ))
        }
        
        # Fit reciprocal-quadratic regression
        reciprocal_quadratic_model <- tryCatch(
          lm(1/subsetex[[2]] ~ I(1/subsetex[[2]]^2) + subsetex[[1]]),
          error = function(e) NULL
        )
        if (!is.null(reciprocal_quadratic_model)) {
          aic_reciprocal_quadratic <- AIC(reciprocal_quadratic_model)
          # Store the results for reciprocal-quadratic regression
          results <- rbind(results, data.frame(
            Peak = Name,
            Direction = dir,
            Year = year,  # Add Year information
            Regression_Type = "Reciprocal-Quadratic",
            Significance = summary(reciprocal_quadratic_model)$coefficients[2, 'Pr(>|t|)'],
            Parameters = paste(coef(reciprocal_quadratic_model), collapse = ", "),
            AIC = aic_reciprocal_quadratic
          ))
        }
      }
    }
  }
}

# Print or save the results data frame
print(results)



write.csv(results, file = "regressions.csv")


######################

