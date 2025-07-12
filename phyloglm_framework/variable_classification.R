# Variable Classification and Type Detection System for PhyloGLM Framework
# This system automatically detects variable types and handles appropriate conversions

library(dplyr)

#' Define variable classifications for the cooperative breeding dataset
#' 
#' @return A list containing variable classifications by type
get_variable_classifications <- function() {
  list(
    # Binary variables (0/1 coding)
    binary = c(
      "HighConfidence_Coop",
      "FemaleSong_Agg01", 
      "Griesser2017FamilialLiving",
      "AnyCoopEqualsCoop",
      "AnyNoncoopEqualsNoncoop",
      "MeanCoopTie2Coop",
      "MeanCoopTie2Noncoop",
      "GeographicRegion_Jetz"  # Holarctic=0, Tropical=1
    ),
    
    # Categorical variables (3+ levels, need factor conversion)
    categorical = c(
      "Territory",              # 3 levels: 1, 2, 3
      "Migration_AVONET",       # 3 levels: 1, 2, 3
      "Territory_12vs3",        # 2 levels but needs factor: 1vs2, 3
      "TerritorialityWeakVsStrong", # 2 levels: 0, 1
      "Habitat_AVONET",         # Multiple habitat types
      "Trophic.Level_AVONET",   # Multiple trophic levels
      "Primary.Lifestyle_AVONET" # Multiple lifestyle categories
    ),
    
    # Continuous variables
    continuous = c(
      "logMass_AVONET",
      "PercentAbsLogWingDimorphism",
      "logMaleFemalePlumageDiffAbs", 
      "Centroid.Latitude_AVONET",
      "Beak.Length_Culmen_AVONET",
      "Wing.Length_AVONET",
      "Tail.Length_AVONET",
      "Range.Size_AVONET"
    ),
    
    # Special continuous (may need transformation)
    continuous_special = c(
      "Centroid.Latitude_AVONET",  # May want absolute value
      "MaleFemalePlumageDiff"      # May want absolute value
    ),
    
    # Control variables (typically included in most models)
    control = c(
      "logMass_AVONET"
    ),
    
    # Response variables (can be dependent variables)
    responses = c(
      "FemaleSong_Agg01",
      "HighConfidence_Coop",
      "Griesser2017FamilialLiving"
    ),
    
    # Primary predictors (main variables of interest)
    predictors = c(
      "HighConfidence_Coop",
      "FemaleSong_Agg01",
      "Griesser2017FamilialLiving",
      "AnyCoopEqualsCoop",
      "MeanCoopTie2Coop"
    )
  )
}

#' Detect variable type automatically based on data characteristics
#' 
#' @param data A dataframe containing the variables
#' @param var_name Character string of variable name
#' @param classifications Optional list of predefined classifications
#' @return Character string indicating variable type
detect_variable_type <- function(data, var_name, classifications = NULL) {
  
  if (!var_name %in% colnames(data)) {
    stop(paste("Variable", var_name, "not found in data"))
  }
  
  # Use predefined classifications if available
  if (!is.null(classifications)) {
    for (type in names(classifications)) {
      if (var_name %in% classifications[[type]]) {
        return(type)
      }
    }
  }
  
  # Automatic detection based on data characteristics
  var_data <- data[[var_name]]
  var_data <- var_data[!is.na(var_data)]  # Remove NAs for analysis
  
  unique_vals <- unique(var_data)
  n_unique <- length(unique_vals)
  
  # Binary detection
  if (n_unique == 2 && all(sort(unique_vals) == c(0, 1))) {
    return("binary")
  }
  
  # Small number of integer values suggests categorical
  if (is.numeric(var_data) && all(var_data == as.integer(var_data)) && n_unique <= 10) {
    return("categorical")
  }
  
  # Character or factor variables are categorical
  if (is.character(var_data) || is.factor(var_data)) {
    return("categorical")
  }
  
  # Otherwise assume continuous
  if (is.numeric(var_data)) {
    return("continuous")
  }
  
  # Default
  return("unknown")
}

#' Prepare variable for analysis based on its type
#' 
#' @param data A dataframe containing the variable
#' @param var_name Character string of variable name
#' @param var_type Character string indicating variable type
#' @param reference_level For categorical variables, which level to use as reference
#' @return The prepared variable as a vector
prepare_variable <- function(data, var_name, var_type, reference_level = NULL) {
  
  var_data <- data[[var_name]]
  
  switch(var_type,
    "binary" = {
      # Ensure 0/1 coding
      if (all(var_data %in% c(0, 1, NA))) {
        return(as.numeric(var_data))
      } else {
        # Convert to 0/1 if needed
        unique_vals <- unique(var_data[!is.na(var_data)])
        if (length(unique_vals) == 2) {
          return(as.numeric(var_data == max(unique_vals, na.rm = TRUE)))
        }
      }
    },
    
    "categorical" = {
      # Convert to factor with meaningful levels
      if (var_name == "Territory") {
        return(factor(var_data, levels = c(1, 2, 3), 
                     labels = c("Weak", "Moderate", "Strong")))
      } else if (var_name == "Migration_AVONET") {
        return(factor(var_data, levels = c(1, 2, 3),
                     labels = c("Sedentary", "Partial", "Full")))
      } else if (var_name == "TerritorialityWeakVsStrong") {
        return(factor(var_data, levels = c(0, 1),
                     labels = c("Weak", "Strong")))
      } else {
        # Generic factor conversion
        f <- factor(var_data)
        if (!is.null(reference_level) && reference_level %in% levels(f)) {
          f <- relevel(f, ref = reference_level)
        }
        return(f)
      }
    },
    
    "continuous" = {
      # Return as numeric, possibly with transformations
      if (var_name == "Centroid.Latitude_AVONET") {
        # Option for absolute latitude
        return(as.numeric(var_data))
      } else if (var_name == "MaleFemalePlumageDiff") {
        # Option for absolute difference
        return(abs(as.numeric(var_data)))
      } else {
        return(as.numeric(var_data))
      }
    },
    
    "continuous_special" = {
      # Handle special transformations
      if (var_name == "Centroid.Latitude_AVONET") {
        return(abs(as.numeric(var_data)))  # Absolute latitude
      } else {
        return(as.numeric(var_data))
      }
    },
    
    # Default: return as-is
    return(var_data)
  )
}

#' Get variable label for plotting
#' 
#' @param var_name Character string of variable name
#' @return Character string with readable label
get_variable_label <- function(var_name) {
  labels <- list(
    "HighConfidence_Coop" = "Cooperative Breeding",
    "FemaleSong_Agg01" = "Female Song",
    "Griesser2017FamilialLiving" = "Familial Living",
    "Territory" = "Territory (3-level)",
    "Migration_AVONET" = "Migration",
    "TerritorialityWeakVsStrong" = "Territory (Weak/Strong)",
    "logMass_AVONET" = "log(Body Mass)",
    "PercentAbsLogWingDimorphism" = "Wing Dimorphism (%)",
    "logMaleFemalePlumageDiffAbs" = "Plumage Dimorphism",
    "Centroid.Latitude_AVONET" = "Latitude",
    "GeographicRegion_Jetz" = "Geographic Region"
  )
  
  if (var_name %in% names(labels)) {
    return(labels[[var_name]])
  } else {
    # Create readable label from variable name
    return(gsub("_", " ", gsub("\\.", " ", var_name)))
  }
}

#' Validate variable for phyloglm analysis
#' 
#' @param data A dataframe containing the variable
#' @param var_name Character string of variable name
#' @param var_type Character string indicating variable type
#' @return List with validation results
validate_variable <- function(data, var_name, var_type) {
  
  var_data <- data[[var_name]]
  n_total <- length(var_data)
  n_missing <- sum(is.na(var_data))
  n_complete <- n_total - n_missing
  
  unique_vals <- unique(var_data[!is.na(var_data)])
  n_unique <- length(unique_vals)
  
  # Basic validation
  issues <- character(0)
  warnings <- character(0)
  
  # Check for sufficient data
  if (n_complete < 50) {
    issues <- c(issues, paste("Very few complete cases:", n_complete))
  } else if (n_complete < 100) {
    warnings <- c(warnings, paste("Limited complete cases:", n_complete))
  }
  
  # Check missing data proportion
  missing_prop <- n_missing / n_total
  if (missing_prop > 0.5) {
    issues <- c(issues, paste("High missing data:", round(missing_prop * 100, 1), "%"))
  } else if (missing_prop > 0.2) {
    warnings <- c(warnings, paste("Moderate missing data:", round(missing_prop * 100, 1), "%"))
  }
  
  # Type-specific validation
  if (var_type == "binary") {
    if (n_unique != 2) {
      issues <- c(issues, paste("Binary variable has", n_unique, "unique values"))
    }
    # Check for unbalanced binary
    if (n_unique == 2) {
      tab <- table(var_data, useNA = "no")
      min_prop <- min(tab) / sum(tab)
      if (min_prop < 0.05) {
        warnings <- c(warnings, paste("Very unbalanced binary variable: min category =", 
                                     round(min_prop * 100, 1), "%"))
      }
    }
  } else if (var_type == "categorical") {
    if (n_unique < 2) {
      issues <- c(issues, "Categorical variable has < 2 levels")
    } else if (n_unique > 10) {
      warnings <- c(warnings, paste("Many categorical levels:", n_unique))
    }
  } else if (var_type == "continuous") {
    if (n_unique < 10) {
      warnings <- c(warnings, paste("Continuous variable has few unique values:", n_unique))
    }
  }
  
  return(list(
    valid = length(issues) == 0,
    issues = issues,
    warnings = warnings,
    n_complete = n_complete,
    n_missing = n_missing,
    missing_proportion = missing_prop,
    n_unique = n_unique,
    unique_values = if (n_unique <= 10) unique_vals else paste(n_unique, "values")
  ))
}

#' Create a summary of all variables in dataset
#' 
#' @param data A dataframe
#' @param variables Optional vector of variable names to analyze
#' @return Dataframe with variable summaries
summarize_variables <- function(data, variables = NULL) {
  
  if (is.null(variables)) {
    variables <- colnames(data)
  }
  
  classifications <- get_variable_classifications()
  
  summaries <- data.frame(
    variable = character(0),
    type = character(0),
    label = character(0),
    n_complete = integer(0),
    n_missing = integer(0),
    missing_prop = numeric(0),
    n_unique = integer(0),
    validation_status = character(0),
    issues = character(0),
    stringsAsFactors = FALSE
  )
  
  for (var in variables) {
    if (var %in% colnames(data)) {
      var_type <- detect_variable_type(data, var, classifications)
      validation <- validate_variable(data, var, var_type)
      
      summaries <- rbind(summaries, data.frame(
        variable = var,
        type = var_type,
        label = get_variable_label(var),
        n_complete = validation$n_complete,
        n_missing = validation$n_missing,
        missing_prop = round(validation$missing_proportion, 3),
        n_unique = validation$n_unique,
        validation_status = if (validation$valid) "OK" else "Issues",
        issues = paste(c(validation$issues, validation$warnings), collapse = "; "),
        stringsAsFactors = FALSE
      ))
    }
  }
  
  return(summaries)
}