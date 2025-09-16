# Configuration Builder for PhyloGLM Framework
# Creates and manages analysis configurations

#library(yaml)

#' Create a single analysis configuration
#' 
#' @param response Response variable name
#' @param predictors Character vector of predictor variables
#' @param controls Character vector of control variables (optional)
#' @param transformations Named list of transformations (optional)
#' @param complexity_levels Which model complexity levels to include
#' @param max_interactions Maximum interaction order
#' @param name Optional name for the analysis
#' @return List with analysis configuration
create_analysis_config <- function(response,
                                  predictors,
                                  controls = NULL,
                                  transformations = NULL,
                                  complexity_levels = c("null", "main", "additive", "twoway", "threeway"),
                                  max_interactions = 3,
                                  name = NULL) {
  
  config <- list(
    response = response,
    predictors = predictors,
    controls = controls,
    transformations = transformations,
    complexity_levels = complexity_levels,
    max_interactions = max_interactions
  )
  
  if (!is.null(name)) {
    config$name <- name
  }
  
  class(config) <- c("phyloglm_config", "list")
  return(config)
}

#' Create multiple analysis configurations from a specification
#' 
#' @param spec List specifying multiple analyses
#' @return List of analysis configurations
create_analysis_configs_from_spec <- function(spec) {
  
  configs <- list()
  
  # If spec contains multiple analyses
  if ("analyses" %in% names(spec)) {
    for (i in seq_along(spec$analyses)) {
      analysis <- spec$analyses[[i]]
      config <- do.call(create_analysis_config, analysis)
      configs[[length(configs) + 1]] <- config
    }
  } else {
    # Single analysis
    configs[[1]] <- do.call(create_analysis_config, spec)
  }
  
  return(configs)
}

#' Create standard analysis configurations for cooperative breeding
#' 
#' @param include_sets Which standard sets to include
#' @return List of analysis configurations
create_standard_configs <- function(include_sets = "all") {
  
  configs <- list()
  
  # ========== FEMALE SONG AS RESPONSE ==========
  
  # 1. FS vs Coop + TerrWS + Mass
  configs$FS_vs_CB_TerrWS_Mass <- create_analysis_config(
    response = "FemaleSong_Agg01",
    predictors = c("HighConfidence_Coop", "TerritorialityWeakVsStrong"),
    controls = c("logMass_AVONET"),
    complexity_levels = c("null", "main", "additive", "twoway", "threeway"),
    max_interactions = 3,
    name = "FS_vs_CB_TerrWS_Mass"
  )
  
  # 2. FS vs Coop + TerrWS + abs(Latitude)
  configs$FS_vs_CB_TerrWS_absLat <- create_analysis_config(
    response = "FemaleSong_Agg01",
    predictors = c("HighConfidence_Coop", "TerritorialityWeakVsStrong"),
    controls = c("Centroid.Latitude_AVONET"),
    transformations = list(Centroid.Latitude_AVONET = "abs"),
    complexity_levels = c("null", "main", "additive", "twoway", "threeway"),
    max_interactions = 3,
    name = "FS_vs_CB_TerrWS_absLat"
  )
  
  # 3. FS vs Coop + TerrWS + PlumageDim
  configs$FS_vs_CB_TerrWS_PlumDim <- create_analysis_config(
    response = "FemaleSong_Agg01",
    predictors = c("HighConfidence_Coop", "TerritorialityWeakVsStrong"),
    controls = c("logMaleFemalePlumageDiffAbs"),
    complexity_levels = c("null", "main", "additive", "twoway", "threeway"),
    max_interactions = 3,
    name = "FS_vs_CB_TerrWS_PlumDim"
  )
  
  # 4. FS vs Coop + TerrWS + WingDim
  configs$FS_vs_CB_TerrWS_WingDim <- create_analysis_config(
    response = "FemaleSong_Agg01",
    predictors = c("HighConfidence_Coop", "TerritorialityWeakVsStrong"),
    controls = c("PercentAbsLogWingDimorphism"),
    complexity_levels = c("null", "main", "additive", "twoway", "threeway"),
    max_interactions = 3,
    name = "FS_vs_CB_TerrWS_WingDim"
  )
  
  # 5. FS vs Coop + Terr3 + Mass
  configs$FS_vs_CB_Terr3_Mass <- create_analysis_config(
    response = "FemaleSong_Agg01",
    predictors = c("HighConfidence_Coop", "Territory_12vs3"),
    controls = c("logMass_AVONET"),
    complexity_levels = c("null", "main", "additive", "twoway", "threeway"),
    max_interactions = 3,
    name = "FS_vs_CB_Terr3_Mass"
  )
  
  # 6. FS vs Coop + Migration + Mass
  configs$FS_vs_CB_Migration_Mass <- create_analysis_config(
    response = "FemaleSong_Agg01",
    predictors = c("HighConfidence_Coop", "Migration_AVONET"),
    controls = c("logMass_AVONET"),
    complexity_levels = c("null", "main", "additive", "twoway", "threeway"),
    max_interactions = 3,
    name = "FS_vs_CB_Migration_Mass"
  )
  
  # 7. FS vs Coop + TerrWS + Region
  configs$FS_vs_CB_TerrWS_Region <- create_analysis_config(
    response = "FemaleSong_Agg01",
    predictors = c("HighConfidence_Coop", "TerritorialityWeakVsStrong", "GeographicRegion_Jetz"),
    complexity_levels = c("null", "main", "additive", "twoway", "threeway"),
    max_interactions = 3,
    name = "FS_vs_CB_TerrWS_Region"
  )
  
  # 8. FS vs Coop + TerrWS + Region + Mass (4 predictors)
  configs$FS_vs_CB_TerrWS_Region_Mass <- create_analysis_config(
    response = "FemaleSong_Agg01",
    predictors = c("HighConfidence_Coop", "TerritorialityWeakVsStrong", "GeographicRegion_Jetz"),
    controls = c("logMass_AVONET"),
    complexity_levels = c("null", "main", "additive", "twoway", "threeway", "fourway"),
    max_interactions = 4,
    name = "FS_vs_CB_TerrWS_Region_Mass"
  )
  
  # 9. FS vs Coop + Fam + TerrWS
  configs$FS_vs_CB_Fam_TerrWS <- create_analysis_config(
    response = "FemaleSong_Agg01",
    predictors = c("HighConfidence_Coop", "Griesser2017FamilialLiving", "TerritorialityWeakVsStrong"),
    complexity_levels = c("null", "main", "additive", "twoway", "threeway"),
    max_interactions = 3,
    name = "FS_vs_CB_Fam_TerrWS"
  )
  
  # 10. FS vs Fam + TerrWS + Mass
  configs$FS_vs_Fam_TerrWS_Mass <- create_analysis_config(
    response = "FemaleSong_Agg01",
    predictors = c("Griesser2017FamilialLiving", "TerritorialityWeakVsStrong"),
    controls = c("logMass_AVONET"),
    complexity_levels = c("null", "main", "additive", "twoway", "threeway"),
    max_interactions = 3,
    name = "FS_vs_Fam_TerrWS_Mass"
  )
  
  # ========== COOPERATIVE BREEDING AS RESPONSE ==========
  
  # 11. CB vs FS + TerrWS + Mass
  configs$CB_vs_FS_TerrWS_Mass <- create_analysis_config(
    response = "HighConfidence_Coop",
    predictors = c("FemaleSong_Agg01", "TerritorialityWeakVsStrong"),
    controls = c("logMass_AVONET"),
    complexity_levels = c("null", "main", "additive", "twoway", "threeway"),
    max_interactions = 3,
    name = "CB_vs_FS_TerrWS_Mass"
  )
  
  # 12. CB vs FS + TerrWS + abs(Latitude)
  configs$CB_vs_FS_TerrWS_absLat <- create_analysis_config(
    response = "HighConfidence_Coop",
    predictors = c("FemaleSong_Agg01", "TerritorialityWeakVsStrong"),
    controls = c("Centroid.Latitude_AVONET"),
    transformations = list(Centroid.Latitude_AVONET = "abs"),
    complexity_levels = c("null", "main", "additive", "twoway", "threeway"),
    max_interactions = 3,
    name = "CB_vs_FS_TerrWS_absLat"
  )
  
  # 13. CB vs FS + TerrWS + PlumageDim
  configs$CB_vs_FS_TerrWS_PlumDim <- create_analysis_config(
    response = "HighConfidence_Coop",
    predictors = c("FemaleSong_Agg01", "TerritorialityWeakVsStrong"),
    controls = c("logMaleFemalePlumageDiffAbs"),
    complexity_levels = c("null", "main", "additive", "twoway", "threeway"),
    max_interactions = 3,
    name = "CB_vs_FS_TerrWS_PlumDim"
  )
  
  # 14. CB vs FS + TerrWS + WingDim
  configs$CB_vs_FS_TerrWS_WingDim <- create_analysis_config(
    response = "HighConfidence_Coop",
    predictors = c("FemaleSong_Agg01", "TerritorialityWeakVsStrong"),
    controls = c("PercentAbsLogWingDimorphism"),
    complexity_levels = c("null", "main", "additive", "twoway", "threeway"),
    max_interactions = 3,
    name = "CB_vs_FS_TerrWS_WingDim"
  )
  
  # 15. CB vs FS + Terr3 + Mass
  configs$CB_vs_FS_Terr3_Mass <- create_analysis_config(
    response = "HighConfidence_Coop",
    predictors = c("FemaleSong_Agg01", "Territory_12vs3"),
    controls = c("logMass_AVONET"),
    complexity_levels = c("null", "main", "additive", "twoway", "threeway"),
    max_interactions = 3,
    name = "CB_vs_FS_Terr3_Mass"
  )
  
  # 16. CB vs FS + Migration + Mass
  configs$CB_vs_FS_Migration_Mass <- create_analysis_config(
    response = "HighConfidence_Coop",
    predictors = c("FemaleSong_Agg01", "Migration_AVONET"),
    controls = c("logMass_AVONET"),
    complexity_levels = c("null", "main", "additive", "twoway", "threeway"),
    max_interactions = 3,
    name = "CB_vs_FS_Migration_Mass"
  )
  
  # 17. CB vs FS + TerrWS + Region
  configs$CB_vs_FS_TerrWS_Region <- create_analysis_config(
    response = "HighConfidence_Coop",
    predictors = c("FemaleSong_Agg01", "TerritorialityWeakVsStrong", "GeographicRegion_Jetz"),
    complexity_levels = c("null", "main", "additive", "twoway", "threeway"),
    max_interactions = 3,
    name = "CB_vs_FS_TerrWS_Region"
  )
  
  # 18. CB vs FS + TerrWS + Region + Mass (4 predictors)
  configs$CB_vs_FS_TerrWS_Region_Mass <- create_analysis_config(
    response = "HighConfidence_Coop",
    predictors = c("FemaleSong_Agg01", "TerritorialityWeakVsStrong", "GeographicRegion_Jetz"),
    controls = c("logMass_AVONET"),
    complexity_levels = c("null", "main", "additive", "twoway", "threeway", "fourway"),
    max_interactions = 4,
    name = "CB_vs_FS_TerrWS_Region_Mass"
  )
  
  # 19. CB vs FS + Fam + TerrWS
  configs$CB_vs_FS_Fam_TerrWS <- create_analysis_config(
    response = "HighConfidence_Coop",
    predictors = c("FemaleSong_Agg01", "Griesser2017FamilialLiving", "TerritorialityWeakVsStrong"),
    complexity_levels = c("null", "main", "additive", "twoway", "threeway"),
    max_interactions = 3,
    name = "CB_vs_FS_Fam_TerrWS"
  )
  
  # 20. CB vs Fam + TerrWS + Mass
  configs$CB_vs_Fam_TerrWS_Mass <- create_analysis_config(
    response = "HighConfidence_Coop",
    predictors = c("Griesser2017FamilialLiving", "TerritorialityWeakVsStrong"),
    controls = c("logMass_AVONET"),
    complexity_levels = c("null", "main", "additive", "twoway", "threeway"),
    max_interactions = 3,
    name = "CB_vs_Fam_TerrWS_Mass"
  )
  
  return(configs)
}

#' Save configurations to YAML file
#' 
#' @param configs List of configurations
#' @param file Path to output YAML file
# save_configs_to_yaml <- function(configs, file) {
#   # Convert to simpler structure for YAML
#   yaml_list <- list(
#     created = Sys.time(),
#     n_analyses = length(configs),
#     analyses = lapply(configs, function(c) {
#       list(
#         name = ifelse(is.null(c$name), "Unnamed", c$name),
#         response = c$response,
#         predictors = c$predictors,
#         controls = c$controls,
#         transformations = c$transformations,
#         complexity_levels = c$complexity_levels,
#         max_interactions = c$max_interactions
#       )
#     })
#   )
#   
#   write_yaml(yaml_list, file)
# }

#' Load configurations from YAML file
#' 
#' @param file Path to YAML file
#' @return List of configurations
# load_configs_from_yaml <- function(file) {
#   yaml_data <- read_yaml(file)
#   
#   configs <- list()
#   for (i in seq_along(yaml_data$analyses)) {
#     analysis <- yaml_data$analyses[[i]]
#     config <- create_analysis_config(
#       response = analysis$response,
#       predictors = analysis$predictors,
#       controls = analysis$controls,
#       transformations = analysis$transformations,
#       complexity_levels = analysis$complexity_levels,
#       max_interactions = analysis$max_interactions,
#       name = analysis$name
#     )
#     configs[[length(configs) + 1]] <- config
#   }
#   
#   return(configs)
# }

#' Create a configuration template file
#' 
#' @param file Path to output file
# create_config_template <- function(file = "phyloglm_config_template.yaml") {
#   
#   template <- '# PhyloGLM Analysis Configuration Template
# # Edit this file to specify your analyses
# 
# # Analysis specifications
# analyses:
#   # Example 1: Basic analysis
#   - name: "FemaleSong_vs_CoopBreeding_Basic"
#     response: "FemaleSong_Agg01"
#     predictors: 
#       - "HighConfidence_Coop"
#     controls:
#       - "logMass_AVONET"
#     complexity_levels:
#       - "null"
#       - "main"
#       - "additive"
#     max_interactions: 2
#     
#   # Example 2: Analysis with territorial moderation
#   - name: "FemaleSong_vs_CoopBreeding_Territorial"
#     response: "FemaleSong_Agg01"
#     predictors:
#       - "HighConfidence_Coop"
#       - "TerritorialityWeakVsStrong"
#     controls:
#       - "logMass_AVONET"
#     complexity_levels:
#       - "null"
#       - "main"
#       - "additive"
#       - "twoway"
#       - "threeway"
#     max_interactions: 3
#     
#   # Example 3: Analysis with transformations
#   - name: "FemaleSong_vs_CoopBreeding_Latitude"
#     response: "FemaleSong_Agg01"
#     predictors:
#       - "HighConfidence_Coop"
#     controls:
#       - "logMass_AVONET"
#       - "Centroid.Latitude_AVONET"
#     transformations:
#       Centroid.Latitude_AVONET: "abs"  # Use absolute latitude
#     complexity_levels:
#       - "main"
#       - "additive"
#       - "twoway"
#     max_interactions: 2
# 
# # Available variables by type:
# # Binary: HighConfidence_Coop, FemaleSong_Agg01, Griesser2017FamilialLiving, GeographicRegion_Jetz
# # Categorical: Territory, Migration_AVONET, TerritorialityWeakVsStrong
# # Continuous: logMass_AVONET, PercentAbsLogWingDimorphism, logMaleFemalePlumageDiffAbs, Centroid.Latitude_AVONET
# 
# # Available transformations:
# # abs: absolute value
# # log: logarithm
# # sqrt: square root
# # scale: standardize (center and scale)
# 
# # Complexity levels:
# # null: intercept only
# # main: single predictor models
# # additive: multiple predictors without interactions
# # twoway: two-way interactions
# # threeway: three-way interactions
# '
#   
#   writeLines(template, file)
#   message("Configuration template saved to:", file)
# }

#' Validate a configuration
#' 
#' @param config Single configuration object
#' @param data Data frame to check variables against
#' @return List with validation results
validate_config <- function(config, data = NULL) {
  
  issues <- character(0)
  warnings <- character(0)
  
  # Check required fields
  if (is.null(config$response)) {
    issues <- c(issues, "No response variable specified")
  }
  
  if (is.null(config$predictors) || length(config$predictors) == 0) {
    issues <- c(issues, "No predictor variables specified")
  }
  
  # Check for response in predictors
  if (config$response %in% config$predictors) {
    issues <- c(issues, "Response variable cannot be a predictor")
  }
  
  # Check complexity levels
  valid_levels <- c("null", "main", "additive", "twoway", "threeway")
  invalid_levels <- setdiff(config$complexity_levels, valid_levels)
  if (length(invalid_levels) > 0) {
    issues <- c(issues, paste("Invalid complexity levels:", 
                             paste(invalid_levels, collapse = ", ")))
  }
  
  # Check max_interactions
  if (config$max_interactions > length(config$predictors) + length(config$controls)) {
    warnings <- c(warnings, "max_interactions exceeds number of predictors")
  }
  
  # Check variables exist in data if provided
  if (!is.null(data)) {
    all_vars <- unique(c(config$response, config$predictors, config$controls))
    missing_vars <- setdiff(all_vars, colnames(data))
    if (length(missing_vars) > 0) {
      issues <- c(issues, paste("Variables not found in data:", 
                               paste(missing_vars, collapse = ", ")))
    }
  }
  
  return(list(
    valid = length(issues) == 0,
    issues = issues,
    warnings = warnings
  ))
}

#' Print method for configurations
#' 
#' @param x Configuration object
#' @param ... Additional arguments
print.phyloglm_config <- function(x, ...) {
  cat("PhyloGLM Analysis Configuration\n")
  cat("-------------------------------\n")
  if (!is.null(x$name)) cat("Name:", x$name, "\n")
  cat("Response:", x$response, "\n")
  cat("Predictors:", paste(x$predictors, collapse = ", "), "\n")
  if (!is.null(x$controls)) {
    cat("Controls:", paste(x$controls, collapse = ", "), "\n")
  }
  if (!is.null(x$transformations)) {
    cat("Transformations:\n")
    for (var in names(x$transformations)) {
      cat("  ", var, ":", x$transformations[[var]], "\n")
    }
  }
  cat("Complexity levels:", paste(x$complexity_levels, collapse = ", "), "\n")
  cat("Max interactions:", x$max_interactions, "\n")
}