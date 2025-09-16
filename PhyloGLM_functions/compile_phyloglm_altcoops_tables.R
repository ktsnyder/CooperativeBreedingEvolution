library(tidyverse)
library(flextable)

# Set the path to the data directory
data_dir <- file.path("Outputs", "PhyloGLM_outputs", "AltCoops")

# Option to use bootstrap means for odds ratios (default: FALSE)
use_bootstrap_OR <- TRUE  # Set to TRUE to recalculate OR from Boot_Mean

# Define custom order for cooperative breeding classifiers
coop_order <- c(
  "MeanCoopTie2Noncoop",
  "MeanCoopTie2Coop",
  "MeanCoopOmitTies",
  "AnyCoopEqualsCoop",
  "AnyNoncoopEqualsNoncoop",
  "CockburnCoop",
  "CockburnInferred",
  "HighConf_Coop_DefaultToCockburnInferred",
  "Griesser2017Coop",
  "DaleCoop",
  "CornwallisCoop",
  "JetzCoopInclCockburn"
)

# Get species counts
data = read.csv("Data_R.csv")
species_n <- list()
for (i in 1:length(coop_order)) {
  tempcoop = coop_order[i]
  datasub = data[complete.cases(data[,c(tempcoop, "FemaleSong_Agg01", "TerritorialityWeakVsStrong", "logMass_normalized")]),]
  species_n[[i]] <- c(tempcoop, length(datasub$species))
}

# Convert species_n list to data frame for easy joining
species_n_df <- do.call(rbind, species_n) %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  setNames(c("Classifier", "n_species")) %>%
  mutate(n_species = as.numeric(n_species))

# Function to extract response variable name from filename
extract_response_var <- function(filename) {
  # Extract the part between "coefficients_OddsRatios_" and "_vs_"
  str_extract(filename, "(?<=coefficients_OddsRatios_).*?(?=_vs_)")
}

# Function to extract cooperative breeding classifier from filename
extract_coop_classifier <- function(filename) {
  # Extract the part between "FS_vs_" and "xTerrWS"
  str_extract(filename, "(?<=FS_vs_)[^x]+(?=xTerrWS)")
}

# Function to read and process CSV files
process_csv_file <- function(filepath, type = "vs_FS") {
  df <- read_csv(filepath, show_col_types = FALSE)
  
  # Extract filename without path
  filename <- basename(filepath)
  
  if (type == "vs_FS") {
    # Extract response variable name for AltCoop as response
    response_var <- extract_response_var(filename)
    
    # Handle bootstrap OR if requested
    if (use_bootstrap_OR) {
      df <- df %>%
        mutate(Odds_Ratio = exp(Boot_Mean))
    }
    
    # Select relevant columns and rename
    result <- df %>%
      select(Parameter, Odds_Ratio, Significance) %>%
      mutate(Response_Variable = response_var) %>%
      pivot_wider(names_from = Parameter, 
                  values_from = c(Odds_Ratio, Significance),
                  names_sep = "_")
    
  } else {  # type == "FS_vs"
    # Extract cooperative breeding classifier
    coop_classifier <- extract_coop_classifier(filename)
    
    # Handle bootstrap OR if requested
    if (use_bootstrap_OR) {
      df <- df %>%
        mutate(Odds_Ratio = exp(Boot_Mean))
    }
    
    # Select relevant columns and rename
    result <- df %>%
      select(Parameter, Odds_Ratio, Significance) %>%
      mutate(Cooperative_Breeding_Classifier = coop_classifier) %>%
      pivot_wider(names_from = Parameter, 
                  values_from = c(Odds_Ratio, Significance),
                  names_sep = "_")
  }
  
  return(result)
}

# Process all "_vs_FS" files (AltCoop as response variable)
vs_fs_files <- list.files(data_dir, 
                          pattern = "coefficients_OddsRatios_.*_vs_FSxTerrWS_boot[0-9]+.csv$", 
                          full.names = TRUE)

# Exclude files that start with "FS_vs"
vs_fs_files <- vs_fs_files[!grepl("FS_vs", basename(vs_fs_files))]

# Process each file and combine
altcoop_response_data <- map_df(vs_fs_files, ~process_csv_file(.x, type = "vs_FS"))

# Create Table 1: AltCoop as response variable
# Reorder and rename columns for the final table
altcoop_table <- altcoop_response_data %>%
  select(
    Response_Variable,
    Intercept = `Odds_Ratio_(Intercept)`,
    Intercept_sig = `Significance_(Intercept)`,
    FemaleSong = Odds_Ratio_FemaleSong_Agg01,
    FemaleSong_sig = Significance_FemaleSong_Agg01,
    Terr = Odds_Ratio_TerritorialityWeakVsStrong,
    Terr_sig = Significance_TerritorialityWeakVsStrong,
    FSxTerrWS = `Odds_Ratio_FemaleSong_Agg01:TerritorialityWeakVsStrong`,
    FSxTerrWS_sig = `Significance_FemaleSong_Agg01:TerritorialityWeakVsStrong`
  ) %>%
  mutate(Response_Variable = factor(Response_Variable, levels = coop_order)) %>%
  arrange(Response_Variable) %>%
  mutate(Response_Variable = as.character(Response_Variable)) %>%
  # Add n_species by joining with species_n_df
  left_join(species_n_df, by = c("Response_Variable" = "Classifier"))

# Note: Mass is NA for these models as it's not included when CB is the response

# Process all "FS_vs_" files (Female Song as response variable)
fs_vs_files <- list.files(data_dir, 
                          pattern = "coefficients_OddsRatios_FS_vs_.*xTerrWS_Mass_boot[0-9]+.csv$", 
                          full.names = TRUE)

# Create Table 2: Female Song as response variable
# Process each file individually to handle different coefficient names
fs_table_list <- list()

for (file in fs_vs_files) {
  df <- read_csv(file, show_col_types = FALSE)
  coop_classifier <- extract_coop_classifier(basename(file))
  
  # Find the cooperative breeding coefficient (varies by model)
  # It's the parameter that's not Intercept, Territoriality, logMass, or an interaction
  cb_param <- df$Parameter[!grepl("\\(Intercept\\)|Territoriality|logMass|:", df$Parameter)]
  
  # Extract values for each parameter
  intercept_row <- df[df$Parameter == "(Intercept)", ]
  cb_row <- df[df$Parameter == cb_param, ]
  terr_row <- df[df$Parameter == "TerritorialityWeakVsStrong", ]
  # Check for either logMass_normalized or logMass_AVONET
  mass_param <- if("logMass_normalized" %in% df$Parameter) "logMass_normalized" else "logMass_AVONET"
  mass_row <- df[df$Parameter == mass_param, ]
  interaction_row <- df[grepl(":", df$Parameter), ]
  
  # Handle bootstrap OR if requested
  if (use_bootstrap_OR) {
    intercept_OR <- exp(intercept_row$Boot_Mean)
    cb_OR <- exp(cb_row$Boot_Mean)
    terr_OR <- exp(terr_row$Boot_Mean)
    mass_OR <- exp(mass_row$Boot_Mean)
    interaction_OR <- exp(interaction_row$Boot_Mean[1])
  } else {
    intercept_OR <- intercept_row$Odds_Ratio
    cb_OR <- cb_row$Odds_Ratio
    terr_OR <- terr_row$Odds_Ratio
    mass_OR <- mass_row$Odds_Ratio
    interaction_OR <- interaction_row$Odds_Ratio[1]
  }
  
  row_data <- data.frame(
    Cooperative_Breeding_Classifier = coop_classifier,
    Response_Variable = "Female Song",
    Intercept = intercept_OR,
    Intercept_sig = intercept_row$Significance,
    Cooperative_Breeding = cb_OR,
    CB_sig = cb_row$Significance,
    TerritorialityWeakVsStrong = terr_OR,
    Terr_sig = terr_row$Significance,
    logMass_normalized = mass_OR,
    Mass_sig = mass_row$Significance,
    CBxTerrWS = interaction_OR,
    CBxTerrWS_sig = interaction_row$Significance[1],
    stringsAsFactors = FALSE
  )
  
  fs_table_list[[length(fs_table_list) + 1]] <- row_data
}

fs_table <- bind_rows(fs_table_list) %>%
  mutate(Cooperative_Breeding_Classifier = factor(Cooperative_Breeding_Classifier, 
                                                  levels = coop_order)) %>%
  arrange(Cooperative_Breeding_Classifier) %>%
  mutate(Cooperative_Breeding_Classifier = as.character(Cooperative_Breeding_Classifier)) %>%
  # Add n_species by joining with species_n_df
  left_join(species_n_df, by = c("Cooperative_Breeding_Classifier" = "Classifier"))

# Function to format values with significance and color
format_with_significance <- function(value, sig, digits = 3) {
  if (is.na(value)) return("")
  
  # Format the value
  formatted_val <- format(round(value, digits), nsmall = digits)
  
  # Add significance stars
  if (!is.na(sig) && sig != "") {
    # Replace "." with "^" for p < 0.1
    sig <- gsub("\\.", "^", sig)
    formatted_val <- paste0(formatted_val, sig)
  }
  
  return(formatted_val)
}

# Create formatted versions of the tables
altcoop_formatted <- altcoop_table %>%
  rowwise() %>%
  mutate(
    Intercept = format_with_significance(Intercept, Intercept_sig),
    FemaleSong = format_with_significance(FemaleSong, FemaleSong_sig),
    Terr = format_with_significance(Terr, Terr_sig),
    FSxTerrWS = format_with_significance(FSxTerrWS, FSxTerrWS_sig)
  ) %>%
  ungroup() %>%
  select(Response_Variable, n_species, Intercept, FemaleSong, Terr, FSxTerrWS)

fs_formatted <- fs_table %>%
  rowwise() %>%
  mutate(
    Intercept = format_with_significance(Intercept, Intercept_sig),
    Cooperative_Breeding = format_with_significance(Cooperative_Breeding, CB_sig),
    TerritorialityWeakVsStrong = format_with_significance(TerritorialityWeakVsStrong, Terr_sig),
    logMass_normalized = format_with_significance(logMass_normalized, Mass_sig),
    CBxTerrWS = format_with_significance(CBxTerrWS, CBxTerrWS_sig)
  ) %>%
  ungroup() %>%
  select(Cooperative_Breeding_Classifier, Response_Variable, n_species, Intercept, 
         Cooperative_Breeding, TerritorialityWeakVsStrong, logMass_normalized, CBxTerrWS)

# Create filename suffix if using bootstrap OR
file_suffix <- if (use_bootstrap_OR) "_OR_from_Boot_Means" else ""

# Save as CSV files
write_csv(altcoop_formatted, file.path(data_dir, paste0("phyloglm_altcoop_response_table", file_suffix, "boot", nBoot,".csv")))
write_csv(fs_formatted, file.path(data_dir, paste0("phyloglm_female_song_response_table", file_suffix, "boot", nBoot, ".csv")))

# Create publication-ready tables with color coding using flextable
# Function to determine cell color based on OR value and significance
get_cell_color <- function(or_val, sig) {
  if (is.na(or_val) || is.na(sig) || sig == "") {
    return("white")
  }
  # Don't color marginally significant results (p < 0.1)
  if (sig == ".") {
    return("white")
  }
  if (or_val > 1) {
    return("#3498db")  # Blue for OR > 1
  } else {
    return("#e74c3c")  # Red for OR < 1
  }
}

# Create flextable for AltCoop as response
ft_altcoop <- flextable(altcoop_formatted) %>%
  theme_vanilla() %>%
  autofit() %>%
  align(align = "center", part = "all") %>%
  bold(part = "header")

# Apply conditional formatting
for (i in 1:nrow(altcoop_table)) {
  for (col in c("Intercept", "FemaleSong", "Terr", "FSxTerrWS")) {
    or_col <- col
    sig_col <- paste0(col, "_sig")
    
    if (or_col %in% names(altcoop_table) && sig_col %in% names(altcoop_table)) {
      or_val <- altcoop_table[[or_col]][i]
      sig_val <- altcoop_table[[sig_col]][i]
      
      if (!is.na(or_val) && !is.na(sig_val) && sig_val != "") {
        cell_color <- get_cell_color(or_val, sig_val)
        if (cell_color != "white") {
          ft_altcoop <- ft_altcoop %>%
            bg(i = i, j = col, bg = cell_color) %>%
            color(i = i, j = col, color = "white")
        }
      }
    }
  }
}

# Create flextable for Female Song as response
ft_fs <- flextable(fs_formatted) %>%
  theme_vanilla() %>%
  autofit() %>%
  align(align = "center", part = "all") %>%
  bold(part = "header")

# Apply conditional formatting
for (i in 1:nrow(fs_table)) {
  for (col in c("Intercept", "Cooperative_Breeding", "TerritorialityWeakVsStrong", 
                "logMass_normalized", "CBxTerrWS")) {
    # Map display column names to data column names
    data_col <- switch(col,
                       "Cooperative_Breeding" = "Cooperative_Breeding",
                       "TerritorialityWeakVsStrong" = "TerritorialityWeakVsStrong",
                       "logMass_normalized" = "logMass_normalized",
                       "CBxTerrWS" = "CBxTerrWS",
                       col)
    
    sig_col <- switch(col,
                      "Cooperative_Breeding" = "CB_sig",
                      "TerritorialityWeakVsStrong" = "Terr_sig",
                      "logMass_normalized" = "Mass_sig",
                      "CBxTerrWS" = "CBxTerrWS_sig",
                      paste0(col, "_sig"))
    
    if (data_col %in% names(fs_table) && sig_col %in% names(fs_table)) {
      or_val <- fs_table[[data_col]][i]
      sig_val <- fs_table[[sig_col]][i]
      
      if (!is.na(or_val) && !is.na(sig_val) && sig_val != "") {
        cell_color <- get_cell_color(or_val, sig_val)
        if (cell_color != "white") {
          ft_fs <- ft_fs %>%
            bg(i = i, j = col, bg = cell_color) %>%
            color(i = i, j = col, color = "white")
        }
      }
    }
  }
}

# Save flextables as HTML for viewing
save_as_html(ft_altcoop, path = file.path(data_dir, paste0("phyloglm_altcoop_response_table", file_suffix, "boot", nBoot, ".html")))
save_as_html(ft_fs, path = file.path(data_dir, paste0("phyloglm_female_song_response_table", file_suffix, "boot", nBoot, ".html")))

# Print tables to console
cat("\n=== Table 1: PhyloGLM AltCoops - AltCoop as Response Variable ===\n")
print(altcoop_formatted)

cat("\n\n=== Table 2: PhyloGLM AltCoops - Female Song as Response Variable ===\n")
print(fs_formatted)

cat("\n\nTables saved as:\n")
cat(paste0("- phyloglm_altcoop_response_table", file_suffix, "boot", nBoot,".csv\n"))
cat(paste0("- phyloglm_female_song_response_table", file_suffix, "boot", nBoot,".csv\n"))
cat(paste0("- phyloglm_altcoop_response_table", file_suffix,"boot", nBoot, ".html (with color coding)\n"))
cat(paste0("- phyloglm_female_song_response_table", file_suffix,"boot", nBoot, ".html (with color coding)\n"))

cat("\nColor coding: Blue = OR > 1, Red = OR < 1 (only for significant values p<0.05)\n")
cat("Significance: ^ p<0.1 (not colored), * p<0.05, ** p<0.01, *** p<0.001\n")

if (use_bootstrap_OR) {
  cat("\nNOTE: Odds ratios calculated from bootstrap means (Boot_Mean)\n")
} else {
  cat("\nNOTE: Odds ratios from original model estimates\n")
}

