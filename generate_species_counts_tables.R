# R script to generate species counts tables from Data_R.csv
# Reproduces tables from Supplement_CooperativeBreeding_NatEcoEvo_Revision1_species_counts_tables_only.pdf

# USER CONFIGURATION
# Optional label to append to all CSV file names (set to "" for no label)
label <- "Data_R_official"  # Example: label <- "_subset1" or label <- "_analysis2024"
label <- "Data_R_2025-09-11"

# Load required libraries
library(dplyr)

# Read the data
data <- read.csv("Data_R.csv", stringsAsFactors = FALSE)
data <- read.csv('/Users/kate/Desktop/CooperativeBreedingEvolution/Source Data Process_CB/Data_R_2025-09-11.csv', stringsAsFactors = F)

# Remove rows with missing species names
data <- data[!is.na(data$species) & data$species != "", ]

cat("Total species in dataset:", nrow(data), "\n\n")

# Helper function to create file names with optional label
create_filename <- function(base_name, label) {
  if (label == "") {
    return(paste0(base_name, ".csv"))
  } else {
    return(paste0(base_name, label, ".csv"))
  }
}

# ============================================================================
# Supplemental Table 1: Binary states of each sociality metric
# ============================================================================

cat("=== Supplemental Table 1: Binary states of each sociality metric ===\n")

# Initialize results dataframe for Table 1
table1_results <- data.frame(
  Variable = character(),
  State_0_Label = character(),
  State_0_Count = integer(),
  State_1_Label = character(),
  State_1_Count = integer(),
  stringsAsFactors = FALSE
)

# Function to create binary count table and store results
create_binary_table <- function(column_name, state_0_label, state_1_label) {
  if (column_name %in% names(data)) {
    counts <- table(data[[column_name]], useNA = "no")
    if (length(counts) >= 2) {
      cat(sprintf("%-40s State 0 (%s): %4d, State 1 (%s): %4d\n", 
                  column_name, state_0_label, counts[1], state_1_label, counts[2]))
      
      # Add to results dataframe
      table1_results <<- rbind(table1_results, data.frame(
        Variable = column_name,
        State_0_Label = state_0_label,
        State_0_Count = as.integer(counts[1]),
        State_1_Label = state_1_label,
        State_1_Count = as.integer(counts[2]),
        stringsAsFactors = FALSE
      ))
    } else {
      cat(sprintf("%-40s Not enough states\n", column_name))
    }
  } else {
    cat(sprintf("%-40s Column not found\n", column_name))
  }
}

# Griesser et al (2023) variables
create_binary_table("Griesser2023.Colonial01", "Non-colonial", "Colonial")
create_binary_table("Griesser2017FamilialLiving", "Non-familial living", "Familial living")
create_binary_table("Griesser2023.Asocial0VsSocial1", "Asocial", "Social")
create_binary_table("Griesser2023.GroupsLargerThanPair", "Groups pair or smaller", "Groups larger than pair")
create_binary_table("Griesser2023.LargestGroupSizes", "Smaller groups", "Largest group sizes")
create_binary_table("Griesser2023.LongSocialBonds", "Short social bonds", "Season or shorter social bond")
create_binary_table("Griesser2023.SeasonOrLongerSocialBonds", "Season or shorter social bonds", "Long social bonds")
create_binary_table("Griesser2023.MoreThanTwoCaretakers", "Two or fewer caretakers", "More than two caretakers")
create_binary_table("Griesser2023.TwoOrMoreCaretakers", "Fewer than two caretakers", "Two or more caretakers")

# Other variables
create_binary_table("Final.polygyny", "Social monogamy", "Social polygyny")
create_binary_table("Territory_12vs3", "No/weak/seasonal territoriality", "Year-round territoriality")
create_binary_table("TerritorialityWeakVsStrong", "Weak or no territoriality", "Territory actively defended")

# Cooperative breeding classifications
create_binary_table("HighConfidence_Coop", "Non-cooperative", "Cooperative")
create_binary_table("MeanCoopOmitTies", "Non-cooperative", "Cooperative")
create_binary_table("MeanCoopTie2Coop", "Non-cooperative", "Cooperative")
create_binary_table("MeanCoopTie2Noncoop", "Non-cooperative", "Cooperative")
create_binary_table("AnyCoopEqualsCoop", "Non-cooperative", "Cooperative")
create_binary_table("AnyNoncoopEqualsNoncoop", "Non-cooperative", "Cooperative")

# Save Table 1 to CSV
write.csv(table1_results, create_filename("supplemental_table1_binary_states", label), row.names = FALSE)
cat("Saved Supplemental Table 1 to:", create_filename("supplemental_table1_binary_states", label), "\n\n")

### ------------------------------------------------------------
## Supplemental Table 2: counts by sociality-state vs key traits
## ------------------------------------------------------------

## INPUT:
## - df: data.frame with the columns listed below
## - Values must be coded 0/1 (NAs allowed; NAs are ignored in counts)

## Example:
## df <- read.csv("your_species_traits.csv", stringsAsFactors = FALSE)
df <- data

## ---- configure trait columns ----
sociality_traits <- c(
  "Griesser2023.Colonial01",
  "Griesser2017FamilialLiving",
  "Griesser2023.Asocial0VsSocial1",
  "Griesser2023.GroupsLargerThanPair",
  "Griesser2023.LargestGroupSizes",
  "Griesser2023.LongSocialBonds",
  "Griesser2023.SeasonOrLongerSocialBonds",
  "Griesser2023.MoreThanTwoCaretakers",
  "Griesser2023.TwoOrMoreCaretakers",
  "Final.polygyny",
  "TerritorialityWeakVsStrong",
  "Territory_12vs3"
)

female_col <- "FemaleSong_Agg01"     # 0/1
coop_col    <- "HighConfidence_Coop" # 0/1

## ---- small helpers ----

# Safely coerce a vector to 0/1 (keeps NAs as NA)
bin_coerce <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  x <- suppressWarnings(as.numeric(x))
  x[!(x %in% c(0,1))] <- NA_real_
  x
}

# Count 0s and 1s in `col` within `idx` subset, ignoring NAs in the counted column
count_zeros_ones <- function(df, idx, col) {
  v <- bin_coerce(df[[col]])[idx]
  c(
    `0` = sum(v == 0, na.rm = TRUE),
    `1` = sum(v == 1, na.rm = TRUE)
  )
}

## ---- sanitize key columns once ----
df[[female_col]] <- bin_coerce(df[[female_col]])
df[[coop_col]]   <- bin_coerce(df[[coop_col]])

## ---- build the table ----
rows_list <- list()

for (tr in sociality_traits) {
  # ensure binary for the sociality column too (keeps NA as NA)
  df[[tr]] <- bin_coerce(df[[tr]])
  
  # states 0 and 1 for this sociality trait
  for (state in c(0, 1)) {
    # subset: rows where this trait == state (and not NA)
    idx <- which(df[[tr]] == state)
    
    # counts of Female Song within this subset
    fs_counts   <- count_zeros_ones(df, idx, female_col)
    # counts of Cooperative Breeding within this subset
    coop_counts <- count_zeros_ones(df, idx, coop_col)
    
    # assemble a row
    row_name <- paste0(tr, "_", state)
    rows_list[[row_name]] <- data.frame(
      Row = row_name,
      FemaleSong_Agg01_0   = unname(fs_counts["0"]),
      FemaleSong_Agg01_1   = unname(fs_counts["1"]),
      HighConfidence_Coop_0 = unname(coop_counts["0"]),
      HighConfidence_Coop_1 = unname(coop_counts["1"]),
      stringsAsFactors = FALSE
    )
  }
}

supp_table2 <- do.call(rbind, rows_list)
row.names(supp_table2) <- NULL

## ---- optional: order rows by trait then state (0 before 1) ----
order_idx <- order(
  match(sub("_\\d+$", "", supp_table2$Row), sociality_traits),
  sub("^.*_(\\d+)$", "\\1", supp_table2$Row)
)
supp_table2 <- supp_table2[order_idx, ]

## ---- show the result ----
print(supp_table2)

## ---- optional: write to CSV ----
# write.csv(supp_table2, "Supplemental_Table2_counts.csv", row.names = FALSE)


# Save Table 2 basic results to CSV
write.csv(supp_table2, 
          paste0("supplemental_table2_female_song_coop_sociality ", label, ".csv"), row.names = FALSE)
cat("Saved Supplemental Table 2 to:", paste0("supplemental_table2_female_song_coop_sociality ", label, ".csv"), "\n\n")

# ============================================================================
# Supplemental Table 9: Social system classifications and female song
# ============================================================================

cat("=== Supplemental Table 9: Social system classifications and female song ===\n")

# Initialize results dataframe for Table 9
table9_results <- data.frame(
  Classification = character(),
  Female_Song = character(),
  Cooperative_Breeding = character(),
  Count = integer(),
  stringsAsFactors = FALSE
)

coop_classifications <- c("HighConfidence_Coop", "MeanCoopTie2Noncoop", "MeanCoopTie2Coop", 
                         "MeanCoopOmitTies", "AnyCoopEqualsCoop", "AnyNoncoopEqualsNoncoop", "BiagoliniCoop", "DowningCoop", "HighConf_Coop_DefaultToCockburnInferred", "CockburnCoop", "CockburnInferred", "JetzCoopInclCockburn", "Griesser2017Coop", "DaleCoop", "CornwallisCoop")

for (coop_col in coop_classifications) {
  if (coop_col %in% names(data)) {
    cat(sprintf("\n%s:\n", coop_col))
    
    # Get the 2x2 table data
    filtered_data <- data[!is.na(data[[coop_col]]) & 
                         !is.na(data$FemaleSong_Agg01), ]
    
    if (nrow(filtered_data) > 0) {
      table_result <- table(filtered_data$FemaleSong_Agg01, filtered_data[[coop_col]])
      
      if (nrow(table_result) == 2 && ncol(table_result) == 2) {
        cat(sprintf("Female Song Absent, Non-cooperative: %d\n", table_result[1,1]))
        cat(sprintf("Female Song Absent, Cooperative: %d\n", table_result[1,2]))
        cat(sprintf("Female Song Present, Non-cooperative: %d\n", table_result[2,1]))
        cat(sprintf("Female Song Present, Cooperative: %d\n", table_result[2,2]))
        
        # Add to table9_results
        new_rows <- data.frame(
          Classification = rep(coop_col, 4),
          Female_Song = c("Absent", "Absent", "Present", "Present"),
          Cooperative_Breeding = c("Non-cooperative", "Cooperative", "Non-cooperative", "Cooperative"),
          Count = c(table_result[1,1], table_result[1,2], table_result[2,1], table_result[2,2]),
          stringsAsFactors = FALSE
        )
        table9_results <- rbind(table9_results, new_rows)
      }
    }
  }
}

# Save Table 9 to CSV
write.csv(table9_results, create_filename("supplemental_table9_social_system_classifications", label), row.names = FALSE)
cat("Saved Supplemental Table 9 to:", create_filename("supplemental_table9_social_system_classifications", label), "\n\n")

# ============================================================================
# Supplemental Table 12: Territorial contexts
# ============================================================================

cat("=== Supplemental Table 12: Territorial contexts ===\n")

# Initialize results dataframe for Table 12
table12_results <- data.frame(
  Territorial_Classification = character(),
  Territorial_State = character(),
  Female_Song = character(),
  Cooperative_Breeding = character(),
  Count = integer(),
  stringsAsFactors = FALSE
)

# Function for territorial analysis
create_territorial_table <- function(territorial_col) {
  if (territorial_col %in% names(data) && "FemaleSong_Agg01" %in% names(data) && 
      "HighConfidence_Coop" %in% names(data)) {
    
    # Remove NAs
    clean_data <- data[!is.na(data[[territorial_col]]) & 
                       !is.na(data$FemaleSong_Agg01) & 
                       !is.na(data$HighConfidence_Coop), ]
    
    if (nrow(clean_data) > 0) {
      territorial_levels <- unique(clean_data[[territorial_col]])
      
      for (terr_level in territorial_levels) {
        subset_data <- clean_data[clean_data[[territorial_col]] == terr_level, ]
        
        cat(sprintf("\n%s = %s:\n", territorial_col, terr_level))
        
        table_result <- table(subset_data$FemaleSong_Agg01, subset_data$HighConfidence_Coop)
        
        if (nrow(table_result) == 2 && ncol(table_result) == 2) {
          cat(sprintf("Female Song Absent, Non-cooperative: %d\n", table_result[1,1]))
          cat(sprintf("Female Song Absent, Cooperative: %d\n", table_result[1,2]))
          cat(sprintf("Female Song Present, Non-cooperative: %d\n", table_result[2,1]))
          cat(sprintf("Female Song Present, Cooperative: %d\n", table_result[2,2]))
          
          # Add to table12_results
          new_rows <- data.frame(
            Territorial_Classification = rep(territorial_col, 4),
            Territorial_State = rep(as.character(terr_level), 4),
            Female_Song = c("Absent", "Absent", "Present", "Present"),
            Cooperative_Breeding = c("Non-cooperative", "Cooperative", "Non-cooperative", "Cooperative"),
            Count = c(table_result[1,1], table_result[1,2], table_result[2,1], table_result[2,2]),
            stringsAsFactors = FALSE
          )
          table12_results <<- rbind(table12_results, new_rows)
        }
      }
    }
  }
}

# Territorial analyses
if ("TerritorialityWeakVsStrong" %in% names(data)) {
  cat("Weak vs Strong Territoriality:\n")
  create_territorial_table("TerritorialityWeakVsStrong")
}

if ("Territory_12vs3" %in% names(data)) {
  cat("\nYear-round Territoriality:\n")
  create_territorial_table("Territory_12vs3")
}

# Save Table 12 to CSV
write.csv(table12_results, create_filename("supplemental_table12_territorial_contexts", label), row.names = FALSE)
cat("Saved Supplemental Table 12 to:", create_filename("supplemental_table12_territorial_contexts", label), "\n\n")

# ============================================================================
# Supplemental Table 25: Geographic region analysis
# ============================================================================

cat("=== Supplemental Table 25: Geographic region analysis ===\n")

# Initialize results dataframe for Table 25
table25_results <- data.frame(
  Geographic_Region = character(),
  Female_Song = character(),
  Cooperative_Breeding = character(),
  Count = integer(),
  stringsAsFactors = FALSE
)

if ("GeographicRegion_Jetz" %in% names(data) || "Realm_Jetz2011" %in% names(data)) {
  # Try different geographic column names
  geo_col <- NULL
  if ("GeographicRegion_Jetz" %in% names(data)) {
    geo_col <- "GeographicRegion_Jetz"
  } else if ("Realm_Jetz2011" %in% names(data)) {
    geo_col <- "Realm_Jetz2011"
  }
  
  if (!is.null(geo_col)) {
    clean_data <- data[!is.na(data[[geo_col]]) & 
                       !is.na(data$FemaleSong_Agg01) & 
                       !is.na(data$HighConfidence_Coop), ]
    
    # Create simplified geographic categories (Holarctic vs Tropical)
    if (nrow(clean_data) > 0) {
      geo_regions <- unique(clean_data[[geo_col]])
      
      for (region in geo_regions) {
        subset_data <- clean_data[clean_data[[geo_col]] == region, ]
        
        cat(sprintf("\n%s:\n", region))
        
        table_result <- table(subset_data$FemaleSong_Agg01, subset_data$HighConfidence_Coop)
        
        if (nrow(table_result) == 2 && ncol(table_result) == 2) {
          cat(sprintf("Female Song Absent, Non-cooperative: %d\n", table_result[1,1]))
          cat(sprintf("Female Song Absent, Cooperative: %d\n", table_result[1,2]))
          cat(sprintf("Female Song Present, Non-cooperative: %d\n", table_result[2,1]))
          cat(sprintf("Female Song Present, Cooperative: %d\n", table_result[2,2]))
          
          # Add to table25_results
          new_rows <- data.frame(
            Geographic_Region = rep(region, 4),
            Female_Song = c("Absent", "Absent", "Present", "Present"),
            Cooperative_Breeding = c("Non-cooperative", "Cooperative", "Non-cooperative", "Cooperative"),
            Count = c(table_result[1,1], table_result[1,2], table_result[2,1], table_result[2,2]),
            stringsAsFactors = FALSE
          )
          table25_results <- rbind(table25_results, new_rows)
        }
      }
    }
  }
}

# Save Table 25 to CSV
write.csv(table25_results, create_filename("supplemental_table25_geographic_regions", label), row.names = FALSE)
cat("Saved Supplemental Table 25 to:", create_filename("supplemental_table25_geographic_regions", label), "\n\n")

# ============================================================================
# Summary statistics
# ============================================================================

cat("\n=== Summary Statistics ===\n")
cat("Total species:", nrow(data), "\n")

# Count non-missing values for key variables
key_vars <- c("HighConfidence_Coop", "FemaleSong_Agg01", "TerritorialityWeakVsStrong")
for (var in key_vars) {
  if (var %in% names(data)) {
    non_missing <- sum(!is.na(data[[var]]))
    cat(sprintf("%s: %d species with data\n", var, non_missing))
  }
}

cat("\nScript completed successfully!\n")
cat("All tables have been saved to CSV files with the label:", ifelse(label == "", "(no label)", label), "\n")


