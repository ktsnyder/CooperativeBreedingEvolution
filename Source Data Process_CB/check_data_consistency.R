# Add abs_Latitude, Migration_num, and logMass_normalized variables
df <- read.csv("Data_R_2025-06-09.csv")
df$Migration_num <- as.numeric(df$Migration_AVONET)-2
df$abs_Latitude <- abs(df$Centroid.Latitude_AVONET)
df$abs_Latitude_normalized <- scale(df$abs_Latitude, center = TRUE, scale = TRUE)
df$logMass_normalized <- scale(df$logMass_AVONET, center = TRUE, scale = TRUE)
df$PercentAbsLogWingDimorphism_normalized <- scale(df$PercentAbsLogWingDimorphism, center = TRUE, scale = TRUE)
df$logMaleFemalePlumageDiffAbs_normalized <- scale(df$logMaleFemalePlumageDiffAbs, center = TRUE, scale = TRUE)
#write.csv(df, "Data_R_2025-07-21.csv", row.names = FALSE)
data <- df

data = read.csv("Data_R_2025-07-21.csv")
# Look for the most recent Data_R file
data_files <- list.files(pattern = "^Data_R_.*\\.csv$", full.names = TRUE)
if (length(data_files) == 0) {
  stop("No Data_R_*.csv files found in current directory")
}
# Use the most recently modified file
latest_file <- data_files[which.max(file.mtime(data_files))]
cat("Using data file:", latest_file, "\n")
data <- read.csv(latest_file)

if (!file.exists("Data_R.csv")) {
  stop("Required file 'Data_R.csv' not found in current directory")
}
data_official <- read.csv("Data_R.csv")

#### Compare data ----## 
## -----------------------------------------------------------
## Compare two nearly-identical datasets: data vs data_official
## Assumptions:
## - Both have a "species" column with identical species names
## - data has 230 columns; data_official has those + 10 extras
## - Compare only the 230 common columns (including "species")
## -----------------------------------------------------------

## OPTIONAL: set a small numeric tolerance (e.g., 1e-9). 0 = exact.
tol <- 0.000001

## --- 1) Basic checks and alignment by species ----------------
stopifnot("species" %in% names(data), "species" %in% names(data_official))

# Ensure species are unique and identical as a set
if (anyDuplicated(data$species) || anyDuplicated(data_official$species)) {
  stop("Duplicate species detected; cannot align uniquely by 'species'.")
}
if (!setequal(data$species, data_official$species)) {
  missing_in_official <- setdiff(data$species, data_official$species)
  missing_in_data     <- setdiff(data_official$species, data$species)
  stop(sprintf(
    "Species sets differ.\nMissing in data_official (n=%d): %s\nMissing in data (n=%d): %s",
    length(missing_in_official), paste(head(missing_in_official, 10), collapse=", "),
    length(missing_in_data),     paste(head(missing_in_data, 10), collapse=", ")
  ))
}

# Align rows by species
o1 <- order(data$species)
o2 <- match(data$species[o1], data_official$species)
d1 <- data[o1, , drop = FALSE]
d2 <- data_official[o2, , drop = FALSE]
row.names(d1) <- row.names(d2) <- d1$species

# Common columns (by name), preserving data’s column order
common_cols <- intersect(names(d1), names(d2))

# If you want to exclude the key from the comparison, keep it out of `cmp_cols`
cmp_cols <- setdiff(common_cols, "species")

## --- 2) Helper to compare vectors with NA==NA and optional tol
compare_vec <- function(x, y, tol = 0) {
  # Coerce factors to their values
  if (is.factor(x)) x <- as.character(x)
  if (is.factor(y)) y <- as.character(y)
  
  both_na <- is.na(x) & is.na(y)
  
  # Numeric near-equality if both are numeric
  if (is.numeric(x) && is.numeric(y)) {
    both_num <- !is.na(x) & !is.na(y)
    near_eq  <- rep(FALSE, length(x))
    if (tol > 0) {
      near_eq[both_num] <- abs(x[both_num] - y[both_num]) <= tol
    } else {
      near_eq[both_num] <- x[both_num] == y[both_num]
    }
    out <- both_na | near_eq
  } else {
    # Fall back to character comparison
    out <- both_na | (as.character(x) == as.character(y))
  }
  out
}

## --- 3) Column class comparison (optional but helpful) --------
class_diff <- vapply(cmp_cols, function(cl) {
  paste(class(d1[[cl]]), collapse = "/") != paste(class(d2[[cl]]), collapse = "/")
}, logical(1))
class_mismatches <- data.frame(
  column = names(class_diff)[class_diff],
  class_in_data = vapply(names(class_diff)[class_diff], function(cl) paste(class(d1[[cl]]), collapse = "/"), character(1)),
  class_in_data_official = vapply(names(class_diff)[class_diff], function(cl) paste(class(d2[[cl]]), collapse = "/"), character(1)),
  row.names = NULL
)

## --- 4) Build a mismatch matrix (rows = species, cols = cmp_cols)
mismatch_mat <- sapply(cmp_cols, function(cl) {
  !compare_vec(d1[[cl]], d2[[cl]], tol = tol)
})
# Ensure matrix shape when one column
mismatch_mat <- as.matrix(mismatch_mat)
row.names(mismatch_mat) <- row.names(d1)

## --- 5) Summaries --------------------------------------------
col_mismatch_counts <- colSums(mismatch_mat, na.rm = TRUE)
row_mismatch_counts <- rowSums(mismatch_mat, na.rm = TRUE)
total_mismatches    <- sum(mismatch_mat, na.rm = TRUE)

## --- 6) Long-form differences (species, column, value1, value2)
if (total_mismatches > 0) {
  idx <- which(mismatch_mat, arr.ind = TRUE)
  diff_report <- data.frame(
    species = row.names(mismatch_mat)[idx[, "row"]],
    column  = colnames(mismatch_mat)[idx[, "col"]],
    value_in_data         = mapply(function(r, c) d1[[c]][r], idx[, "row"], colnames(mismatch_mat)[idx[, "col"]]),
    value_in_data_official= mapply(function(r, c) d2[[c]][r], idx[, "row"], colnames(mismatch_mat)[idx[, "col"]]),
    row.names = NULL,
    check.names = FALSE
  )
} else {
  diff_report <- data.frame(
    species = character(0),
    column  = character(0),
    value_in_data = character(0),
    value_in_data_official = character(0)
  )
}

## --- 7) Output/prints ----------------------------------------
cat("\n=== Dataset comparison (common columns only) ===\n")
cat(sprintf("Common columns compared (excluding 'species'): %d\n", length(cmp_cols)))
cat(sprintf("Total mismatches found: %d\n", total_mismatches))

if (nrow(class_mismatches) > 0) {
  cat("\nColumns with differing classes between data and data_official:\n")
  print(class_mismatches, row.names = FALSE)
}

if (total_mismatches == 0) {
  cat("\nAll compared cells match (NA==NA treated as equal).\n")
} else {
  cat("\nTop 20 columns by mismatch count:\n")
  print(sort(col_mismatch_counts, decreasing = TRUE)[1:min(20, length(col_mismatch_counts))])
  
  cat("\nTop 20 species (rows) by mismatch count:\n")
  print(sort(row_mismatch_counts, decreasing = TRUE)[1:min(20, length(row_mismatch_counts))])
  
  cat("\nFirst 20 detailed differences:\n")
  print(utils::head(diff_report, 20), row.names = FALSE)
}

## --- 8) (Optional) write CSV summaries -----------------------
# utils::write.csv(data.frame(column = names(col_mismatch_counts),
#                             mismatches = as.integer(col_mismatch_counts)),
#                  "column_mismatch_counts.csv", row.names = FALSE)
# utils::write.csv(data.frame(species = names(row_mismatch_counts),
#                             mismatches = as.integer(row_mismatch_counts)),
#                  "row_mismatch_counts.csv", row.names = FALSE)
# utils::write.csv(diff_report, "cell_level_differences.csv", row.names = FALSE)
# 
