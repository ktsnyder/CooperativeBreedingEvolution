### Run bias analyses and calculate downsampling ----

# Supplemental Table 13, 14 - Run bias analyses and calculate downsampling numbers ----
source(file.path("Bias_Test_functions","generate_bias_report.R"))
source(file.path("Bias_Test_functions","run_bias_tests.R"))
source(file.path("Bias_Test_functions","run_dimorphism_bias_tests.R"))
source(file.path("Bias_Test_functions","calculate_stratified_downsampling_with_territoriality.R"))
source(file.path("Phylopath_functions","run_phylopath_fxns.R") )

newdata = "Data_R.csv"
treefile = "ConsensusPasserineTreeHackett4_1000_OscineSubset.nex"
df_bias <- dfIn_phylo <- read.csv("Data_R.csv")
tree <- tree_phylo <- read.nexus("ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")
phylopath_output_dir = file.path("Outputs","PhylopathDownsampled", paste(all_traits_phylopath_label, "models"))
n_iterations = 10 # Set number of iterations (use 500 for publication-level results, 10 for testing)

# 1. Generate comprehensive bias report (all 6 tests)
cat("\nGenerating comprehensive bias report...\n")
generate_bias_report(
  df = df_bias,
  tree = tree,
  output_file = "Outputs/Bias_Test_Results.md",
  output_dir = "Outputs/"
)

# 2. Run standard bias tests for downsampling recommendations
cat("\nRunning bias tests for downsampling...\n")
bias_results <- run_bias_tests(
  df = df_bias,
  tree = tree,
  output_dir = "Outputs/BiasTests",
  geographic_col = "GeographicRegion_Jetz",
  territoriality_cols = c("TerritorialityWeakVsStrong", "Territory_12vs3"),
  save_plots = TRUE
)

# 3. Run dimorphism bias tests (optional - creates violin plots)
cat("\nRunning dimorphism bias tests...\n")
dimorphism_results <- run_dimorphism_bias_tests(
  df = df_bias,
  tree = tree,
  output_dir = file.path("Outputs","DimorphismBias")
)

# 4. Calculate stratified downsampling (for detailed bias correction)
cat("\nCalculating stratified downsampling for bias correction...\n")
downsampling_results <- calculate_stratified_downsampling(
  df = df_bias,
  stratify_vars = c("GeographicRegion_Jetz", "HighConfidence_Coop"),
  data_col = "FemaleSong_Agg01",
  output_file = file.path("Outputs", "Stratified_Downsampling_Calculations.md")
)

cat("Downsampling calculations saved to: Outputs/Stratified_Downsampling_Calculations.md\n")
cat("\nSummary of species to remove:\n")
if (!is.null(downsampling_results$downsampling$holarctic_noncoop)) {
  cat("- Holarctic non-cooperative:", downsampling_results$downsampling$holarctic_noncoop$n_to_remove, "species\n")
} # should be 83
if (!is.null(downsampling_results$downsampling$tropical_coop)) {
  cat("- Tropical cooperative:", downsampling_results$downsampling$tropical_coop$n_to_remove, "species\n")
} # should be 24
if (!is.null(downsampling_results$downsampling$global_coop)) {
  cat("- Global cooperative:", downsampling_results$downsampling$global_coop$n_to_remove, "species\n")
} # should be 15
if (!is.null(downsampling_results$downsampling$territoriality_bias)) {
  cat("- Territoriality bias:", downsampling_results$downsampling$territoriality_bias$calculation, "\n")
} # should be 266
if (!is.null(downsampling_results$downsampling$territory_12vs3_bias)) {
  cat("- Territory_12vs3 bias:", downsampling_results$downsampling$territory_12vs3_bias$calculation, "\n")
} # should be 155
