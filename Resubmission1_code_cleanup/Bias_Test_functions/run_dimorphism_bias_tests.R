# run_dimorphism_bias_tests.R
# Functions for testing sexual dimorphism biases with visualizations
# Kate Snyder
# Created: 2025-06-08
#
# This script tests whether sexually dimorphic species are more likely to have
# female song data than monomorphic species, using both plumage and size dimorphism

library(phytools)
library(dplyr)
library(ggplot2)
library(cowplot)

#' Test for dimorphism bias in data availability
#'
#' @param df Data frame with species data
#' @param tree Phylogenetic tree
#' @param dimorphism_col Column name for dimorphism measure
#' @param data_col Column to test for data availability (e.g., "FemaleSong_Agg01")
#' @param output_dir Directory for saving plots
#' @return List with test results and plots
test_dimorphism_bias <- function(df,
                                tree,
                                dimorphism_col,
                                data_col = "FemaleSong_Agg01",
                                output_dir = "Outputs/DimorphismBias") {
  
  # Create output directory if needed
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Filter to species in tree with dimorphism data
  df_work <- df[df$species %in% tree$tip.label & !is.na(df[[dimorphism_col]]), ]
  
  # Remove infinite values
  df_work <- df_work[is.finite(df_work[[dimorphism_col]]), ]
  
  # Add indicators for having data
  df_work$HaveFSData <- !is.na(df_work[[data_col]])
  df_work$HaveCBData <- !is.na(df_work$HighConfidence_Coop)
  df_work$HaveBothData <- df_work$HaveFSData & df_work$HaveCBData
  
  # Convert logical to factor for better plotting
  df_work$HaveFSData_factor <- factor(df_work$HaveFSData, 
                                     levels = c(FALSE, TRUE),
                                     labels = c("No FS Data", "Has FS Data"))
  
  df_work[[data_col]] <- factor(df_work[[data_col]],
                                levels = c(0, 1),
                                labels = c("FS Absent", "FS Present"))
  
  results <- list()
  
  # 1. Test with HaveFSData
  cat("\n=== Testing", dimorphism_col, "bias on having Female Song data ===\n")
  
  # Perform t-test
  t_test_have_data <- t.test(df_work[[dimorphism_col]] ~ df_work$HaveFSData)
  results$t_test_have_data <- t_test_have_data
  
  cat("Welch Two Sample t-test:\n")
  cat("t =", round(t_test_have_data$statistic, 3), 
      ", df =", round(t_test_have_data$parameter, 1),
      ", p-value =", format.pval(t_test_have_data$p.value, digits = 3), "\n")
  cat("Mean dimorphism (No FS Data):", round(t_test_have_data$estimate[1], 3), "\n")
  cat("Mean dimorphism (Has FS Data):", round(t_test_have_data$estimate[2], 3), "\n\n")
  
  # Perform Wilcoxon rank-sum test on the same data
  wilcox_test_have_data <- wilcox.test(df_work[[dimorphism_col]] ~ df_work$HaveFSData)
  results$wilcox_test_have_data <- wilcox_test_have_data
  
  cat("Wilcoxon rank sum test:\n")
  cat("W =", wilcox_test_have_data$statistic, 
      ", p-value =", format.pval(wilcox_test_have_data$p.value, digits = 3), "\n")
  cat("Median dimorphism (No FS Data):", round(median(df_work[[dimorphism_col]][df_work$HaveFSData == FALSE], na.rm = TRUE), 3), "\n")
  cat("Median dimorphism (Has FS Data):", round(median(df_work[[dimorphism_col]][df_work$HaveFSData == TRUE], na.rm = TRUE), 3), "\n\n")
  
  # Create violin plot for HaveFSData
  p1 <- ggplot(df_work, aes(x = HaveFSData_factor, y = .data[[dimorphism_col]], 
                           fill = HaveFSData_factor)) +
    geom_violin(alpha = 0.7, scale = "width") +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
    scale_fill_manual(values = c("No FS Data" = "#E41A1C", "Has FS Data" = "#377EB8")) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          plot.title = element_text(size = 16, face = "bold")) +
    labs(x = "Female Song Data Availability",
         y = dimorphism_col,
         title = paste("Sexual Dimorphism by Data Availability"),
         subtitle = paste0("t-test p = ", format.pval(t_test_have_data$p.value, digits = 3),
                          " | wilcox p = ",  format.pval(wilcox_test_have_data$p.value, digits = 3), "  (n = ", nrow(df_work), " species)"))
  
  # 2. Test with actual FemaleSong values (excluding NA)
  df_fs_only <- df_work[!is.na(df_work[[data_col]]), ]
  
  if (nrow(df_fs_only) > 10) {
    cat("\n=== Testing", dimorphism_col, "by Female Song presence/absence ===\n")
    
    t_test_fs_value <- t.test(df_fs_only[[dimorphism_col]] ~ df_fs_only[[data_col]])
    results$t_test_fs_value <- t_test_fs_value
    
    cat("Welch Two Sample t-test:\n")
    cat("t =", round(t_test_fs_value$statistic, 3), 
        ", df =", round(t_test_fs_value$parameter, 1),
        ", p-value =", format.pval(t_test_fs_value$p.value, digits = 3), "\n")
    cat("Mean dimorphism (FS Absent):", round(t_test_fs_value$estimate[1], 3), "\n")
    cat("Mean dimorphism (FS Present):", round(t_test_fs_value$estimate[2], 3), "\n\n")
    
    # Create violin plot for FemaleSong values
    p2 <- ggplot(df_fs_only, aes(x = .data[[data_col]], y = .data[[dimorphism_col]], 
                                fill = .data[[data_col]])) +
      geom_violin(alpha = 0.7, scale = "width") +
      geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
      stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
      scale_fill_manual(values = c("FS Absent" = "#FEE08B", "FS Present" = "#3288BD")) +
      theme_minimal() +
      theme(legend.position = "none",
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 14),
            plot.title = element_text(size = 16, face = "bold")) +
      labs(x = "Female Song",
           y = dimorphism_col,
           title = paste("Sexual Dimorphism by Female Song Presence"),
           subtitle = paste0("p = ", format.pval(t_test_fs_value$p.value, digits = 3),
                            " (n = ", nrow(df_fs_only), " species with FS data)"))
  } else {
    p2 <- NULL
  }
  
  # 3. Summary statistics
  summary_stats <- df_work %>%
    group_by(HaveFSData) %>%
    summarise(
      n = n(),
      mean_dimorphism = mean(.data[[dimorphism_col]], na.rm = TRUE),
      sd_dimorphism = sd(.data[[dimorphism_col]], na.rm = TRUE),
      median_dimorphism = median(.data[[dimorphism_col]], na.rm = TRUE),
      .groups = "drop"
    )
  
  results$summary_stats <- summary_stats
  
  cat("\nSummary statistics:\n")
  print(summary_stats)
  
  # Combine plots
  if (!is.null(p2)) {
    combined_plot <- plot_grid(p1, p2, ncol = 2, align = "h")
  } else {
    combined_plot <- p1
  }
  
  # Save plots
  plot_filename <- paste0(dimorphism_col, "_bias_plots")
  ggsave(file.path(output_dir, paste0(plot_filename, ".pdf")), 
         combined_plot, width = 12, height = 6)
  ggsave(file.path(output_dir, paste0(plot_filename, ".png")), 
         combined_plot, width = 12, height = 6, dpi = 150)
  
  results$plots <- list(p1 = p1, p2 = p2, combined = combined_plot)
  
  return(results)
}

#' Run comprehensive dimorphism bias tests
#'
#' @param df Data frame with species data
#' @param tree Phylogenetic tree
#' @param output_dir Directory for saving results
#' @return List of all test results
run_dimorphism_bias_tests <- function(df,
                                    tree,
                                    output_dir = "Outputs/DimorphismBias") {
  
  results <- list()
  
  # Test plumage dimorphism
  if ("logMaleFemalePlumageDiffAbs" %in% colnames(df)) {
    cat("\n==== TESTING PLUMAGE DIMORPHISM BIAS ====\n")
    results$plumage <- test_dimorphism_bias(
      df = df,
      tree = tree,
      dimorphism_col = "logMaleFemalePlumageDiffAbs",
      data_col = "FemaleSong_Agg01",
      output_dir = output_dir
    )
  }
  
  # Test size dimorphism
  if ("PercentAbsLogWingDimorphism" %in% colnames(df)) {
    cat("\n==== TESTING SIZE DIMORPHISM BIAS ====\n")
    results$size <- test_dimorphism_bias(
      df = df,
      tree = tree,
      dimorphism_col = "PercentAbsLogWingDimorphism",
      data_col = "FemaleSong_Agg01",
      output_dir = output_dir
    )
  }
  
  # Save full results
  saveRDS(results, file.path(output_dir, "dimorphism_bias_test_results.rds"))
  
  # Create summary report
  cat("\n\n=== DIMORPHISM BIAS TEST SUMMARY ===\n")
  
  if (!is.null(results$plumage)) {
    cat("\nPlumage Dimorphism:\n")
    cat("  Effect on having FS data: p =", 
        format.pval(results$plumage$t_test_have_data$p.value, digits = 3), "\n")
    if (!is.null(results$plumage$t_test_fs_value)) {
      cat("  Effect on FS presence (when data available): p =", 
          format.pval(results$plumage$t_test_fs_value$p.value, digits = 3), "\n")
    }
  }
  
  if (!is.null(results$size)) {
    cat("\nSize Dimorphism:\n")
    cat("  Effect on having FS data: p =", 
        format.pval(results$size$t_test_have_data$p.value, digits = 3), "\n")
    if (!is.null(results$size$t_test_fs_value)) {
      cat("  Effect on FS presence (when data available): p =", 
          format.pval(results$size$t_test_fs_value$p.value, digits = 3), "\n")
    }
  }
  
  return(results)
}

# Example usage:
if (FALSE) {
  # Load data and tree
  df <- read.csv("Data_R_2025-06-07.csv")
  tree <- read.nexus("2022-03-16ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")
  
  # Run dimorphism bias tests
  dimorphism_results <- run_dimorphism_bias_tests(df, tree)
}