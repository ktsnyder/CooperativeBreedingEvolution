# Enhanced version of create_bootstrap_forest_plot that includes model names
# This version extracts the base model name from the comprehensive table

create_bootstrap_forest_plot <- function(
  results_dir = "Outputs/PhyloglmResults/stepwise_boot500_20250716",
  output_file = "main_coefficients_comparison.png",
  aic_threshold = 2
) {
  
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
  library(grid)
  
  # Load the main results
  results_file <- list.files(results_dir, pattern = "expansion_results_boot.*\\.rds$", full.names = TRUE)[1]
  expansion_results <- readRDS(results_file)
  
  # Also load the main effects summary for easier access
  summary_file <- list.files(results_dir, pattern = "main_effects_changes_boot.*_with_ci\\.csv$", full.names = TRUE)[1]
  main_effects_summary <- read.csv(summary_file)
  
  # Load comprehensive table to get base model formulas
  comp_table_file <- list.files(results_dir, pattern = "stepwise_boot.*_comprehensive_table\\.csv$", full.names = TRUE)[1]
  comp_table <- read.csv(comp_table_file)
  
  # Extract unique base model names for each analysis
  base_model_info <- comp_table %>%
    select(Analysis_Name, Base_Model_AIC) %>%
    distinct() %>%
    group_by(Analysis_Name) %>%
    slice(1) %>%
    ungroup()
  
  # Try to get model names from the original RDS file
  model_names <- list()
  for (analysis in names(expansion_results)) {
    if (!is.null(expansion_results[[analysis]]$base_model)) {
      # The base_model object should have the formula
      formula_obj <- formula(expansion_results[[analysis]]$base_model)
      formula_str <- paste(deparse(formula_obj), collapse = " ")
      
      # Create a simplified model name from the formula
      # e.g., "FemaleSong_Agg01 ~ HighConfidence_Coop * TerritorialityWeakVsStrong + logMass_AVONET"
      # becomes something like "MainPred*Terr + Mass"
      model_name <- formula_str
      model_name <- gsub("FemaleSong_Agg01", "FS", model_name)
      model_name <- gsub("HighConfidence_Coop", "CB", model_name)
      model_name <- gsub("TerritorialityWeakVsStrong", "Terr", model_name)
      model_name <- gsub("Territory_12vs3", "Terr3", model_name)
      model_name <- gsub("logMass_AVONET", "Mass", model_name)
      model_name <- gsub(" ~ ", " ~ ", model_name)
      model_name <- gsub("\\s+", " ", model_name)  # Clean up extra spaces
      
      # Further simplification for display
      if (grepl("CB \\* Terr", model_name) || grepl("FS \\* Terr", model_name)) {
        if (grepl("Mass", model_name)) {
          model_name <- "MainPred×Terr + Mass"
        } else {
          model_name <- "MainPred×Terr"
        }
      } else if (grepl("CB \\+ Terr", model_name) || grepl("FS \\+ Terr", model_name)) {
        if (grepl("Mass", model_name)) {
          model_name <- "MainPred + Terr + Mass"
        } else {
          model_name <- "MainPred + Terr"
        }
      } else if (grepl("Mass", model_name)) {
        if (grepl("Terr", model_name)) {
          model_name <- "Terr + Mass"
        } else {
          model_name <- "MainPred + Mass"
        }
      } else {
        model_name <- "MainPred"
      }
      
      model_names[[analysis]] <- model_name
    }
  }
  
  # Prepare data for plotting
  plot_data <- main_effects_summary %>%
    filter(!is.na(Main_Effect_Base_Coef)) %>%
    mutate(
      Direction = ifelse(grepl("^FS_vs_CB", Analysis), "CB→FS", "FS→CB"),
      Analysis_Clean = case_when(
        Analysis == "FS_vs_CB_TerrWS_Mass" ~ "FS vs CB (Terr W/S)",
        Analysis == "CB_vs_FS_TerrWS_Mass" ~ "CB vs FS (Terr W/S)",
        Analysis == "FS_vs_CB_Terr3_Mass" ~ "FS vs CB (Terr 12vs3)",
        Analysis == "CB_vs_FS_Terr3_Mass" ~ "CB vs FS (Terr 12vs3)"
      ),
      Is_Significant = AIC_Improvement > aic_threshold,
      # Order analyses for consistent display
      Analysis_Order = factor(Analysis, levels = c("CB_vs_FS_TerrWS_Mass", "CB_vs_FS_Terr3_Mass",
                                                  "FS_vs_CB_TerrWS_Mass", "FS_vs_CB_Terr3_Mass"))
    ) %>%
    arrange(Analysis_Order, desc(AIC_Improvement))
  
  # Create individual plot function
  create_panel <- function(data, analysis_name, direction, model_name = NULL) {
    # Check if there's data to plot
    if (nrow(data) == 0) {
      panel_title <- ifelse(is.null(model_name), 
                           analysis_name, 
                           paste0(analysis_name, "\n(Base: ", model_name, ")"))
      return(
        ggplot() + 
          theme_void() +
          annotate("text", x = 0.5, y = 0.5, 
                  label = "No main coefficient\nto track in base model",
                  size = 4, hjust = 0.5, vjust = 0.5, color = "gray50") +
          labs(title = panel_title) +
          theme_minimal() +
          theme(
            plot.title = element_text(size = 10, face = "bold"),
            panel.border = element_rect(fill = NA, color = "gray80")
          )
      )
    }
    
    panel_title <- ifelse(is.null(model_name), 
                         analysis_name, 
                         paste0(analysis_name, "\n(Base: ", model_name, ")"))
    
    p <- ggplot(data, aes(x = reorder(Predictor_Name, AIC_Improvement))) +
      # Base model estimates
      geom_point(aes(y = Main_Effect_Base_Coef), shape = 1, size = 3, color = "gray40") +
      geom_errorbar(aes(ymin = Main_Effect_Base_CI_Lower, 
                       ymax = Main_Effect_Base_CI_Upper),
                   width = 0.2, alpha = 0.5, color = "gray40") +
      # Expanded model estimates  
      geom_point(aes(y = Main_Effect_Expanded_Coef), shape = 16, size = 3, color = "black") +
      geom_errorbar(aes(ymin = Main_Effect_Expanded_CI_Lower, 
                       ymax = Main_Effect_Expanded_CI_Upper),
                   width = 0.2, color = "black") +
      # Zero line
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
      # Styling
      coord_flip() +
      labs(x = "", 
           y = paste0("Coefficient (", direction, ")"),
           title = panel_title) +
      theme_minimal() +
      theme(
        axis.text.y = element_text(
          face = ifelse(data$Is_Significant, "bold", "plain"),
          size = 9
        ),
        plot.title = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 9),
        panel.border = element_rect(fill = NA, color = "gray80")
      )
    
    return(p)
  }
  
  # Create plots for each analysis
  plot_list <- list()
  
  # Order: CB_vs_FS_TerrWS, CB_vs_FS_Terr3, FS_vs_CB_TerrWS, FS_vs_CB_Terr3
  analysis_order <- c("CB_vs_FS_TerrWS_Mass", "CB_vs_FS_Terr3_Mass", 
                     "FS_vs_CB_TerrWS_Mass", "FS_vs_CB_Terr3_Mass")
  
  for (analysis in analysis_order) {
    analysis_data <- plot_data %>% filter(Analysis == analysis)
    
    # Get the base model name
    model_name <- model_names[[analysis]]
    
    # Get analysis clean name and direction
    if (analysis == "CB_vs_FS_TerrWS_Mass") {
      analysis_clean <- "CB vs FS (Terr W/S)"
      direction <- "FS→CB"
    } else if (analysis == "CB_vs_FS_Terr3_Mass") {
      analysis_clean <- "CB vs FS (Terr 12vs3)"
      direction <- "FS→CB"
    } else if (analysis == "FS_vs_CB_TerrWS_Mass") {
      analysis_clean <- "FS vs CB (Terr W/S)"
      direction <- "CB→FS"
    } else if (analysis == "FS_vs_CB_Terr3_Mass") {
      analysis_clean <- "FS vs CB (Terr 12vs3)"
      direction <- "CB→FS"
    }
    
    plot_list[[analysis]] <- create_panel(analysis_data, analysis_clean, direction, model_name)
  }
  
  # Arrange in 2x2 grid
  final_grid <- arrangeGrob(
    grobs = plot_list,
    ncol = 2,
    nrow = 2,
    top = textGrob("Main Association Coefficients: Base (○) vs Expanded (●) Models", 
                   gp = gpar(fontsize = 16, fontface = "bold")),
    bottom = textGrob("95% bootstrap confidence intervals shown. Bold parameters indicate AIC improvement > 2.", 
                     gp = gpar(fontsize = 10)),
    left = textGrob("Parameter added to expanded model", rot = 90, 
                   gp = gpar(fontsize = 12))
  )
  
  # Save plot
  ggsave(file.path(results_dir, output_file),
         final_grid, width = 14, height = 10, dpi = 300)
  
  message("Plot saved to: ", file.path(results_dir, output_file))
  
  # Also create a simpler version with patchwork if available
  if (requireNamespace("patchwork", quietly = TRUE)) {
    library(patchwork)
    
    # Create combined plot with patchwork
    combined_plot <- (plot_list[[1]] | plot_list[[2]]) / 
                     (plot_list[[3]] | plot_list[[4]]) +
      plot_annotation(
        title = "Main Association Coefficients: Base (○) vs Expanded (●) Models",
        subtitle = "95% bootstrap confidence intervals shown. Bold parameters indicate AIC improvement > 2.",
        caption = "Y-axis: Parameter added to expanded model"
      ) &
      theme(plot.title = element_text(size = 14, face = "bold"),
            plot.subtitle = element_text(size = 11))
    
    ggsave(file.path(results_dir, "main_coefficients_comparison_patchwork.png"),
           combined_plot, width = 14, height = 10, dpi = 300)
    
    message("Alternative plot saved to: ", file.path(results_dir, "main_coefficients_comparison_patchwork.png"))
  }
  
  return(plot_data)
}

# Test the function
# plot_data <- create_bootstrap_forest_plot(
#   results_dir = "Outputs/PhyloglmResults/stepwise_boot500_20250716",
#   output_file = "main_coefficients_comparison_with_models.png"
# )