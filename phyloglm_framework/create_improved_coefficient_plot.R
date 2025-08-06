# Standalone script to create improved coefficient comparison plot
library(ggplot2)
library(dplyr)
library(patchwork)

# Read the comprehensive table
comprehensive_table <- read.csv("Outputs/PhyloglmResults/stepwise_top2_models_20250715/top2_models_comprehensive_table1.csv")

# Read model deltaAIC info
all_results <- readRDS("Outputs/PhyloglmResults/combined_all_analyses_batch1-4_boot1000/all_results.rds")

# Build deltaAIC table
model_deltaAIC <- data.frame()
for (analysis_name in c("FS_vs_CB_TerrWS_Mass", "CB_vs_FS_TerrWS_Mass", 
                        "FS_vs_CB_Terr3_Mass", "CB_vs_FS_Terr3_Mass")) {
  model_comparison <- all_results[[analysis_name]]$comparison$comparison
  for (i in 1:2) {
    model_deltaAIC <- rbind(model_deltaAIC, data.frame(
      Analysis = analysis_name,
      Model_Rank = i,
      Model_Name = model_comparison$Model[i],
      deltaAIC = model_comparison$deltaAIC[i]
    ))
  }
}

# Prepare plot data
plot_data <- comprehensive_table %>%
  filter(!is.na(Base_Subset_Main_Coef)) %>%
  left_join(model_deltaAIC, by = c("Analysis", "Model_Rank")) %>%
  mutate(
    Direction = ifelse(grepl("^FS_vs_CB", Analysis), "CB→FS", "FS→CB"),
    Model_Label = paste0(Base_Model_Name, "\n(ΔAIC=", round(deltaAIC, 1), ")"),
    Is_Significant = AIC_Improvement > 2,
    Analysis_Clean = case_when(
      Analysis == "FS_vs_CB_TerrWS_Mass" ~ "FS vs CB (Terr W/S)",
      Analysis == "CB_vs_FS_TerrWS_Mass" ~ "CB vs FS (Terr W/S)",
      Analysis == "FS_vs_CB_Terr3_Mass" ~ "FS vs CB (Terr 1/3)",
      Analysis == "CB_vs_FS_Terr3_Mass" ~ "CB vs FS (Terr 1/3)"
    )
  )

# Create function to make individual plot
make_rank_plot <- function(data, rank, title) {
  p <- ggplot(data, aes(x = reorder(Predictor, Base_Subset_Main_Coef))) +
    # Base model estimates
    geom_point(aes(y = Base_Subset_Main_Coef), shape = 1, size = 3, color = "gray40") +
    geom_errorbar(aes(ymin = Base_Subset_Main_CI_Lower, 
                     ymax = Base_Subset_Main_CI_Upper),
                 width = 0.2, alpha = 0.5, color = "gray40") +
    # Expanded model estimates
    geom_point(aes(y = Expanded_Main_Coef), shape = 16, size = 3, color = "black") +
    geom_errorbar(aes(ymin = Expanded_Main_CI_Lower, 
                     ymax = Expanded_Main_CI_Upper),
                 width = 0.2, color = "black") +
    # Zero line
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    # Faceting
    facet_grid(Analysis_Clean ~ ., scales = "free_y", space = "free_y") +
    # Styling
    coord_flip() +
    labs(title = title, x = "", y = unique(data$Direction)) +
    theme_minimal() +
    theme(
      strip.text.y = element_text(size = 10, angle = 0),
      axis.text.y = element_text(
        face = ifelse(data$Is_Significant, "bold", "plain")
      ),
      panel.spacing = unit(0.5, "lines"),
      plot.title = element_text(size = 12, face = "bold")
    )
  
  # Add subtitle to panels
  p <- p + 
    ggtitle(title, subtitle = unique(data$Model_Label))
  
  return(p)
}

# Create plots for each rank
p_rank1 <- plot_data %>%
  filter(Model_Rank == 1) %>%
  make_rank_plot(rank = 1, title = "Best Models")

p_rank2 <- plot_data %>%
  filter(Model_Rank == 2) %>%
  make_rank_plot(rank = 2, title = "Second-Best Models")

# Add direction labels to x-axis
add_direction_labels <- function(p, direction_label) {
  p + scale_y_continuous(sec.axis = sec_axis(~., name = direction_label))
}

# Combine plots
final_plot <- (p_rank1 | p_rank2) +
  plot_annotation(
    title = "Main Association Coefficients: Base (○) vs Expanded (●) Models",
    subtitle = "95% confidence intervals shown. Bold parameters indicate AIC improvement > 2.",
    caption = "Y-axis: Parameter added to expanded model"
  ) &
  theme(plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11))

# Save plot
ggsave("Outputs/PhyloglmResults/stepwise_top2_models_20250715/main_coefficients_comparison_improved.png",
       final_plot, width = 14, height = 10, dpi = 300)

# Also create a simpler version with manual layout
library(gridExtra)
library(grid)

# Function to add boldface to significant predictors
make_label <- function(predictor, is_sig) {
  if (is_sig) {
    bquote(bold(.(predictor)))
  } else {
    predictor
  }
}

# Create individual plots arranged in correct order
# Order: CB_vs_FS_TerrWS, CB_vs_FS_Terr3, FS_vs_CB_TerrWS, FS_vs_CB_Terr3
analysis_order <- c("CB_vs_FS_TerrWS_Mass", "CB_vs_FS_Terr3_Mass", 
                   "FS_vs_CB_TerrWS_Mass", "FS_vs_CB_Terr3_Mass")

# Create a 4x2 grid of plots (8 positions total)
plot_grid <- list()

for (i in 1:8) {
  plot_grid[[i]] <- NULL
}

# Fill in the grid in the correct positions
for (analysis_idx in 1:4) {
  analysis_name <- analysis_order[analysis_idx]
  
  for (rank in 1:2) {
    # Calculate grid position (row-major order)
    grid_position <- (analysis_idx - 1) * 2 + rank
    
    # Check if this combination exists in the data
    analysis_data <- plot_data %>% 
      filter(Analysis == analysis_name & Model_Rank == rank)
    
    if (nrow(analysis_data) > 0) {
      # Create labels with boldface
      analysis_data$Predictor_Label <- mapply(
        make_label, 
        analysis_data$Predictor, 
        analysis_data$Is_Significant,
        SIMPLIFY = FALSE
      )
      
      p <- ggplot(analysis_data, aes(x = reorder(Predictor, AIC_Improvement))) +
        geom_point(aes(y = Base_Subset_Main_Coef), shape = 1, size = 3) +
        geom_errorbar(aes(ymin = Base_Subset_Main_CI_Lower, 
                         ymax = Base_Subset_Main_CI_Upper),
                     width = 0.2, alpha = 0.5) +
        geom_point(aes(y = Expanded_Main_Coef), shape = 16, size = 3) +
        geom_errorbar(aes(ymin = Expanded_Main_CI_Lower, 
                         ymax = Expanded_Main_CI_Upper),
                     width = 0.2) +
        geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
        coord_flip() +
        scale_x_discrete(labels = analysis_data$Predictor_Label) +
        labs(x = "", 
             y = paste0("Coefficient (", unique(analysis_data$Direction), ")"),
             title = unique(analysis_data$Analysis_Clean),
             subtitle = unique(analysis_data$Model_Label)) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 11),
          plot.subtitle = element_text(size = 9),
          axis.title.x = element_text(size = 9, color = "gray40")
        )
      
      plot_grid[[grid_position]] <- p
    } else {
      # Create empty plot for missing data
      empty_plot <- ggplot() + 
        theme_void() +
        annotate("text", x = 0.5, y = 0.5, 
                label = "No CB→FS coefficient\nto track in base model",
                size = 4, hjust = 0.5, vjust = 0.5, color = "gray50") +
        labs(title = ifelse(analysis_name == "FS_vs_CB_Terr3_Mass", 
                           "FS vs CB (Terr 1/3)",
                           unique(plot_data$Analysis_Clean[plot_data$Analysis == analysis_name])),
             subtitle = "Terr_Mass (ΔAIC=0.0)") +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 11),
          plot.subtitle = element_text(size = 9)
        )
      
      plot_grid[[grid_position]] <- empty_plot
    }
  }
}

# Arrange in grid with proper layout
final_grid <- arrangeGrob(
  grobs = plot_grid,
  ncol = 2,
  nrow = 4,
  top = textGrob("Main Association Coefficients: Base (○) vs Expanded (●) Models", 
                 gp = gpar(fontsize = 16, fontface = "bold")),
  bottom = textGrob("95% confidence intervals shown. Bold parameters indicate AIC improvement > 2.", 
                   gp = gpar(fontsize = 10)),
  left = textGrob("Parameter added to expanded model", rot = 90, 
                 gp = gpar(fontsize = 12))
)

ggsave("Outputs/PhyloglmResults/stepwise_top2_models_20250715/main_coefficients_comparison_grid.png",
       final_grid, width = 16, height = 12, dpi = 300)

cat("Improved plots saved to:\n")
cat("- main_coefficients_comparison_improved.png\n")
cat("- main_coefficients_comparison_grid.png\n")