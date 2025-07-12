### Repeatedly downsample species to account for biases in data
# Kate Snyder
# May 15, 2025
# Edited 6/9/2025 -  added downsample_dimorphism_bias()
# Also edited 6/9/2025 - added functions per Update_Run_Analyses_Instructions.md
# 6/10/2025 - new functions run_phylopath_dimorphism_correction(), create_downsampled_plots_flexible() per claude_code_sessions/phylopath_dimorphism_improvements.R
# 7/11/2025 - run_phylopath_dimorphism_correction(): added args to change which vars used, output to folder with trait names in folder name; moved actual "output_dir" and "phylopath_output_dir" variable creation to inside run_multiple_phylopath() and run_phylopath_dimorphism_correction(), respectively - args will no longer do anything.

library(phylopath)
library(phytools)


#### Main Phylopath Runner Functions ----

# Function to run multiple iterations and aggregate results
run_multiple_phylopath <- function(dfIn, tree, downsample_columns, downsample_values, numToRemove,
                                   n_iterations = 500, 
                                   female_song_var = "FemaleSong_Agg01", 
                                   coop_breeding_var = "HighConfidence_Coop", 
                                   territoriality_var = "TerritorialityWeakVsStrong", 
                                   mass_var = "logMass_AVONET",
                                   save_conditional_plots = TRUE,
                                   save_path_coefficients = TRUE,
                                   save_downsampling_plots = TRUE,
                                   plotlabel = "",
                                   output_dir = "obsolete" # now defaults to Outputs/PhylopathDownsampled/[trait names]
                                   ) {
  
  require(dplyr)
  require(ggplot2)
  require(gridExtra)
  
  all_traits_phylopath_label <- paste(female_song_var, coop_breeding_var, territoriality_var, mass_var)
  
  output_dir = file.path("Outputs","PhylopathDownsampled", paste(all_traits_phylopath_label, "models"))
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Create a data frame to store results
  results_df <- data.frame(
    seed = integer(n_iterations),
    nSub2dCIC_Models = integer(n_iterations),
    CICsub2_models_string = character(n_iterations),
    nSpecies = integer(n_iterations),
    best_model = character(n_iterations),
    best_CICc = numeric(n_iterations),
    stringsAsFactors = FALSE
  )
  
  # Create lists to store additional results
  conditional_average_plots <- list()
  path_coefficients <- list()
  all_models_frequency <- list()
  
  # For the detailed model dataframe
  detailed_models_list <- list()
  
  # Start time for progress tracking
  start_time <- Sys.time()
  
  # Run iterations
  for (i in 1:n_iterations) {
    # Set seed for this iteration
    current_seed <- i + 1000  # Using a base offset to avoid low seeds
    
    # Run the phylopath analysis
    cat(sprintf("Running iteration %d/%d (seed: %d)...\n", i, n_iterations, current_seed))
    current_output <- downsample_run_phylopath(
      dfIn = dfIn, 
      tree = tree, 
      downsample_columns = downsample_columns, 
      downsample_values = downsample_values, 
      numToRemove = numToRemove, 
      seed = current_seed,
      female_song_var = female_song_var, 
      coop_breeding_var = coop_breeding_var, 
      territoriality_var = territoriality_var, 
      mass_var = mass_var
    )
    
    # Extract summary from the result
    s <- summary(current_output$result)
    
    # Store basic results
    results_df$seed[i] <- current_seed
    results_df$nSub2dCIC_Models[i] <- length(current_output$CICsub2_models)
    results_df$CICsub2_models_string[i] <- paste(current_output$CICsub2_models, collapse = ", ")
    results_df$nSpecies[i] <- current_output$nSpecies
    results_df$best_model[i] <- s$model[1]  # First model is the best
    results_df$best_CICc[i] <- s$CICc[1]    # CICc of the best model
    
    # Store the conditional average plot if requested
    if (save_conditional_plots) {
      conditional_average_plots[[i]] <- current_output$conditional_average_plot + 
        ggtitle(paste("Seed:", current_seed, "| nSpecies:", current_output$nSpecies)) +
        theme(plot.title = element_text(size = 10, face = "bold"))
    }
    
    # Create detailed model information for each model with delta_CICc < 2
    sub2_models <- current_output$CICsub2_models
    for (model_idx in which(s$delta_CICc < 2)) {
      model_name <- s$model[model_idx]
      model_CICc <- s$CICc[model_idx]
      model_delta_CICc <- s$delta_CICc[model_idx]
      
      # Get the chosen model to extract path coefficients
      chosen_model <- choice(current_output$result, model_name)
      
      # Create a row for this model
      model_row <- data.frame(
        seed = current_seed,
        model = model_name,
        CICc = model_CICc,
        delta_CICc = model_delta_CICc,
        nSpecies = current_output$nSpecies,
        stringsAsFactors = FALSE
      )
      
      # Extract edge information properly from the coefficient matrix
      edges_coef <- chosen_model$coef
      edges_se <- chosen_model$se
      
      # Check if we have the coefficient matrix
      if (!is.null(edges_coef) && is.matrix(edges_coef)) {
        # Process coefficient matrix which is in adjacency matrix form
        for (from_idx in 1:nrow(edges_coef)) {
          for (to_idx in 1:ncol(edges_coef)) {
            # Get the coefficient value
            coef_value <- edges_coef[from_idx, to_idx]
            
            # Only process non-zero coefficients (actual paths)
            if (coef_value != 0) {
              # Get the variable names
              from_name <- rownames(edges_coef)[from_idx]
              to_name <- colnames(edges_coef)[to_idx]
              
              # Get standard error
              se_value <- edges_se[from_idx, to_idx]
              
              # Calculate approximate p-value (z-test)
              # Note: This is an approximation, might want to extract actual p-values if available
              z_value <- coef_value / se_value
              p_value <- 2 * (1 - pnorm(abs(z_value)))
              
              # Create column name for this edge
              edge_name <- paste0(from_name, "_to_", to_name)
              edge_est_name <- paste0(edge_name, "_est")
              edge_se_name <- paste0(edge_name, "_se")
              edge_p_name <- paste0(edge_name, "_p")
              
              # Add to model row
              model_row[[edge_est_name]] <- coef_value
              model_row[[edge_se_name]] <- se_value
              model_row[[edge_p_name]] <- p_value
            }
          }
        }
      }
      
      # Add to detailed models list
      detailed_models_list[[length(detailed_models_list) + 1]] <- model_row
    }
    
    # Extract and store path coefficients if requested
    if (save_path_coefficients) {
      # For each model, extract the path coefficients
      for (model_name in sub2_models) {
        # Extract the d_sep object for the model
        d_sep_obj <- current_output$result$d_sep[[model_name]]
        
        if (!is.null(d_sep_obj) && nrow(d_sep_obj) > 0) {
          # Create a simplified data frame without the complex list column
          d_sep_simple <- data.frame(
            d_sep = d_sep_obj$d_sep,
            p_value = d_sep_obj$p,
            seed = current_seed,
            model_name = model_name,
            stringsAsFactors = FALSE
          )
          
          # Add to the path_coefficients list
          path_coefficients[[length(path_coefficients) + 1]] <- d_sep_simple
        }
      }
      
      # Update frequency of all models
      for (model_name in s$model) {
        if (model_name %in% names(all_models_frequency)) {
          all_models_frequency[[model_name]] <- all_models_frequency[[model_name]] + 1
        } else {
          all_models_frequency[[model_name]] <- 1
        }
      }
    }
    
    # Progress update
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    estimated_total <- elapsed * (n_iterations / i)
    remaining <- estimated_total - elapsed
    cat(sprintf("  Progress: %.1f%% (estimated time remaining: %.1f minutes)\n", 
                100 * i / n_iterations, remaining / 60))
  }
  
  # Combine detailed models into a single data frame
  detailed_models_df <- bind_rows(detailed_models_list)
  
  # Save the detailed models data frame
  if (!is.null(detailed_models_df) && nrow(detailed_models_df) > 0) {
    write.csv(detailed_models_df, file.path(output_dir, paste0("detailed_models_", plotlabel, "_", n_iterations, "_", Sys.Date(), ".csv")), row.names = FALSE)
    cat("Detailed model information saved to CSV file\n")
  }
  
  # Convert all_models_frequency to a data frame
  model_freq_df <- data.frame(
    model_name = names(all_models_frequency),
    frequency = unlist(all_models_frequency),
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(frequency))
  
  # Combine all path coefficients into a single data frame if any were collected
  if (save_path_coefficients && length(path_coefficients) > 0) {
    path_coefficients_df <- do.call(rbind, path_coefficients)
    
    # Analyze path significance and directions
    path_summary <- path_coefficients_df %>%
      group_by(d_sep) %>%
      summarize(
        n_occurrences = n(),
        n_significant = sum(p_value < 0.05, na.rm = TRUE),
        prop_significant = n_significant / n_occurrences,
        mean_p = mean(p_value, na.rm = TRUE),
        median_p = median(p_value, na.rm = TRUE)
      ) %>%
      arrange(mean_p)
    
    # Save path coefficients to CSV
   # write.csv(path_coefficients_df, paste0("path_coefficients_", n_iterations, "_",Sys.Date(), ".csv"), row.names = FALSE)
  #  write.csv(path_summary, paste0("path_summary_", n_iterations, "_", Sys.Date(), ".csv"), row.names = FALSE)
    
   # cat("Path coefficients and summary saved to CSV files\n")
  } else {
    path_coefficients_df <- NULL
    path_summary <- NULL
  }
  
  # Save conditional average plots to a PDF if requested
  if (save_conditional_plots && length(conditional_average_plots) > 0) {
    pdf_name <- file.path(output_dir, paste0("conditional_average_plots_", plotlabel, "_", n_iterations, "_", Sys.Date(), ".pdf"))
    
    # Increased height from 10 to 15
    pdf(pdf_name, width = 12, height = 15)
    
    # Plot multiple plots per page
    plots_per_page <- 6
    num_plots <- length(conditional_average_plots)
    num_pages <- ceiling(num_plots / plots_per_page)
    
    for (page in 1:num_pages) {
      start_idx <- (page - 1) * plots_per_page + 1
      end_idx <- min(page * plots_per_page, num_plots)
      page_plots <- conditional_average_plots[start_idx:end_idx]
      
      # Arrange plots in a grid
      do.call(grid.arrange, c(page_plots, ncol = 2))
    }
    
    dev.off()
    cat(sprintf("Conditional average plots saved to %s\n", pdf_name))
  }
  
  # Calculate and add the frequency of each model appearing in the CICsub2_models
  all_sub2_models <- unlist(strsplit(results_df$CICsub2_models_string, ", "))
  if (length(all_sub2_models) > 0) {
    sub2_model_freq <- as.data.frame(table(all_sub2_models))
    names(sub2_model_freq) <- c("model_name", "frequency_in_sub2")
    sub2_model_freq$proportion_in_sub2 <- sub2_model_freq$frequency_in_sub2 / n_iterations
    
    # Merge with overall model frequency
    model_freq_df <- merge(model_freq_df, sub2_model_freq, by = "model_name", all = TRUE)
    model_freq_df[is.na(model_freq_df)] <- 0
    model_freq_df <- model_freq_df %>% arrange(desc(frequency_in_sub2))
  }
  
  # Save model frequencies to CSV
  write.csv(model_freq_df, file.path(output_dir, paste0("model_frequencies_", plotlabel, "_", n_iterations, "_", Sys.Date(), ".csv")), row.names = FALSE)
  
  # Return the results
  results <- list(
    results_df = results_df,
    model_frequencies = model_freq_df,
    path_summary = path_summary,
    path_coefficients = path_coefficients_df,
    detailed_models = detailed_models_df
  )
  
  # Create summary plot for top models if there are any
  if (nrow(model_freq_df) > 0 && "frequency_in_sub2" %in% names(model_freq_df)) {
    top_n <- min(10, nrow(model_freq_df))
    top_models <- model_freq_df %>% 
      arrange(desc(frequency_in_sub2)) %>% 
      head(top_n)
    
    # Create a bar plot of top models
    top_models_plot <- ggplot(top_models, aes(x = reorder(model_name, frequency_in_sub2), y = frequency_in_sub2)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      coord_flip() +
      labs(
        title = paste("Top", top_n, "Models by Frequency in CICsub2"),
        x = "Model",
        y = "Frequency"
      ) +
      theme_minimal()
    
    # Save the plot
    pdf(file.path(output_dir, paste0("top_models_", plotlabel, "_", n_iterations, "_", Sys.Date(), ".pdf")), width = 10, height = 8)
    print(top_models_plot)
    dev.off()
  }
  
  # Analyze edge coefficients across models
  if (!is.null(detailed_models_df) && nrow(detailed_models_df) > 0) {
    # Get all column names that end with "_est"
    est_cols <- grep("_est$", names(detailed_models_df), value = TRUE)
    
    # Create a summary of edge estimates
    edge_summary <- data.frame(
      edge = gsub("_est$", "", est_cols),
      stringsAsFactors = FALSE
    )
    
    for (col in est_cols) {
      # Get the corresponding p-value column
      p_col <- gsub("_est$", "_p", col)
      
      # Only use cases where the column exists
      if (p_col %in% names(detailed_models_df)) {
        # Extract data without NA values
        valid_data <- detailed_models_df[!is.na(detailed_models_df[[col]]), ]
        
        if (nrow(valid_data) > 0) {
          edge_index <- which(edge_summary$edge == gsub("_est$", "", col))
          
          # Calculate summary statistics
          edge_summary$n_occurrences[edge_index] <- nrow(valid_data)
          edge_summary$mean_estimate[edge_index] <- mean(valid_data[[col]], na.rm = TRUE)
          edge_summary$median_estimate[edge_index] <- median(valid_data[[col]], na.rm = TRUE)
          edge_summary$n_significant[edge_index] <- sum(valid_data[[p_col]] < 0.05, na.rm = TRUE)
          edge_summary$prop_significant[edge_index] <- edge_summary$n_significant[edge_index] / edge_summary$n_occurrences[edge_index]
          edge_summary$mean_p[edge_index] <- mean(valid_data[[p_col]], na.rm = TRUE)
        }
      }
    }
    
    # Sort by significance proportion
    edge_summary <- edge_summary %>% 
      arrange(desc(prop_significant))
    
    # Save edge summary
    write.csv(edge_summary, file.path(output_dir, paste0("edge_summary_", plotlabel, "_", n_iterations, "_", Sys.Date(), ".csv")), row.names = FALSE)
    
    # Add to results
    results$edge_summary <- edge_summary
  }
  
  return(results)
}



# female_song_var <- "FemaleSong_Agg01"
# coop_breeding_var <- "HighConfidence_Coop"
# territoriality_var <- "TerritorialityWeakVsStrong"
# mass_var <- "logMass_AVONET"

## example usage run_CB_FS_Terr_phylopath():
#output <- run_CB_FS_Terr_phylopath(dfIn = df_subset, female_song_var = "FemaleSong_Agg01", coop_breeding_var = "HighConfidence_Coop", territoriality_var = "TerritorialityWeakVsStrong", mass_var = "logMass_AVONET", plots2pdf = TRUE)

run_CB_FS_Terr_phylopath <- function(dfIn, tree, female_song_var, coop_breeding_var, territoriality_var, mass_var = "none", include_terr_response = FALSE, plots2pdf = FALSE,                                      plots2png = FALSE, output_dir = "Outputs/PhylopathPlots") {
  require(phylopath)
  require(dplyr)
  require(ggplot2)
  
  if (mass_var != "none") {
    dfIn_clean <- dfIn[complete.cases(dfIn[, c(female_song_var, coop_breeding_var, territoriality_var, mass_var)]), ]
  } else {
    dfIn_clean <- dfIn[complete.cases(dfIn[, c(female_song_var, coop_breeding_var, territoriality_var)]), ]  
  }
  
  tree_clean <- drop.tip(tree, setdiff(tree$tip.label, dfIn_clean$species))
  
  dfIn_clean[,female_song_var] <- as.factor(dfIn_clean[,female_song_var])
  dfIn_clean[,coop_breeding_var] <- as.factor(dfIn_clean[,coop_breeding_var])
  dfIn_clean[,territoriality_var] <- as.factor(dfIn_clean[,territoriality_var])
  nSpecies = length(tree_clean$tip.label)
  rownames(dfIn_clean) <- dfIn_clean$species
  
  # Create a variable mapping dictionary
  var_map <- list(
    "FS" = female_song_var,
    "COOP" = coop_breeding_var,
    "TERR" = territoriality_var,
    "MASS" = mass_var
  )
  
  
  #### Model sets ----
  model_patterns <- list(
    ## SERIES A: Basic models without mass, no direct FS-COOP interaction
    "A0_TERR→FS" = c("FS ~ TERR"),
    "A0_TERR→COOP" = c("COOP ~ TERR"),
    "A0_TERR→FS_TERR→COOP" = c("FS ~ TERR", "COOP ~ TERR"),
    
    ## SERIES N: Null models
    "N0_INDEP" = c("FS ~ 1", "COOP ~ 1", "TERR ~ 1"),
    
    ## SERIES X: COOP influences FS (direction 1)
    "X1_COOP→FS" = c("FS ~ COOP"),
    "X1_COOP→FS_TERR→FS" = c("FS ~ TERR + COOP"),
    "X1_COOP→FS_TERR→COOP" = c("FS ~ COOP", "COOP ~ TERR"),
    
    ## SERIES Y: FS influences COOP (direction 2)
    "Y2_FS→COOP" = c("COOP ~ FS"),
    "Y2_FS→COOP_TERR→COOP" = c("COOP ~ TERR + FS"),
    "Y2_TERR→FS_FS→COOP" = c("FS ~ TERR", "COOP ~ FS")
  )
  
  model_patterns_with_MASS <- list(
    ## SERIES B: Direct mass effects, no direct FS-COOP interaction
    "B0_TERR→FS_MASS→FS" = c("FS ~ TERR + MASS"),
    "B0_TERR→COOP_MASS→COOP" = c("COOP ~ TERR + MASS"),
    "B0_TERR→FS_TERR→COOP_MASS→FS_MASS→COOP" = c("FS ~ TERR + MASS", "COOP ~ TERR + MASS"),
    "B0_MASS→FS" = c("FS ~ MASS"),
    "B0_MASS→COOP" = c("COOP ~ MASS"),
    "B0_MASS→FS_MASS→COOP" = c("FS ~ MASS", "COOP ~ MASS"),
    
    ## SERIES C: Mass affects only COOP, no direct FS-COOP interaction
    "C0_TERR→FS_TERR→COOP_MASS→COOP" = c("FS ~ TERR", "COOP ~ TERR + MASS"),
    "C0_TERR→FS_MASS→COOP" = c("FS ~ TERR", "COOP ~ MASS"),
    
    ## SERIES D: Mass affects only FS, no direct FS-COOP interaction
    "D0_TERR→FS_TERR→COOP_MASS→FS" = c("FS ~ TERR + MASS", "COOP ~ TERR"),
    "D0_MASS→FS_TERR→COOP" = c("FS ~ MASS", "COOP ~ TERR"),
    
    ## SERIES E: Mass affects TERR which affects others, no direct FS-COOP interaction
    "E0_TERR→FS_TERR→COOP_MASS→TERR" = c("FS ~ TERR", "COOP ~ TERR", "TERR ~ MASS"),
    "E0_TERR→FS_TERR→COOP_MASS→TERR_MASS→FS" = c("FS ~ TERR + MASS", "COOP ~ TERR", "TERR ~ MASS"),
    "E0_TERR→FS_MASS→TERR" = c("FS ~ TERR", "TERR ~ MASS"),
    "E0_TERR→COOP_MASS→TERR_MASS→COOP" = c("COOP ~ TERR + MASS", "TERR ~ MASS"),
    "E0_TERR→COOP_MASS→TERR" = c("COOP ~ TERR", "TERR ~ MASS"),
    
    ## SERIES N: Null models
    "N0_MASS→FS_MASS→COOP_MASS→TERR" = c("FS ~ MASS", "COOP ~ MASS", "TERR ~ MASS"),
    
    ## SERIES X: COOP influences FS (direction 1)
    "X1_COOP→FS_MASS→FS" = c("FS ~ COOP + MASS"),
    "X1_COOP→FS_TERR→FS_MASS→FS" = c("FS ~ COOP + TERR + MASS"),
    "X1_COOP→FS_TERR→COOP_MASS→FS_MASS→COOP" = c("FS ~ COOP + MASS", "COOP ~ TERR + MASS"),
    "X1_COOP→FS_TERR→COOP_MASS→COOP" = c("FS ~ COOP", "COOP ~ TERR + MASS"),
    "X1_COOP→FS_TERR→COOP_MASS→FS" = c("FS ~ COOP + MASS", "COOP ~ TERR"),
    "X1_COOP→FS_TERR→FS_TERR→COOP_MASS→FS" = c("FS ~ TERR + COOP + MASS", "COOP ~ TERR"),
    "X1_COOP→FS_TERR→FS_TERR→COOP_MASS→COOP" = c("FS ~ TERR + COOP", "COOP ~ TERR + MASS"),
    "X1_COOP→FS_TERR→FS_TERR→COOP_MASS→FS_MASS→COOP" = c("FS ~ TERR + COOP + MASS", "COOP ~ TERR + MASS"),
    "X1_COOP→FS_TERR→COOP_MASS→FS_MASS→TERR" = c("FS ~ COOP + MASS", "COOP ~ TERR", "TERR ~ MASS"),
    "X1_COOP→FS_TERR→FS_TERR→COOP_MASS→FS_MASS→TERR" = c("FS ~ TERR + COOP + MASS", "COOP ~ TERR", "TERR ~ MASS"),
    "X1_COOP→FS_TERR→FS_TERR→COOP" = c("FS ~ TERR + COOP", "COOP ~ TERR"), # doesn't contain MASS, but "fully connected" if run in a model set without the 4th trait
    
    ## SERIES Y: FS influences COOP (direction 2)
    "Y2_FS→COOP_MASS→COOP" = c("COOP ~ FS + MASS"),
    "Y2_FS→COOP_TERR→COOP_MASS→COOP" = c("COOP ~ FS + TERR + MASS"),
    "Y2_TERR→FS_FS→COOP_MASS→FS_MASS→COOP" = c("FS ~ TERR + MASS", "COOP ~ FS + MASS"),
    "Y2_TERR→FS_FS→COOP_MASS→COOP" = c("FS ~ TERR", "COOP ~ FS + MASS"),
    "Y2_TERR→FS_FS→COOP_MASS→FS" = c("FS ~ TERR + MASS", "COOP ~ FS"),
    "Y2_TERR→FS_FS→COOP_TERR→COOP_MASS→FS" = c("FS ~ TERR + MASS", "COOP ~ FS + TERR"),
    "Y2_TERR→FS_FS→COOP_TERR→COOP_MASS→COOP" = c("FS ~ TERR", "COOP ~ FS + TERR + MASS"),
    "Y2_TERR→FS_FS→COOP_TERR→COOP_MASS→FS_MASS→COOP" = c("FS ~ TERR + MASS", "COOP ~ FS + TERR + MASS"),
    "Y2_TERR→FS_FS→COOP_MASS→FS_MASS→TERR" = c("FS ~ TERR + MASS", "COOP ~ FS", "TERR ~ MASS"),
    "Y2_TERR→FS_FS→COOP_TERR→COOP_MASS→FS_MASS→TERR" = c("FS ~ TERR + MASS", "COOP ~ FS + TERR", "TERR ~ MASS"),
    "Y2_TERR→FS_FS→COOP_TERR→COOP" = c("FS ~ TERR", "COOP ~ FS + TERR") # doesn't contain MASS, but "fully connected" if run in a model set without the 4th trait
  )
  
  model_patterns_with_TERR_response <- list(
    ## SERIES F: Basic models with TERR as response or an intermediate factor (no direct FS-COOP interaction)
    "F0_COOP→TERR" = c("TERR ~ COOP"),
    "F0_FS→TERR" = c("TERR ~ FS"),
    "F0_FS→TERR_COOP→TERR" = c("TERR ~ FS + COOP"),
    "F0_FS→TERR_TERR→COOP" = c("COOP ~ TERR", "TERR ~ FS"),
    "F0_COOP→TERR_TERR→FS" = c("FS ~ TERR", "TERR ~ COOP"),
    # Series F containing MASS:
    "F0_FS→TERR_TERR→COOP_MASS→FS" = c("COOP ~ TERR", "TERR ~ FS + MASS"),
    "F0_COOP→TERR_TERR→FS_MASS→COOP" = c("FS ~ TERR", "TERR ~ COOP + MASS"),
    "F0_FS→TERR_TERR→COOP_MASS→COOP" = c("COOP ~ TERR + MASS", "TERR ~ FS"),
    "F0_COOP→TERR_TERR→FS_MASS→FS" = c("FS ~ TERR + MASS", "TERR ~ COOP"),
    "F0_FS→TERR_TERR→COOP_MASS→FS_MASS→COOP" = c("COOP ~ TERR + MASS", "TERR ~ FS + MASS"),
    "F0_COOP→TERR_TERR→FS_MASS→FS_MASS→COOP" = c("FS ~ TERR + MASS", "TERR ~ COOP + MASS"),
    "F0_MASS→FS_FS→TERR_TERR→COOP" = c("COOP ~ TERR", "TERR ~ FS", "FS ~ MASS"),
    "F0_MASS→COOP_COOP→TERR_TERR→FS" = c("FS ~ TERR", "TERR ~ COOP", "COOP ~ MASS"),
    "F0_FS→TERR_MASS→TERR" = c("TERR ~ FS + MASS"),
    "F0_COOP→TERR_MASS→TERR" = c("TERR ~ COOP + MASS"),
    "F0_FS→TERR_COOP→TERR_MASS→TERR" = c("TERR ~ FS + COOP + MASS"),
    
    ## SERIES G: FS influences COOP and one/both influence TERR (FS → COOP)
    "G2_FS→COOP_COOP→TERR" = c("TERR ~ COOP", "COOP ~ FS"),
    "G2_FS→COOP_FS→TERR" = c("TERR ~ FS", "COOP ~ FS"),
    ## SERIES G containing MASS
    "G2_FS→COOP_COOP→TERR_MASS→FS" = c("TERR ~ COOP", "COOP ~ FS + MASS"),
    "G2_FS→COOP_COOP→TERR_MASS→TERR" = c("TERR ~ COOP + MASS", "COOP ~ FS"),
    "G2_FS→COOP_COOP→TERR_MASS→FS_MASS→TERR" = c("TERR ~ COOP + MASS", "COOP ~ FS + MASS"),
    "G2_MASS→FS_FS→COOP_COOP→TERR" = c("TERR ~ COOP", "COOP ~ FS", "FS ~ MASS"),
    "G2_FS→COOP_FS→TERR_COOP→TERR" = c("TERR ~ FS + COOP", "COOP ~ FS"), # doesn't contain MASS, but "fully connected" if run in a model set without the 4th trait
    
    ## SERIES H: COOP influences FS and one/both influence TERR (COOP → FS)
    "H1_COOP→FS_FS→TERR" = c("TERR ~ FS", "FS ~ COOP"),
    "H1_COOP→FS_COOP→TERR" = c("TERR ~ COOP", "FS ~ COOP"),
    ## SERIES H containing MASS
    "H1_COOP→FS_FS→TERR_MASS→COOP" = c("TERR ~ FS", "FS ~ COOP + MASS"),
    "H1_COOP→FS_FS→TERR_MASS→TERR" = c("TERR ~ FS + MASS", "FS ~ COOP"),
    "H1_COOP→FS_FS→TERR_MASS→COOP_MASS→TERR" = c("TERR ~ FS + MASS", "FS ~ COOP + MASS"),
    "H1_MASS→COOP_COOP→FS_FS→TERR" = c("TERR ~ FS", "FS ~ COOP", "COOP ~ MASS"),
    "H1_COOP→FS_FS→TERR_COOP→TERR" = c("TERR ~ FS + COOP", "FS ~ COOP")
  )
  
  #### select model sets ----
  
  if (mass_var != "none") {
    model_patterns = c(model_patterns, model_patterns_with_MASS)
    if (include_terr_response) {
      model_patterns = c(model_patterns, model_patterns_with_TERR_response)
    }
  }
  
  # Generate all the models
  models_list <- lapply(names(model_patterns), function(model_name) {
    create_formula_set(model_patterns[[model_name]], var_map)
  })
  names(models_list) <- names(model_patterns)
  
  # Create the final model set
  models <- do.call(define_model_set, models_list)
  
  
  #### run phylo_path
  result <- phylo_path(
    models, 
    data = dfIn_clean, 
    tree = tree_clean, 
    model = 'lambda',
    na.rm = TRUE
  )
  
  s <- summary(result)
  
  s_plot <- plot(s)
  
  # prep to plot maps
  phylopath_map_positions = data.frame("name" = c(female_song_var, coop_breeding_var, territoriality_var, mass_var), "x" = c(5, 8, 5, 2), "y" = c(5, 9, 1, 9) )
  
  best_model <- best(result)
  best_model_plot <- plot(best_model, text_size = 3, manual_layout = phylopath_map_positions)
  full_average <- average(result, cut_off = 2, avg_method = "full")
  full_average_plot <- plot(full_average, text_size = 3, manual_layout = phylopath_map_positions)
  conditional_average <- average(result, cut_off = 2, avg_method = "conditional")
  conditional_average_plot <- plot(conditional_average, text_size = 3, manual_layout = phylopath_map_positions) 
  
  CICsub2_models = s$model[which(s$delta_CICc<2)]
  print(paste("There are", length(CICsub2_models), "models with delta CICc < 2"))
  
  
  topModelPlotList <- list()
  topModelPlotTitles <- list()
  for (i in 1:length(CICsub2_models)) {
    tempmodel = s$model[which(s$delta_CICc<2)][i]
    tempCIC = s$CICc[which(s$delta_CICc<2)][i]
    temp_delta_CIC = s$delta_CICc[which(s$delta_CICc<2)][i]
    tempmodel_clean <- gsub("→", "to", tempmodel)
    tempmodel_title <- paste(tempmodel_clean, "\nCICc:", round(tempCIC, 2), "| delta_CICc:", round(temp_delta_CIC, 2))
    tempplot <- plot(choice(result, choice = tempmodel), text_size = 2, manual_layout = phylopath_map_positions)
    topModelPlotList[[tempmodel_clean]] <- tempplot
    topModelPlotTitles[[tempmodel_clean]] <- tempmodel_title
  }
  
  
  if (plots2pdf) {
    pdfname = paste("phylopath plots", female_song_var, coop_breeding_var, territoriality_var, mass_var, Sys.Date(), ".pdf")
    
    # Required library
    require(gridExtra)
    require(ggpubr)
    
    # Create a PDF device
    pdf(pdfname, width = 15, height = 8)
    
    # Plot the summary plot
    print(s_plot)
    
    # Plot top models with titles, up to 6 per page
    num_plots <- length(topModelPlotList)
    plot_counter <- 0
    plots_per_page <- 4
    page_plots <- list()
    
    for (i in 1:num_plots) {
      model_name <- names(topModelPlotList)[i]
      current_plot <- topModelPlotList[[model_name]]
      title <- topModelPlotTitles[[model_name]]
      
      # Add title to the plot
      titled_plot <- current_plot + 
        ggtitle(title) + 
        theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
              plot.margin = margin(0, 10, 0, 10)) +  # Larger margins
        coord_cartesian(xlim = c(-1, 11), ylim = c(-1, 11), expand = TRUE)  # Expand plot area 
      
      plot_counter <- plot_counter + 1
      page_plots[[plot_counter]] <- titled_plot
      
      # When we reach plots_per_page or the last plot, print the page
      if (plot_counter == plots_per_page || i == num_plots) {
        # Arrange and print the plots
        grid_layout <- grid.arrange(grobs = page_plots, ncol = 2)
        print(grid_layout)
        
        # Reset for next page
        plot_counter <- 0
        page_plots <- list()
      }
    }
    
    # Print full average model with title
    print(full_average_plot + 
            ggtitle("Full averaged model") + 
            theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
                  plot.margin = margin(0, 10, 0, 10)) +  # Larger margins
            coord_cartesian(xlim = c(-1, 11), ylim = c(-1, 11), expand = TRUE))  # Expand plot area 
    
    # Print conditional average model with title
    print(conditional_average_plot + 
            ggtitle("Conditional averaged model") + 
            theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
                  plot.margin = margin(0, 10, 0, 10)) +  # Larger margins
            coord_cartesian(xlim = c(-1, 11), ylim = c(-1, 11), expand = TRUE))  # Expand plot area 
    
    # Close the PDF device
    dev.off()
  }
  
  # Create CICc bar plot
  cicbar_plot <- plot_cicbar(s)
  
  # Add better versions of the DAG plots using the new function
  best_model_plot_enhanced <- plot_phylopath_dag(best_model, var_map, "Best Model (Lowest CICc)")
  full_average_plot_enhanced <- plot_phylopath_dag(full_average, var_map, "Fully-Averaged Models (Δ CICc < 2)")
  conditional_average_plot_enhanced <- plot_phylopath_dag(conditional_average, var_map, "Conditional-Averaged Models (Δ CICc < 2)")
  
  # Save to PNG if requested
  if (plots2png) {
    # Create file base name
    var_string <- paste(c(female_song_var, coop_breeding_var, territoriality_var, mass_var), collapse = " ")
    file_base <- paste0("phylopath ", var_string, " ", nSpecies, "species")
    
    enhanced_plots <- list(
      best_model_plot = best_model_plot_enhanced,
      full_average_plot = full_average_plot_enhanced,
      conditional_average_plot = conditional_average_plot_enhanced,
      cicbar_plot = cicbar_plot,
      summary_plot = s_plot
    )
    
    save_phylopath_plots_png(enhanced_plots, file_base, output_dir)
  }
  
  
  outlist <- list(
    var_map = var_map,
    result = result,
    summary_plot = s_plot,
    best_model_plot = best_model_plot,
    full_average_plot = full_average_plot,
    conditional_average_plot = conditional_average_plot,
    best_model_plot_enhanced = best_model_plot_enhanced,  # NEW
    full_average_plot_enhanced = full_average_plot_enhanced,  # NEW
    conditional_average_plot_enhanced = conditional_average_plot_enhanced,  # NEW
    cicbar_plot = cicbar_plot,  # NEW
    topModelPlotList = topModelPlotList,
    topModelPlotTitles = topModelPlotTitles,
    CICsub2_models = CICsub2_models,
    nSpecies = nSpecies
  )
  
  return(outlist)
  
} # end function run_CB_FS_Terr_phylopath



## Function to downsample dataset and run run_CB_FS_Terr_phylopath() - updated to match downsample_run_phylopath_fixed()
downsample_run_phylopath <- function(dfIn, tree, downsample_columns, downsample_values, numToRemove, seed, female_song_var = "FemaleSong_Agg01", coop_breeding_var = "HighConfidence_Coop", territoriality_var = "TerritorialityWeakVsStrong", mass_var = "none") {
  
  # Add geographic regions if not already present
  if (!("GeographicRegion_Jetz" %in% colnames(dfIn))) {
    dfIn$GeographicRegion_Jetz = NA
    dfIn$GeographicRegion_Jetz[which(dfIn$Realm_Jetz2011 %in% c("PA", "NeA"))] <- "Holarctic"
    dfIn$GeographicRegion_Jetz[which(dfIn$Realm_Jetz2011 %in% c("AT", "NT", "IM", "AA", "OC"))] <- "Tropical"
  }
  
  # Subset to just Oscine species
  dfOs = dfIn[which(dfIn$species %in% tree$tip.label),]
  
  # Add data columns if not present
  if (!("HaveFSData" %in% colnames(dfOs))) {
    dfOs$HaveFSData = !is.na(dfOs[[female_song_var]])
  }
  if (!("HaveCBData" %in% colnames(dfOs))) {
    dfOs$HaveCBData = !is.na(dfOs[[coop_breeding_var]])
  }
  if (!("HaveFSCBData" %in% colnames(dfOs))) {
    dfOs$HaveData = !is.na(dfOs[[female_song_var]]) & !is.na(dfOs[[coop_breeding_var]])
  }
  
  # First, identify the rows that meet all criteria
  if (length(downsample_columns) == 2) {
    matching_rows <- which(dfOs[[downsample_columns[1]]] == downsample_values[1] & 
                             dfOs[[downsample_columns[2]]] == downsample_values[2])
  } else if (length(downsample_columns) == 3) {
    matching_rows <- which(dfOs[[downsample_columns[1]]] == downsample_values[1] & 
                             dfOs[[downsample_columns[2]]] == downsample_values[2] & 
                             dfOs[[downsample_columns[3]]] == downsample_values[3])
  } else {
    stop("downsample_columns must have length 2 or 3")
  }
  
  # Check how many matching rows we have
  num_matching <- length(matching_rows)
  print(paste("Number of rows matching criteria:", num_matching))
  
  # Make sure there are at least numToRemove rows to remove
  if (num_matching < numToRemove) {
    stop(paste("There are only", num_matching, "rows matching the criteria, but trying to remove", numToRemove))
  } else {
    # Randomly select numToRemove rows to drop
    set.seed(seed)  # For reproducibility
    rows_to_drop <- sample(matching_rows, numToRemove)
    
    # Create the downsampled dataframe by removing the selected rows
    dfOs_downsampled <- dfOs[-rows_to_drop, ]
  }
  
  output <- run_CB_FS_Terr_phylopath(dfIn = dfOs_downsampled, tree = tree, female_song_var = female_song_var, coop_breeding_var = coop_breeding_var, territoriality_var = territoriality_var, mass_var = mass_var, plots2pdf = FALSE)
  return(output)
  
}

#' Run phylopath with dimorphism bias correction and create detailed models dataframe
#'
#' This function replaces the code in Run_Analyses.R section "run phylopath with dimorphism bias correction/downsampling"
#' It properly creates a detailed_models_df object that can be used with create_downsampled_plots()
#'
#' @param dfIn_phylo Full dataset with species data
#' @param tree Phylogenetic tree
#' @param dim_info List with dimorphism information (col, label, description)
#' @param n_iterations Number of iterations to run
#' @param phylopath_output_dir Output directory for plots
#' @param save_outputs Whether to save CSV files and plots
#' @return List containing results and detailed_models_df
run_phylopath_dimorphism_correction <- function(dfIn_phylo, 
                                                tree,
                                                dim_info,
                                                n_iterations = 500,
                                                save_outputs = TRUE,
                                                female_song_var = "FemaleSong_Agg01",
                                                coop_breeding_var = "HighConfidence_Coop",
                                                territoriality_var = "TerritorialityWeakVsStrong",
                                                mass_var = "logMass_AVONET",
                                                phylopath_output_dir = "obsolete"
                                                ) {
  
  require(gridExtra)
  
  
  all_traits_phylopath_label <- paste(female_song_var, coop_breeding_var, territoriality_var, mass_var)
  
  phylopath_output_dir = file.path("Outputs","PhylopathDownsampled", paste(all_traits_phylopath_label, "models"))
  
  # Create output directory specific to this dimorphism type
  dim_output_dir <- file.path(phylopath_output_dir, dim_info$label)
  if (!dir.exists(dim_output_dir)) {
    dir.create(dim_output_dir, recursive = TRUE)
  }
  
  cat("\n\n========================================\n")
  cat("Running", dim_info$description, "bias correction...\n")
  cat("========================================\n")
  
  # Check if the column exists
  if (!dim_info$col %in% colnames(dfIn_phylo)) {
    cat("Warning: Column", dim_info$col, "not found. Skipping...\n")
    return(NULL)
  }
  
  # Run the dimorphism downsampling analysis
  dimorphism_downsample <- downsample_dimorphism_bias(
    df = dfIn_phylo,
    dimorphism_col = dim_info$col,
    data_col = "FemaleSong_Agg01",
    n_iterations = n_iterations
  )
  
  # Create distribution plot if requested
  if (save_outputs) {
    prefix_dimorphism <- paste0("Remove", 
                                dimorphism_downsample$n_to_remove, 
                                "High", 
                                dim_info$label)
    
    dist_plot_result <- create_downsampled_dimorphism_distribution_plot( 
      dfIn_phylo = dfIn_phylo,
      dimorphism_downsample = dimorphism_downsample,
      dim_info = dim_info,
      output_dir = dim_output_dir,
      prefix = prefix_dimorphism,
      save_plot = TRUE
    )
    
    # Create the violin plot (multiple iterations view)
    violin_plot_result <- create_dimorphism_downsamples_violin_plot(
      dfIn_phylo = dfIn_phylo,
      dimorphism_downsample = dimorphism_downsample,
      dim_info = dim_info,
      output_dir = dim_output_dir,
      prefix = prefix_dimorphism,
      save_plot = TRUE,
      show_iterations = min(20, n_iterations)
    )
  }
  
  cat("Need to remove", dimorphism_downsample$n_to_remove, "species to correct", dim_info$description, "bias\n")
  cat("Current mean:", round(dimorphism_downsample$current_stats$current_mean, 3), "\n")
  cat("Target mean:", round(dimorphism_downsample$target_stats$target_mean, 3), "\n")
  
  # Run phylopath on each iteration
  all_results <- list()
  all_model_summaries <- list()
  detailed_models_list <- list()
  
  for (i in 1:n_iterations) {
    cat(sprintf("Running iteration %d/%d...\r", i, n_iterations))
    
    # Get the downsampled dataset for this iteration
    df_downsampled <- dimorphism_downsample$iterations[[i]]$remaining_data
    
    # Filter to species in tree
    df_downsampled <- df_downsampled[df_downsampled$species %in% tree$tip.label, ]
    rownames(df_downsampled) <- df_downsampled$species
    
    # Check if we have enough species
    if (nrow(df_downsampled) < 50) {
      warning(paste("Iteration", i, "has only", nrow(df_downsampled), "species. Skipping..."))
      next
    }
    
    # Run phylopath on this iteration
    result_i <- run_CB_FS_Terr_phylopath(
      dfIn = df_downsampled,
      tree = tree,
      female_song_var = female_song_var,
      coop_breeding_var = coop_breeding_var,
      territoriality_var = territoriality_var,
      mass_var = mass_var,  
      plots2pdf = FALSE
    )
    
    all_results[[i]] <- result_i
    
    # Extract model summaries and create detailed models entries
    if (!is.null(result_i$result)) {
      summary_i <- summary(result_i$result)
      
      if (!is.null(summary_i)) {
        # For each model with delta_CICc < 2
        for (model_idx in which(summary_i$delta_CICc < 2)) {
          model_name <- summary_i$model[model_idx]
          model_CICc <- summary_i$CICc[model_idx]
          model_delta_CICc <- summary_i$delta_CICc[model_idx]
          
          # Get the chosen model to extract path coefficients
          chosen_model <- choice(result_i$result, model_name)
          
          # Create a row for this model
          model_row <- data.frame(
            seed = i,  # Using iteration number as seed
            model = model_name,
            CICc = model_CICc,
            delta_CICc = model_delta_CICc,
            nSpecies = result_i$nSpecies,
            stringsAsFactors = FALSE
          )
          
          # Extract edge information from the coefficient matrix
          if (!is.null(chosen_model$coef) && is.matrix(chosen_model$coef)) {
            edges_coef <- chosen_model$coef
            edges_se <- chosen_model$se
            
            # Process coefficient matrix
            for (from_idx in 1:nrow(edges_coef)) {
              for (to_idx in 1:ncol(edges_coef)) {
                coef_value <- edges_coef[from_idx, to_idx]
                
                # Only process non-zero coefficients (actual paths)
                if (coef_value != 0) {
                  from_name <- rownames(edges_coef)[from_idx]
                  to_name <- colnames(edges_coef)[to_idx]
                  
                  # Get standard error
                  se_value <- edges_se[from_idx, to_idx]
                  
                  # Calculate approximate p-value
                  z_value <- coef_value / se_value
                  p_value <- 2 * (1 - pnorm(abs(z_value)))
                  
                  # Create column names for this edge
                  edge_name <- paste0(from_name, "_to_", to_name)
                  edge_est_name <- paste0(edge_name, "_est")
                  edge_se_name <- paste0(edge_name, "_se")
                  edge_p_name <- paste0(edge_name, "_p")
                  
                  # Add to model row
                  model_row[[edge_est_name]] <- coef_value
                  model_row[[edge_se_name]] <- se_value
                  model_row[[edge_p_name]] <- p_value
                }
              }
            }
          }
          
          # Add to detailed models list
          detailed_models_list[[length(detailed_models_list) + 1]] <- model_row
        }
        
        # Also create summary for aggregated results
        model_summary_i <- as.data.frame(summary_i)
        model_summary_i$iteration <- i
        model_summary_i$model_name <- rownames(model_summary_i)
        all_model_summaries[[i]] <- model_summary_i
      }
    }
  }
  
  cat("\n")  # New line after progress indicator
  
  # Combine detailed models into a single data frame
  detailed_models_df <- bind_rows(detailed_models_list)
  
  # Aggregate results across iterations
  aggregated_summaries <- if (length(all_model_summaries) > 0) {
    do.call(rbind, all_model_summaries)
  } else {
    data.frame()
  }
  
  # Calculate average model performance
  avg_model_performance <- if (nrow(aggregated_summaries) > 0) {
    aggregated_summaries %>%
      group_by(model_name) %>%
      summarise(
        mean_CICc = mean(CICc, na.rm = TRUE),
        sd_CICc = sd(CICc, na.rm = TRUE),
        times_best = sum(delta_CICc == 0, na.rm = TRUE),
        times_in_top = sum(delta_CICc < 2, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(mean_CICc)
  } else {
    data.frame()
  }
  
  # Calculate model frequencies for compatibility with create_downsampled_plots
  model_frequencies <- if (nrow(detailed_models_df) > 0) {
    detailed_models_df %>%
      group_by(model) %>%
      summarise(
        frequency_in_sub2 = n(),
        .groups = "drop"
      ) %>%
      rename(model_name = model) %>%
      arrange(desc(frequency_in_sub2))
  } else {
    data.frame()
  }
  
  # Create a summary phylopath output structure compatible with run_multiple_phylopath output
  result_dimorphism <- list(
    results_df = data.frame(
      seed = 1:n_iterations,
      nSub2dCIC_Models = sapply(all_model_summaries, function(x) sum(x$delta_CICc < 2, na.rm = TRUE)),
      stringsAsFactors = FALSE
    ),
    model_frequencies = model_frequencies,
    detailed_models = detailed_models_df,
    model_summaries = aggregated_summaries,
    avg_performance = avg_model_performance,
    n_iterations = n_iterations,
    downsampling_info = dimorphism_downsample,
    all_results = all_results,
    dimorphism_type = dim_info$label
  )
  
  # Save outputs if requested
  if (save_outputs) {
    # Save detailed models CSV
    detailed_models_file <- file.path(dim_output_dir, 
                                      paste0("detailed_models_", 
                                             "Remove", dimorphism_downsample$n_to_remove,
                                             "High", dim_info$label,
                                             "_n", n_iterations, "_", Sys.Date(), ".csv"))
    write.csv(detailed_models_df, detailed_models_file, row.names = FALSE)
    cat("Detailed models saved to:", detailed_models_file, "\n")
    
    # Save model frequencies
    write.csv(model_frequencies, 
              file.path(dim_output_dir, 
                        paste0("model_frequencies_", dim_info$label, "_n", n_iterations, "_", Sys.Date(), ".csv")),
              row.names = FALSE)
  }
  
  cat("\n", dim_info$description, "bias correction complete!\n")
  cat("Results saved to:", dim_output_dir, "\n")
  
  return(result_dimorphism)
}


#### Helper functions to build formulas from patterns ----
build_formula <- function(pattern, var_map) {
  for (var_name in names(var_map)) {
    pattern <- gsub(var_name, var_map[[var_name]], pattern, fixed = TRUE)
  }
  return(as.formula(pattern))
}

# Helper function to create a set of formulas
create_formula_set <- function(patterns, var_map) {
  formulas <- lapply(patterns, function(pattern) {
    build_formula(pattern, var_map)
  })
  return(do.call(c, formulas))
}

# Function to calculate Akaike weights from delta CICc values
calculate_akaike_weights <- function(delta_CICc, cutoff = 2) {
  # Filter models within cutoff
  valid_models <- delta_CICc <= cutoff
  filtered_deltas <- delta_CICc[valid_models]
  
  # Calculate raw weights: exp(-0.5 * delta_CICc)
  raw_weights <- exp(-0.5 * filtered_deltas)
  
  # Normalize to sum to 1 (relative weights)
  rel_weights <- raw_weights / sum(raw_weights)
  
  return(list(
    valid_indices = which(valid_models),
    rel_weights = rel_weights
  ))
}

# Helper function for null coalescing
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) y else x
}

# Downsample based on probabilities weighted by the magnitude of dimorphism so more-dimorphic species are more likely to be removed
downsample_dimorphism_bias <- function(df, 
                                       dimorphism_col = "logMaleFemalePlumageDiffAbs",
                                       data_col = "FemaleSong_Agg01",
                                       n_iterations = 500) {
  
  # Calculate target statistics from no-data group
  no_data_stats <- df %>%
    filter(is.na(!!sym(data_col))) %>%
    summarise(
      target_mean = mean(!!sym(dimorphism_col), na.rm = TRUE),
      target_sd = sd(!!sym(dimorphism_col), na.rm = TRUE),
      target_median = median(!!sym(dimorphism_col), na.rm = TRUE)
    )
  
  # Calculate current statistics for has-data group
  has_data_stats <- df %>%
    filter(!is.na(!!sym(data_col))) %>%
    summarise(
      current_mean = mean(!!sym(dimorphism_col), na.rm = TRUE),
      current_sd = sd(!!sym(dimorphism_col), na.rm = TRUE),
      n_species = n()
    )
  
  # Calculate how many species to remove
  # This could be refined based on distribution matching
  current_mean = has_data_stats$current_mean
  target_mean = no_data_stats$target_mean
  excess_dimorphism <- has_data_stats$current_mean - no_data_stats$target_mean
  n_to_remove <- round(has_data_stats$n_species * (excess_dimorphism / has_data_stats$current_mean))
  
  # Calculate removal probabilities
  df_with_data <- df %>%
    filter(!is.na(!!sym(data_col))) %>%
    mutate(
      dimorphism_excess = !!sym(dimorphism_col) - no_data_stats$target_mean,
      # Use logistic function to convert to probabilities
      removal_weight = case_when(
        dimorphism_excess > 0 ~ 1 + dimorphism_excess,  # Higher weight for more dimorphic
        TRUE ~ exp(dimorphism_excess)  # Lower weight for less dimorphic
      )
    )
  
  # Those species that did not have dimorphism data get the median removal weight
  num_NA_dimorphism = sum(is.na(df_with_data$removal_weight))
  median_removal_weight = median(df_with_data$removal_weight, na.rm = T)
  df_with_data$removal_weight[which(is.na(df_with_data$removal_weight))] <- median_removal_weight
  print(paste("There are", num_NA_dimorphism, "species with", data_col, "data but without", dimorphism_col, "data. These species were assigned the median removal weight (before normalization) of", median_removal_weight))
  
  # Normalize weights
  df_with_data$removal_prob <- df_with_data$removal_weight / sum(df_with_data$removal_weight)
  
  plot(df_with_data[,dimorphism_col],df_with_data$removal_prob)
  plot(density(df_with_data$removal_prob))
  
  # Plot distributions
  n_species <- length(df_with_data[,1])
  
  # Create top plot (scatter plot)
  p1 <- ggplot(df_with_data, aes_string(x = dimorphism_col, y = "removal_prob")) +
    geom_point() +
    labs(
      title = paste("Removal Probabilities vs", dimorphism_col),
      subtitle = paste("N =", n_species, "species | Current mean:", round(current_mean, 3), 
                       "| Target mean:", round(target_mean, 3), "| To remove:", n_to_remove),
      x = dimorphism_col,
      y = "Removal Probability"
    ) +
    theme_minimal()
  
  # Create middle plot (density plot of the dimorphism_col)
  p_dimorph <- ggplot(df_with_data, aes_string(x = dimorphism_col)) +
    geom_density(fill = "lightblue", alpha = 0.7) +
    labs(
      title = "Dimorphism distribution",
      subtitle = paste("Distribution for", n_species, "species"),
      x = dimorphism_col,
      y = "Density"
    ) +
    theme_minimal()
  
  # Create bottom plot (density plot)
  p2 <- ggplot(df_with_data, aes(x = removal_prob)) +
    geom_density(fill = "lightblue", alpha = 0.7) +
    labs(
      title = "Density of Removal Probabilities",
      subtitle = paste("Distribution for", n_species, "species"),
      x = "Removal Probability",
      y = "Density"
    ) +
    theme_minimal()
  
  # Combine plots
  combined_plot <- grid.arrange(p1, p_dimorph, p2, ncol = 1)
  
  # Create filename
  weighted_probs_filename <- paste0("phylopath ", data_col, " ", dimorphism_col, 
                     " downsampling weighted removal probabilities for ", 
                     n_species, " species.png")
  
  # Save to Outputs folder
  ggsave(
    filename = file.path("Outputs", "PhylopathDownsampled", weighted_probs_filename),
    plot = combined_plot,
    width = 10,
    height = 12,
    dpi = 300
  )
  
  
  # Perform iterations
  results <- list()
  for (i in 1:n_iterations) {
    sampled_species <- sample(df_with_data$species,
                              size = n_to_remove,
                              prob = df_with_data$removal_prob,
                              replace = FALSE)
    
    results[[i]] <- list(
      removed_species = sampled_species,
      remaining_data = df[!df$species %in% sampled_species, ]
    )
  }
  
  return(list(
    n_to_remove = n_to_remove,
    target_stats = no_data_stats,
    current_stats = has_data_stats,
    iterations = results
  ))
}

#' Calculate territoriality bias downsampling # this is done elsewhere too... maybe calculate_stratified_downsampling_with_territoriality.R?
#'
#' @param df Data frame with species data
#' @param territoriality_col Column name for territoriality data
#' @param territory_value_high Value representing high territoriality
#' @param territory_value_low Value representing low territoriality
#' @param data_cols Columns that must have data (e.g., c("FemaleSong_Agg01", "HighConfidence_Coop"))
#' @return List with downsampling information
calculate_territoriality_downsampling <- function(df,
                                                  territoriality_col = "TerritorialityWeakVsStrong",
                                                  territory_value_high = "1",
                                                  territory_value_low = "0",
                                                  data_cols = c("FemaleSong_Agg01", "HighConfidence_Coop")) {
  
  # Check if required columns exist
  missing_cols <- setdiff(c(territoriality_col, data_cols), colnames(df))
  if (length(missing_cols) > 0) {
    stop("Missing columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Add HaveData indicator (all specified data columns must be non-NA)
  df$HaveData <- apply(df[data_cols], 1, function(x) all(!is.na(x)))
  
  # Convert territoriality to character for comparison
  df[[territoriality_col]] <- as.character(df[[territoriality_col]])
  
  # Calculate proportions for high and low territoriality
  high_terr <- df %>%
    filter(!!sym(territoriality_col) == territory_value_high) %>%
    summarise(
      species_with_data = sum(HaveData),
      total_species = n(),
      proportion = ifelse(n() > 0, species_with_data / total_species, 0)
    )
  
  low_terr <- df %>%
    filter(!!sym(territoriality_col) == territory_value_low) %>%
    summarise(
      species_with_data = sum(HaveData),
      total_species = n(),
      proportion = ifelse(n() > 0, species_with_data / total_species, 0)
    )
  
  # Check if we have data for both groups
  if (nrow(high_terr) == 0 || high_terr$total_species == 0) {
    warning(paste("No species found with", territoriality_col, "=", territory_value_high))
    high_terr <- data.frame(species_with_data = 0, total_species = 0, proportion = 0)
  }
  
  if (nrow(low_terr) == 0 || low_terr$total_species == 0) {
    warning(paste("No species found with", territoriality_col, "=", territory_value_low))
    low_terr <- data.frame(species_with_data = 0, total_species = 0, proportion = 0)
  }
  
  # Create proportions data frame
  proportions <- data.frame(
    territoriality = c("High", "Low"),
    species_with_data = c(high_terr$species_with_data, low_terr$species_with_data),
    total_species = c(high_terr$total_species, low_terr$total_species),
    proportion = c(high_terr$proportion, low_terr$proportion)
  )
  
  # Check if we can perform downsampling
  if (high_terr$total_species == 0 && low_terr$total_species == 0) {
    return(list(
      proportions = proportions,
      target_proportion = 0,
      n_to_remove = 0,
      downsample_info = list(
        group_to_downsample = "None - no data available",
        territory_value = NA
      ),
      calculation = "No species found with specified territoriality values"
    ))
  }
  
  # Determine which group has higher proportion and needs downsampling
  if (high_terr$proportion > low_terr$proportion) {
    # Need to downsample high territoriality group
    target_proportion <- low_terr$proportion
    n_to_remove <- high_terr$species_with_data - round(high_terr$total_species * target_proportion)
    group_to_downsample <- "High territoriality"
    territory_value <- territory_value_high
  } else {
    # Need to downsample low territoriality group
    target_proportion <- high_terr$proportion
    n_to_remove <- low_terr$species_with_data - round(low_terr$total_species * target_proportion)
    group_to_downsample <- "Low territoriality"
    territory_value <- territory_value_low
  }
  
  # Ensure we don't try to remove more species than available
  n_to_remove <- max(0, n_to_remove)
  
  return(list(
    proportions = proportions,
    target_proportion = target_proportion,
    n_to_remove = n_to_remove,
    downsample_info = list(
      group_to_downsample = group_to_downsample,
      territory_value = territory_value
    ),
    calculation = paste("Remove", n_to_remove, "species from", group_to_downsample, 
                        "to match proportion of", 
                        ifelse(group_to_downsample == "High territoriality", "Low", "High"),
                        "territoriality group")
  ))
}


#### Single phylopath plotting functions ----
# Function to create CICc bar plot - new 2025-06-09
plot_cicbar <- function(summary_result, n_models = 20) {
  model_data <- summary_result %>%
    mutate(is_top = delta_CICc < 2) %>%
    arrange(CICc)
  
  n_models <- min(n_models, nrow(model_data))
  model_data_plot <- model_data[1:n_models, ]
  
  ggplot(model_data_plot, 
         aes(x = reorder(model, -CICc),
             y = CICc,
             fill = is_top)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_text(aes(label = round(delta_CICc, 2)), 
              hjust = -0.2, size = 3) +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "#2E86AB", "FALSE" = "#A7A7A7"),
                      labels = c("TRUE" = "Δ CICc < 2", "FALSE" = "Δ CICc ≥ 2"),
                      name = "Model Set") +
    labs(
      title = "Model Comparison by CICc",
      subtitle = paste("Showing top", n_models, "models. Values show Δ CICc"),
      x = NULL,
      y = "CICc"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(size = 9),
      legend.position = "bottom"
    )
}

# Function to create a single DAG plot with better layout - new 2025-06-09
plot_phylopath_dag <- function(model, var_map, title = "", 
                               text_size = 3, box_x = 22, box_y = 16) {
  # Set up layout positions based on var_map
  n_vars <- length(var_map)
  
  if (n_vars == 4) {
    phylopath_map_positions <- data.frame(
      name = unlist(var_map),
      x = c(5, 8, 5, 2),
      y = c(5, 9, 1, 9),
      stringsAsFactors = FALSE
    )
  } else {
    angles <- seq(0, 2*pi, length.out = n_vars + 1)[1:n_vars]
    phylopath_map_positions <- data.frame(
      name = unlist(var_map),
      x = 5 + 3 * cos(angles),
      y = 5 + 3 * sin(angles),
      stringsAsFactors = FALSE
    )
  }
  
  plot(model, 
       manual_layout = phylopath_map_positions,
       text_size = text_size,
       box_x = box_x, box_y = box_y) +
    ggtitle(title) +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          plot.margin = margin(0, 0, 0, 20)) +
    coord_cartesian(xlim = c(-1, 11), ylim = c(-1, 11), expand = TRUE)
}

# Function to save plots as PNG - new 2025-06-09
save_phylopath_plots_png <- function(plots, file_base, output_dir = "Outputs/PhylopathPlots") {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save individual plots
  if (!is.null(plots$best_model_plot)) {
    ggsave(file.path(output_dir, paste0(file_base, " best_model.png")), 
           plots$best_model_plot, width = 15, height = 10, dpi = 300)
  }
  
  if (!is.null(plots$full_average_plot)) {
    ggsave(file.path(output_dir, paste0(file_base, " full_avg.png")), 
           plots$full_average_plot, width = 15, height = 10, dpi = 300)
  }
  
  if (!is.null(plots$conditional_average_plot)) {
    ggsave(file.path(output_dir, paste0(file_base, " cond_avg.png")), 
           plots$conditional_average_plot, width = 15, height = 10, dpi = 300)
  }
  
  if (!is.null(plots$cicbar_plot)) {
    ggsave(file.path(output_dir, paste0(file_base, " cicbar.png")), 
           plots$cicbar_plot, width = 10, height = 8, dpi = 300)
  }
  
  if (!is.null(plots$summary_plot)) {
    ggsave(file.path(output_dir, paste0(file_base, " summary.png")), 
           plots$summary_plot, width = 12, height = 10, dpi = 300)
  }
  
  cat("PNG files saved to:", output_dir, "\n")
}

#### Downsampled phylopath analyses and plots ----

library(dplyr)
library(ggplot2)
library(tidyr)
library(forcats)

# Function 1: Create path summary statistics
create_path_summary <- function(detailed_models_df) {
  # Get all columns that contain coefficients (ending with "_est")
  coef_cols <- grep("_est$", names(detailed_models_df), value = TRUE)
  
  # Initialize an empty list to store path information
  path_data_list <- list()
  
  # For each coefficient column, extract the information
  for (col in coef_cols) {
    # Get the path name from the column name
    path_name <- gsub("_est$", "", col)
    # Get the p-value column
    p_col <- gsub("_est$", "_p", col)
    # Get the standard error column
    se_col <- gsub("_est$", "_se", col)
    
    # Only include rows where this path exists (non-NA coefficient)
    valid_rows <- !is.na(detailed_models_df[[col]])
    if (any(valid_rows)) {
      # Extract the data for this path
      path_df <- data.frame(
        seed = detailed_models_df$seed[valid_rows],
        model = detailed_models_df$model[valid_rows],
        nSpecies = detailed_models_df$nSpecies[valid_rows], 
        path = path_name,
        from = gsub("_to_.*$", "", path_name),
        to = gsub("^.*_to_", "", path_name),
        coefficient = detailed_models_df[[col]][valid_rows],
        std_error = detailed_models_df[[se_col]][valid_rows],
        p_value = detailed_models_df[[p_col]][valid_rows],
        significant = detailed_models_df[[p_col]][valid_rows] < 0.05,
        stringsAsFactors = FALSE
      )
      
      # Add to the list
      path_data_list[[length(path_data_list) + 1]] <- path_df
    }
  }
  
  # Combine all path data
  if (length(path_data_list) > 0) {
    path_coefficients_df <- bind_rows(path_data_list)
    
    # Create path summary
    path_summary <- path_coefficients_df %>%
      group_by(path, from, to) %>%
      summarize(
        n_occurrences = n(),
        mean_coefficient = mean(coefficient, na.rm = TRUE),
        median_coefficient = median(coefficient, na.rm = TRUE),
        sd_coefficient = sd(coefficient, na.rm = TRUE),
        min_coefficient = min(coefficient, na.rm = TRUE),
        max_coefficient = max(coefficient, na.rm = TRUE),
        q25_coefficient = quantile(coefficient, 0.25, na.rm = TRUE),
        q75_coefficient = quantile(coefficient, 0.75, na.rm = TRUE),
        n_significant = sum(significant, na.rm = TRUE),
        prop_significant = n_significant / n_occurrences,
        mean_p = mean(p_value, na.rm = TRUE),
        median_p = median(p_value, na.rm = TRUE),
        mean_nSpecies = mean(nSpecies, na.rm = TRUE),
        min_nSpecies = min(nSpecies, na.rm = TRUE),
        max_nSpecies = max(nSpecies, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(desc(prop_significant), desc(abs(mean_coefficient)))
    
    return(list(
      path_coefficients = path_coefficients_df,
      path_summary = path_summary
    ))
  } else {
    return(NULL)
  }
}

# Function 2: Create model summary statistics
create_model_summary <- function(detailed_models_df) {
  model_summary <- detailed_models_df %>%
    group_by(model) %>%
    summarize(
      n_iterations = n(),
      mean_CICc = mean(CICc, na.rm = TRUE),
      sd_CICc = sd(CICc, na.rm = TRUE),
      mean_delta_CICc = mean(delta_CICc, na.rm = TRUE),
      sd_delta_CICc = sd(delta_CICc, na.rm = TRUE),
      prop_best_model = sum(delta_CICc == 0) / n(),
      .groups = "drop"
    ) %>%
    arrange(desc(n_iterations), mean_delta_CICc)
  
  return(model_summary)
}

# Function 3: Create coefficient violin/box plots
plot_coefficient_distributions <- function(path_data, show_seeds = FALSE, 
                                           plot_type = "violin", jitter_alpha = 0.6) {

  
  path_data$path_clean <- path_data$path
  
  if ("mean_coefficient" %in% colnames(path_data) & !"coefficient" %in% colnames(path_data)) {
    # Base plot
    p <- ggplot(path_data, aes(x = path_clean, y = mean_coefficient))
    subtitleText = paste("Each point represents the mean coefficient across top models within each seed")
  } else {
    # Base plot
    p <- ggplot(path_data, aes(x = path_clean, y = coefficient))
    subtitleText = paste("Each point represents a path coefficient from one model iteration")
  }
  
  # Add violin or box plot
  if (plot_type == "violin") {
    p <- p + geom_violin(alpha = 0.7, fill = "lightblue", color = "darkblue")
  } else if (plot_type == "box") {
    p <- p + geom_boxplot(alpha = 0.7, fill = "lightblue", color = "darkblue", 
                          outlier.shape = NA)  # Don't show outliers since we'll add points
  }
  
  # Add horizontal line at zero
  p <- p + geom_hline(yintercept = 0, linetype = "dashed", color = "gray50")
  
  # Formatting
  p <- p + 
    labs(
      title = "Distribution of Path Coefficients",
      subtitle = subtitleText,
      x = "Path",
      y = "Coefficient Estimate"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "gray60")
    )
  
  return(p)
}

# Function 4: Create p-value violin/box plots
plot_pvalue_distributions <- function(path_data, show_seeds = FALSE, 
                                      plot_type = "violin", jitter_alpha = 0.6) {

  path_data$path_clean = path_data$path
  
  # Base plot
  p <- ggplot(path_data, aes(x = path_clean, y = p_value))
  
  # Add violin or box plot
  if (plot_type == "violin") {
    p <- p + geom_violin(alpha = 0.7, fill = "lightcoral", color = "darkred")
  } else if (plot_type == "box") {
    p <- p + geom_boxplot(alpha = 0.7, fill = "lightcoral", color = "darkred", 
                          outlier.shape = NA)  # Don't show outliers since we'll add points
  }
  
  # Add horizontal line at 0.05
  p <- p + geom_hline(yintercept = 0.05, linetype = "dashed", color = "gray50")
  
  # Formatting
  p <- p + 
    labs(
      title = "Distribution of Path P-values",
      subtitle = "Dashed line shows α = 0.05",
      x = "Path",
      y = "P-value"
    ) +
    scale_y_continuous(trans = "log10", 
                       breaks = c(0.001, 0.01, 0.05, 0.1, 0.5, 1),
                       labels = c("0.001", "0.01", "0.05", "0.1", "0.5", "1.0")) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "gray60")
    )
  
  return(p)
}

# Function 5: Create a heatmap showing proportion of significant paths
plot_significance_heatmap <- function(path_summary, avg_method = "conditional") {
  # Set column names based on averaging method
  if (avg_method == "conditional") {
    fill_col <- "mean_prop_sig_conditional"
    title_suffix <- "(Conditional Averaging)"
  } else {
    fill_col <- "mean_prop_sig_full"
    title_suffix <- "(Full Averaging)"
  }
  
  # Create the plot
  p <- ggplot(path_summary, aes(x = from, y = to, fill = .data[[fill_col]])) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = sprintf("%.2f\n(%d seeds)", .data[[fill_col]], n_seeds_with_path)),
              color = "white", size = 3, fontface = "bold") +
    scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred", 
                         midpoint = 0.5, name = "Mean Prop.\nSignificant") +
    labs(
      title = paste("Mean Proportion Significant", title_suffix),
      subtitle = "Numbers show: mean proportion (n_seeds_with_path)",
      x = "From Variable",
      y = "To Variable"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "gray60")
    )
  
  return(p)
}

# Function 6: Create coefficient magnitude heatmap
plot_coefficient_heatmap <- function(path_summary, avg_method = "conditional") {
  # Set column names based on averaging method
  if (avg_method == "conditional") {
    fill_col <- "mean_coef_conditional"
    title_suffix <- "(Conditional Averaging)"
  } else if (avg_method == "full") {
    fill_col <- "mean_coef_full"
    title_suffix <- "(Full Averaging)"
  } else {
    fill_col <- "mean_coefficient"
    title_suffix <- "(??)"
  }
  
  p <- ggplot(path_summary, aes(x = from, y = to, fill = .data[[fill_col]])) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = sprintf("%.3f", .data[[fill_col]])),
              color = ifelse(abs(path_summary[[fill_col]]) > 0.5, "white", "black"), 
              size = 3, fontface = "bold") +
    scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred", 
                         midpoint = 0, name = "Mean\nCoefficient") +
    labs(
      title = paste("Mean Coefficient Estimates", title_suffix),
      x = "From Variable",
      y = "To Variable"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 14, face = "bold")
    )
  
  return(p)
}

# Function to create top models bar plot 
plot_top_models <- function(model_frequencies, top_n = 10) {
  if (nrow(model_frequencies) > 0 && "frequency_in_sub2" %in% names(model_frequencies)) {
    top_n <- min(top_n, nrow(model_frequencies))
    top_models <- model_frequencies %>% 
      arrange(desc(frequency_in_sub2)) %>% 
      head(top_n)
    
    # Create a bar plot of top models
    top_models_plot <- ggplot(top_models, aes(x = reorder(model, frequency_in_sub2), y = frequency_in_sub2)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      coord_flip() +
      labs(
        title = paste("Top", top_n, "Models by Frequency in CICsub2"),
        x = "Model",
        y = "Frequency"
      ) +
      theme_minimal()
    
    return(top_models_plot)
  } else {
    return(NULL)
  }
}

# Improved aggregate_by_seed function with proper weighting
aggregate_by_seed <- function(detailed_models_df, cutoff = 2) {
  # Get all path coefficient columns
  coef_cols <- grep("_est$", names(detailed_models_df), value = TRUE)
  path_names <- gsub("_est$", "", coef_cols)
  
  # Get unique seeds
  all_seeds <- unique(detailed_models_df$seed)
  
  # Initialize results list
  seed_results <- list()
  
  # Process each seed separately
  for (seed in all_seeds) {
    seed_data <- detailed_models_df[detailed_models_df$seed == seed, ]
    
    # Calculate Akaike weights for this seed
    weight_info <- calculate_akaike_weights(seed_data$delta_CICc, cutoff = cutoff)
    valid_models <- seed_data[weight_info$valid_indices, ]
    weights <- weight_info$rel_weights
    
    # Initialize seed-level results
    seed_result <- list(
      seed = seed,
      nSpecies = seed_data$nSpecies[1],
      n_models_total = nrow(seed_data),
      n_models_used = length(weight_info$valid_indices),
      conditional = list(),
      full = list()
    )
    
    # Process each path
    for (i in seq_along(path_names)) {
      path_name <- path_names[i]
      coef_col <- coef_cols[i]
      se_col <- gsub("_est$", "_se", coef_col)
      p_col <- gsub("_est$", "_p", coef_col)
      
      # Extract coefficients and standard errors for valid models
      coefs <- valid_models[[coef_col]]
      ses <- valid_models[[se_col]]
      ps <- valid_models[[p_col]]
      
      # Convert "NA" strings to actual NA values
      coefs[coefs == "NA" | is.na(coefs)] <- NA
      ses[ses == "NA" | is.na(ses)] <- NA
      ps[ps == "NA" | is.na(ps)] <- NA
      
      # Convert to numeric
      coefs <- as.numeric(coefs)
      ses <- as.numeric(ses)
      ps <- as.numeric(ps)
      
      # CONDITIONAL AVERAGING: only average over models where path exists
      non_na_indices <- !is.na(coefs)
      if (sum(non_na_indices) > 0) {
        # Renormalize weights for conditional averaging
        cond_weights <- weights[non_na_indices]
        cond_weights <- cond_weights / sum(cond_weights)
        
        # Calculate weighted average (conditional)
        cond_coef <- stats::weighted.mean(coefs[non_na_indices], w = cond_weights)
        
        # For standard error, use the phylopath approach with variance averaging
        cond_vars <- ses[non_na_indices]^2
        # Weighted average of variances + variance of estimates
        avg_var <- stats::weighted.mean(cond_vars + coefs[non_na_indices]^2, w = cond_weights) - cond_coef^2
        cond_se <- sqrt(max(0, avg_var))
        
        # Average p-values (for reference)
        cond_p <- stats::weighted.mean(ps[non_na_indices], w = cond_weights, na.rm = TRUE)
        
        seed_result$conditional[[path_name]] <- list(
          est = cond_coef,
          se = cond_se,
          p = cond_p,
          n_models = sum(non_na_indices)
        )
      } else {
        seed_result$conditional[[path_name]] <- list(
          est = 0,
          se = 0,
          p = NA,
          n_models = 0
        )
      }
      
      # FULL AVERAGING: treat missing paths as 0
      coefs_full <- coefs
      ses_full <- ses
      coefs_full[is.na(coefs_full)] <- 0
      ses_full[is.na(ses_full)] <- 0
      
      # Calculate weighted average (full)
      full_coef <- stats::weighted.mean(coefs_full, w = weights)
      
      # For SE, use phylopath approach
      full_vars <- ses_full^2
      avg_var <- stats::weighted.mean(full_vars + coefs_full^2, w = weights) - full_coef^2
      full_se <- sqrt(max(0, avg_var))
      
      # Average p-values
      ps_full <- ps
      ps_full[is.na(ps_full)] <- 1 # Set missing p-values to 1 (non-significant)
      full_p <- stats::weighted.mean(ps_full, w = weights)
      
      seed_result$full[[path_name]] <- list(
        est = full_coef,
        se = full_se,
        p = full_p,
        n_models = length(weights)
      )
    }
    
    seed_results[[as.character(seed)]] <- seed_result
  }
  
  # Create summary across seeds
  summary_results <- summarize_across_seeds(seed_results, path_names)
  
  return(list(
    seed_level_results_list = seed_results,
    across_seeds_coefficient_summary = summary_results
  ))
}


# Function to summarize results across seeds
summarize_across_seeds <- function(seed_results, path_names) {
  # Initialize summary dataframes
  conditional_summary <- data.frame()
  full_summary <- data.frame()
  
  for (path_name in path_names) {
    # Extract estimates for this path across all seeds
    cond_ests <- sapply(seed_results, function(x) x$conditional[[path_name]]$est %||% 0)
    cond_ses <- sapply(seed_results, function(x) x$conditional[[path_name]]$se %||% 0)
    cond_ps <- sapply(seed_results, function(x) x$conditional[[path_name]]$p %||% 1)  # Added this line
    cond_n_models <- sapply(seed_results, function(x) x$conditional[[path_name]]$n_models %||% 0)
    
    full_ests <- sapply(seed_results, function(x) x$full[[path_name]]$est %||% 0)
    full_ses <- sapply(seed_results, function(x) x$full[[path_name]]$se %||% 0)
    full_ps <- sapply(seed_results, function(x) x$full[[path_name]]$p %||% 1)  # Added this line
    full_n_models <- sapply(seed_results, function(x) x$full[[path_name]]$n_models %||% 0)
    
    # Calculate summary statistics for conditional averaging
    cond_summary_row <- data.frame(
      path = path_name,
      from = sub("_to_.*$", "", path_name),
      to = sub("^.*_to_", "", path_name),
      method = "conditional",
      mean_coef_conditional = mean(cond_ests),  # FIXED: was mean_est
      sd_est = sd(cond_ests),
      min_est = min(cond_ests),
      max_est = max(cond_ests),
      mean_se = mean(cond_ses),
      mean_n_models = mean(cond_n_models),
      n_seeds_with_path = sum(cond_n_models > 0),
      prop_seeds_with_path = sum(cond_n_models > 0) / length(seed_results),
      mean_prop_sig_conditional = mean(cond_ps < 0.05, na.rm = TRUE),  # ADDED this line
      stringsAsFactors = FALSE
    )
    
    # Calculate summary statistics for full averaging
    full_summary_row <- data.frame(
      path = path_name,
      from = sub("_to_.*$", "", path_name),
      to = sub("^.*_to_", "", path_name),
      method = "full",
      mean_coef_full = mean(full_ests),  # FIXED: was mean_est
      sd_est = sd(full_ests),
      min_est = min(full_ests),
      max_est = max(full_ests),
      mean_se = mean(full_ses),
      mean_n_models = mean(full_n_models),
      n_seeds_with_path = sum(cond_n_models > 0), # Same as conditional
      prop_seeds_with_path = sum(cond_n_models > 0) / length(seed_results),
      mean_prop_sig_full = mean(full_ps < 0.05, na.rm = TRUE),  # ADDED this line
      stringsAsFactors = FALSE
    )
    
    conditional_summary <- rbind(conditional_summary, cond_summary_row)
    full_summary <- rbind(full_summary, full_summary_row)
  }
  
  # Sort by most common paths first
  conditional_summary <- conditional_summary %>%
    arrange(desc(prop_seeds_with_path), desc(abs(mean_coef_conditional)))
  
  full_summary <- full_summary %>%
    arrange(desc(prop_seeds_with_path), desc(abs(mean_coef_full)))

  
  return(list(
    conditionalAverage_coefficient_acrossSeeds = conditional_summary,
    fullAverage_coefficient_acrossSeeds = full_summary
  ))
}

# Function to convert seed-level results from aggregate_by_seed() to dataframe format
convert_seed_results_to_dataframe <- function(seed_results, method = "conditional") {
  # Validate method
  if (!method %in% c("conditional", "full")) {
    stop("Method must be either 'conditional' or 'full'")
  }
  
  # Initialize empty list to store rows
  result_rows <- list()
  
  # Extract all unique path names from any seed
  all_paths <- unique(unlist(lapply(seed_results, function(x) names(x[[method]]))))
  
  # Iterate through each seed
  for (seed_name in names(seed_results)) {
    seed_data <- seed_results[[seed_name]]
    
    # Get basic seed information
    seed_info <- list(
      seed = seed_data$seed,
      nSpecies = seed_data$nSpecies,
      n_models_total = seed_data$n_models_total,
      n_models_used = seed_data$n_models_used
    )
    
    # For each path, create a row
    for (path_name in all_paths) {
      path_result <- seed_data[[method]][[path_name]]
      
      # Create a row for this seed-path combination
      if (!is.null(path_result)) {
        row <- c(
          seed_info,
          list(
            path = path_name,
            from = sub("_to_.*$", "", path_name),
            to = sub("^.*_to_", "", path_name),
            mean_coefficient = path_result$est,
            std_error = path_result$se,
            p_value = path_result$p,
            significant = if(!is.na(path_result$p)) path_result$p < 0.05 else NA,
            n_models_with_path = path_result$n_models
          )
        )
      } else {
        # Path doesn't exist for this seed
        row <- c(
          seed_info,
          list(
            path = path_name,
            from = sub("_to_.*$", "", path_name),
            to = sub("^.*_to_", "", path_name),
            mean_coefficient = if(method == "full") 0 else NA,
            std_error = if(method == "full") 0 else NA,
            p_value = if(method == "full") 1 else NA,
            significant = if(method == "full") FALSE else NA,
            n_models_with_path = 0
          )
        )
      }
      
      result_rows[[length(result_rows) + 1]] <- row
    }
  }
  
  # Convert list of rows to dataframe
  if (length(result_rows) > 0) {
    # Convert to dataframe
    result_df <- do.call(rbind, lapply(result_rows, function(x) {
      # Convert each row to a dataframe row
      data.frame(
        seed = as.integer(x$seed),
        nSpecies = as.integer(x$nSpecies),
        n_models_total = as.integer(x$n_models_total),
        n_models_used = as.integer(x$n_models_used),
        path = as.character(x$path),
        from = as.character(x$from),
        to = as.character(x$to),
        mean_coefficient = as.numeric(x$mean_coefficient),
        std_error = as.numeric(x$std_error),
        p_value = as.numeric(x$p_value),
        significant = as.logical(x$significant),
        n_models_with_path = as.integer(x$n_models_with_path),
        method = method,
        stringsAsFactors = FALSE
      )
    }))
    
    # Clean up the dataframe
    result_df <- result_df %>%
      arrange(seed, path)
    
    return(result_df)
  } else {
    return(data.frame())
  }
}

# Wrapper function to get both conditional and full results as dataframes
convert_all_seed_results_to_dataframes <- function(aggregate_results) {
  seed_results <- aggregate_results$seed_level_results_list
  
  # Convert both conditional and full results
  conditional_df <- convert_seed_results_to_dataframe(seed_results, method = "conditional")
  full_df <- convert_seed_results_to_dataframe(seed_results, method = "full")
  
  # Combine both methods into one dataframe
  combined_df <- rbind(conditional_df, full_df)
  
  return(list(
    conditionalAverage_coefficient_perSeed = conditional_df,
    fullAverage_coefficient_perSeed = full_df,
    combined_coefficients_perSeed = combined_df
  ))
}


# Updated create_all_plots() function to work with new aggregate_by_seed()
create_all_plots <- function(detailed_models_df, save_plots = TRUE, file_prefix = "phylopath_analysis", cutoff = 2) {
  # Extract path data using the original function (for model-level analysis)
  cat("Creating path summary...\n")
  path_results <- create_path_summary(detailed_models_df)
  
  if (is.null(path_results)) {
    cat("No path data found!\n")
    return(NULL)
  }
  
  path_data <- path_results$path_coefficients
  path_summary <- path_results$path_summary
  
  # Create model summary and top models plot
  cat("Creating model summary...\n")
  model_summary <- create_model_summary(detailed_models_df)
  
  # Calculate model frequencies from detailed_models_df
  # slightly redundant since this is basically the same as one column of model_summary
  model_freq_df <- detailed_models_df %>%
    group_by(model) %>%
    summarize(
      frequency_in_sub2 = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(frequency_in_sub2))
  
  top_models_plot <- plot_top_models(model_freq_df)
  
  # aggregate_by_seed function with weighting
  cat("Creating seed-level aggregation with proper weighting...\n")
  aggregate_results <- aggregate_by_seed(detailed_models_df, cutoff = cutoff)
  
  # Convert seed-level results to dataframe format for plotting
  seed_dataframes <- convert_all_seed_results_to_dataframes(aggregate_results)
  seed_level_conditional <- seed_dataframes$conditionalAverage_coefficient_perSeed
  seed_level_full <- seed_dataframes$fullAverage_coefficient_perSeed
  seed_level_combined <- seed_dataframes$combined_coefficients_perSeed
  
  # Extract the summary results
  conditionalAverage_coefficient_acrossSeeds <- path_summary_conditional <- aggregate_results$across_seeds_coefficient_summary$conditionalAverage_coefficient_acrossSeeds
  fullAverage_coefficient_acrossSeeds <- aggregate_results$across_seeds_coefficient_summary$fullAverage_coefficient_acrossSeeds

  
  # Create plots using original (model-level) data
  cat("Creating coefficient distribution plots (model-level)...\n")
  coef_violin_models <- plot_coefficient_distributions(path_data, show_seeds = FALSE, plot_type = "violin")
  coef_violin_models <- coef_violin_models + labs(subtitle = "Based on all models (some seeds overrepresented)")
  
  cat("Creating p-value distribution plots (model-level)...\n")
  pval_violin_models <- plot_pvalue_distributions(path_data, show_seeds = FALSE, plot_type = "violin")
  pval_violin_models <- pval_violin_models + labs(subtitle = "Based on all models (some seeds overrepresented)")
  
  # Create plots using seed-aggregated data (conditional averaging)
  cat("Creating coefficient distribution plots (conditional averaging)...\n")
  coef_violin_conditional <- plot_coefficient_distributions(seed_level_conditional, show_seeds = FALSE, plot_type = "violin")
  coef_violin_conditional <- coef_violin_conditional + labs(subtitle = "Conditional averaging: weighted within seeds, equal weight per seed")
  
  # Create plots using seed-aggregated data (full averaging)
  cat("Creating coefficient distribution plots (full averaging)...\n")
  coef_violin_full <- plot_coefficient_distributions(seed_level_full, show_seeds = FALSE, plot_type = "violin")
  coef_violin_full <- coef_violin_full + labs(subtitle = "Full averaging: weighted within seeds, missing paths = 0")
  
  cat("Creating p-value distribution plots (seed-level conditional)...\n")
  pval_violin_conditional <- plot_pvalue_distributions(seed_level_conditional, show_seeds = FALSE, plot_type = "violin")
  pval_violin_conditional <- pval_violin_conditional + labs(subtitle = "Conditional averaging: weighted p-values within seeds")
  
  cat("Creating p-value distribution plots (seed-level full)...\n")
  pval_violin_full <- plot_pvalue_distributions(seed_level_full, show_seeds = FALSE, plot_type = "violin")
  pval_violin_full <- pval_violin_full + labs(subtitle = "Full averaging: weighted p-values within seeds")
  
  # Create heatmaps using the summary results
  cat("Creating heatmaps...\n")
  
  # Model-level heatmaps (original)
#  sig_heatmap_models <- plot_significance_heatmap(path_summary) + 
#    labs(title = "Proportion of Significant Paths (Averaged across all models and seeds - some seeds overrepresented)")
  coef_heatmap_models <- plot_coefficient_heatmap(path_summary, avg_method = "none") + 
    labs(title = "Mean Coefficient Estimates (Averaged across all models and seeds - some seeds overrepresented)")
  
  # Seed-level heatmaps (conditional averaging)
  sig_heatmap_conditional <- plot_significance_heatmap(  conditionalAverage_coefficient_acrossSeeds, avg_method = "conditional") + 
    labs(title = "Proportion of Significant Paths (Conditional Averaging)")
  coef_heatmap_conditional <- plot_coefficient_heatmap(  conditionalAverage_coefficient_acrossSeeds, avg_method = "conditional") + 
    labs(title = "Mean Coefficient Estimates (Conditional Averaging)")
  
  # Seed-level heatmaps (full averaging)
  sig_heatmap_full <- plot_significance_heatmap(fullAverage_coefficient_acrossSeeds, avg_method = "full") + 
    labs(title = "Proportion of Significant Paths (Full Averaging)")
  coef_heatmap_full <- plot_coefficient_heatmap(fullAverage_coefficient_acrossSeeds, avg_method = "full") + 
    labs(title = "Mean Coefficient Estimates (Full Averaging)")
  
  # Save plots if requested
  if (save_plots) {
    cat("Saving plots...\n")
    
    pdf(paste0(file_prefix, "_all_plots_", Sys.Date(), ".pdf"), width = 12, height = 8)
    
    # 1. Top models plot first
    if (!is.null(top_models_plot)) {
      print(top_models_plot)
    }
    
    # 2. Coefficient distributions comparison
    print(coef_violin_models)
    print(coef_violin_conditional)
    print(coef_violin_full)
    
    # 3. P-value distributions comparison
    print(pval_violin_models)
    print(pval_violin_conditional)
    print(pval_violin_full)
    
    # 4. Model-level heatmaps (for reference)
    #print(sig_heatmap_models)
    print(coef_heatmap_models)
    
    # 5. Seed-level heatmaps - conditional averaging (main results)
    print(sig_heatmap_conditional)
    print(coef_heatmap_conditional)
    
    # 6. Seed-level heatmaps - full averaging (conservative estimates)
    print(sig_heatmap_full)
    print(coef_heatmap_full)
    
    dev.off()
    
    # # Save all summaries to CSV
    # write.csv(path_summary_conditional, paste0(file_prefix, "_path_summary_conditional_", Sys.Date(), ".csv"), row.names = FALSE)
    # write.csv(path_summary_full, paste0(file_prefix, "_path_summary_full_", Sys.Date(), ".csv"), row.names = FALSE)
    # write.csv(path_summary_combined, paste0(file_prefix, "_path_summary_combined_", Sys.Date(), ".csv"), row.names = FALSE)
    # write.csv(seed_level_conditional, paste0(file_prefix, "_seed_level_conditional_", Sys.Date(), ".csv"), row.names = FALSE)
    # write.csv(seed_level_full, paste0(file_prefix, "_seed_level_full_", Sys.Date(), ".csv"), row.names = FALSE)
    # write.csv(model_summary, paste0(file_prefix, "_model_summary_", Sys.Date(), ".csv"), row.names = FALSE)
    # write.csv(path_data, paste0(file_prefix, "_path_coefficients_", Sys.Date(), ".csv"), row.names = FALSE)
    # 
    # cat("All files saved successfully!\n")
  }
  
  # Return all results
  return(list(
    # Raw data
    path_data = path_data,
    model_summary = model_summary,
    model_frequencies = model_freq_df,
    
    # Aggregated results
    aggregate_results = aggregate_results,
    seed_level_conditional = seed_level_conditional,
    seed_level_full = seed_level_full,
    
    # Path summaries (across seeds)
    path_summary_conditional = path_summary_conditional,
    path_summary_full = path_summary_full,
    path_summary_combined = path_summary_combined,
    
    # Plots
    plots = list(
      top_models = top_models_plot,
      coef_violin_models = coef_violin_models,
      coef_violin_conditional = coef_violin_conditional,
      coef_violin_full = coef_violin_full,
      pval_violin_models = pval_violin_models,
      pval_violin_conditional = pval_violin_conditional,
      pval_violin_full = pval_violin_full,
      sig_heatmap_models = sig_heatmap_models,
      coef_heatmap_models = coef_heatmap_models,
      sig_heatmap_conditional = sig_heatmap_conditional,
      coef_heatmap_conditional = coef_heatmap_conditional,
      sig_heatmap_full = sig_heatmap_full,
      coef_heatmap_full = coef_heatmap_full
    )
  ))
}



# Function to create coefficient violin plot for downsampled results
plot_downsampled_coefficients <- function(seed_level_data, min_occurrences = 10) {
  # Clean path names
  seed_level_data$path_clean <- gsub("_to_", " → ", seed_level_data$path)
  seed_level_data$path_clean <- gsub("FemaleSong_Agg01", "Female Song", seed_level_data$path_clean)
  seed_level_data$path_clean <- gsub("HighConfidence_Coop", "Cooperation", seed_level_data$path_clean)
  seed_level_data$path_clean <- gsub("TerritorialityWeakVsStrong", "Territoriality", seed_level_data$path_clean)
  seed_level_data$path_clean <- gsub("logMass_AVONET", "Body Mass", seed_level_data$path_clean)
  
  # Filter to paths that appear frequently
  path_counts <- seed_level_data %>%
    group_by(path) %>%
    summarise(n = n()) %>%
    filter(n >= min_occurrences)
  
  seed_level_plot <- seed_level_data %>%
    filter(path %in% path_counts$path)
  
  ggplot(seed_level_plot, aes(x = path_clean, y = mean_coefficient)) +
    geom_violin(fill = "lightblue", alpha = 0.7) +
    geom_boxplot(width = 0.1, alpha = 0.8, outlier.shape = NA) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    coord_flip() +
    labs(
      title = "Distribution of Path Coefficients Across Iterations",
      subtitle = paste("Based on", length(unique(seed_level_plot$seed)), 
                       "downsampling iterations (conditional averaging)"),
      x = NULL,
      y = "Mean Coefficient"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(size = 10)
    )
}

# Function to create model frequency plot (downsampled runs)
plot_model_frequencies <- function(model_frequencies, n_iterations, top_n = 15) {
  model_freq <- model_frequencies %>%
    arrange(desc(frequency_in_sub2)) %>%
    head(top_n) %>%
    mutate(prop_sub2 = frequency_in_sub2 / n_iterations)
  
  ggplot(model_freq, 
         aes(x = reorder(model_name, frequency_in_sub2),
             y = frequency_in_sub2)) +
    geom_bar(stat = "identity", fill = "#2E86AB", alpha = 0.8) +
    geom_text(aes(label = paste0(frequency_in_sub2, " (", 
                                 round(prop_sub2 * 100, 1), "%)")),
              hjust = -0.1, size = 3) +
    coord_flip() +
    labs(
      title = "Top Models Across Downsampling Iterations",
      subtitle = paste("Models appearing in Δ CICc < 2 set across", n_iterations, "iterations"),
      x = NULL,
      y = "Frequency in Top Model Set"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(size = 9)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
}

#' Create dimorphism distribution plot showing bias correction
#'
#' @param dfIn_phylo Original dataframe with all species
#' @param dimorphism_downsample Output from downsample_dimorphism_bias()
#' @param dim_info List with dimorphism variable information (col, label, description)
#' @param output_dir Directory to save plots
#' @param prefix Output filename prefix
#' @param save_plot Whether to save the plot to file
#' @param show_iterations How many iterations to show (default 20 for visibility)
#' @return ggplot object
create_downsampled_dimorphism_distribution_plot <- function(dfIn_phylo,
                                                dimorphism_downsample,
                                                dim_info,
                                                output_dir = NULL,
                                                prefix = NULL,
                                                save_plot = TRUE,
                                                show_iterations = 20) {
  
  library(ggplot2)
  library(dplyr)
  
  # Calculate the mean dimorphism after correction from all iterations
  remaining_means <- sapply(dimorphism_downsample$iterations, function(iter) {
    mean(iter$remaining_data[[dim_info$col]][!is.na(iter$remaining_data$FemaleSong_Agg01)], 
         na.rm = TRUE)
  })
  mean_dimorphism_after <- mean(remaining_means, na.rm = TRUE)
  sd_dimorphism_after <- sd(remaining_means, na.rm = TRUE)
  
  # Base data for original distributions
  plot_data <- rbind(
    data.frame(
      group = "No FS Data",
      dimorphism = dfIn_phylo[[dim_info$col]][is.na(dfIn_phylo$FemaleSong_Agg01)],
      iteration = NA
    ),
    data.frame(
      group = "Has FS Data (Original)",
      dimorphism = dfIn_phylo[[dim_info$col]][!is.na(dfIn_phylo$FemaleSong_Agg01)],
      iteration = NA
    )
  )
  
  # Add multiple iterations of corrected data
  n_iter_to_show <- min(show_iterations, length(dimorphism_downsample$iterations))
  iterations_to_plot <- sample(1:length(dimorphism_downsample$iterations), n_iter_to_show)
  
  for (i in iterations_to_plot) {
    iter_data <- dimorphism_downsample$iterations[[i]]$remaining_data
    plot_data <- rbind(plot_data,
                       data.frame(
                         group = "Has FS Data (Corrected)",
                         dimorphism = iter_data[[dim_info$col]][!is.na(iter_data$FemaleSong_Agg01)],
                         iteration = i
                       )
    )
  }
  
  # Remove NA values
  plot_data <- plot_data[!is.na(plot_data$dimorphism), ]
  
  # Set factor levels to control order
  plot_data$group <- factor(plot_data$group, 
                            levels = c("No FS Data", "Has FS Data (Original)", "Has FS Data (Corrected)"))
  
  # Create density plot with transparency for iterations
  p <- ggplot(plot_data, aes(x = dimorphism, fill = group)) +
    # Original distributions
    geom_density(data = subset(plot_data, group != "Has FS Data (Corrected)"), 
                 alpha = 0.6) +
    # Corrected distributions (multiple iterations with lower alpha)
    geom_density(data = subset(plot_data, group == "Has FS Data (Corrected)"), 
                 aes(group = iteration), 
                 alpha = 0.1, color = NA) +
    # Average corrected distribution
    geom_density(data = subset(plot_data, group == "Has FS Data (Corrected)"), 
                 alpha = 0, color = "#4DAF4A", size = 1.5) +
    # Reference lines
    geom_vline(aes(xintercept = target_mean, linetype = "Target mean"), 
               data = data.frame(target_mean = dimorphism_downsample$target_stats$target_mean),
               color = "red", size = 1) +
    geom_vline(aes(xintercept = current_mean, linetype = "Original mean"), 
               data = data.frame(current_mean = dimorphism_downsample$current_stats$current_mean),
               color = "blue", size = 1) +
    geom_vline(aes(xintercept = corrected_mean, linetype = "Corrected mean"), 
               data = data.frame(corrected_mean = mean_dimorphism_after),
               color = "green", size = 1) +
    scale_fill_manual(values = c("No FS Data" = "#E41A1C", 
                                 "Has FS Data (Original)" = "#377EB8", 
                                 "Has FS Data (Corrected)" = "#4DAF4A")) +
    scale_linetype_manual(name = "Reference lines",
                          values = c("Target mean" = "dashed", 
                                     "Original mean" = "dotted", 
                                     "Corrected mean" = "solid")) +
    labs(
      title = paste(dim_info$description, "Distribution Before and After Bias Correction"),
      subtitle = paste("Removed", dimorphism_downsample$n_to_remove, 
                       "species using propensity weighting (", 
                       dimorphism_downsample$n_to_remove, "/", 
                       dimorphism_downsample$current_stats$n_species,
                       " = ", round(100 * dimorphism_downsample$n_to_remove / 
                                      dimorphism_downsample$current_stats$n_species, 1), "%)\n",
                       "Showing ", n_iter_to_show, " of ", length(dimorphism_downsample$iterations), 
                       " iterations", sep = ""),
      x = paste(dim_info$description, "(log-transformed)"),
      y = "Density",
      fill = "Group"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  
  # Add text annotations with means
  y_position <- max(density(plot_data$dimorphism)$y) * 0.9
  p <- p + 
    annotate("text", x = dimorphism_downsample$target_stats$target_mean, 
             y = y_position, label = paste("Target:", round(dimorphism_downsample$target_stats$target_mean, 3)),
             color = "red", hjust = -0.1, size = 3) +
    annotate("text", x = dimorphism_downsample$current_stats$current_mean, 
             y = y_position * 0.85, label = paste("Original:", round(dimorphism_downsample$current_stats$current_mean, 3)),
             color = "blue", hjust = 1.1, size = 3) +
    annotate("text", x = mean_dimorphism_after, 
             y = y_position * 0.95, label = paste("Mean corrected:", round(mean_dimorphism_after, 3)),
             color = "green", hjust = -0.1, size = 3)
  
  # Save if requested
  if (save_plot && !is.null(output_dir) && !is.null(prefix)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    filename <- file.path(output_dir, paste0(prefix, "_distribution_correction.png"))
    ggsave(
      filename = filename,
      plot = p,
      width = 10,
      height = 6,
      dpi = 300
    )
    cat("Dimorphism distribution plot saved to:", filename, "\n")
  }
  
  # Calculate n_species for each group
  n_no_data <- sum(is.na(dfIn_phylo$FemaleSong_Agg01) & !is.na(dfIn_phylo[[dim_info$col]]))
  n_has_data <- sum(!is.na(dfIn_phylo$FemaleSong_Agg01) & !is.na(dfIn_phylo[[dim_info$col]]))
  n_corrected <- n_has_data - dimorphism_downsample$n_to_remove
  
  # Create summary statistics
  summary_stats <- data.frame(
    Group = c("No FS Data", "Has FS Data (Original)", "Has FS Data (Corrected)"),
    N = c(n_no_data, n_has_data, n_corrected),
    Mean = round(c(dimorphism_downsample$target_stats$target_mean,
                   dimorphism_downsample$current_stats$current_mean,
                   mean_dimorphism_after), 4),
    SD = round(c(dimorphism_downsample$target_stats$target_sd,
                 dimorphism_downsample$current_stats$current_sd,
                 sd_dimorphism_after), 4),
    Median = round(c(dimorphism_downsample$target_stats$target_median,
                     median(dfIn_phylo[[dim_info$col]][!is.na(dfIn_phylo$FemaleSong_Agg01)], na.rm = TRUE),
                     median(remaining_means)), 4)
  )
  
  # Add convergence info
  convergence <- abs(mean_dimorphism_after - dimorphism_downsample$target_stats$target_mean)
  cat("Mean dimorphism after correction:", round(mean_dimorphism_after, 3), "\n")
  cat("SD of means across iterations:", round(sd_dimorphism_after, 4), "\n")
  cat("Convergence (distance from target):", round(convergence, 4), "\n")
  
  if (save_plot && !is.null(output_dir) && !is.null(prefix)) {
    write.csv(
      summary_stats,
      file = file.path(output_dir, paste0(prefix, "_summary_statistics.csv")),
      row.names = FALSE
    )
  }
  
  # Return both plot and summary stats
  return(list(
    plot = p,
    summary_stats = summary_stats,
    mean_after = mean_dimorphism_after,
    convergence = convergence
  ))
}


#' Create violin plot showing dimorphism distributions across iterations
#'
#' @param dfIn_phylo Original dataframe with all species
#' @param dimorphism_downsample Output from downsample_dimorphism_bias()
#' @param dim_info List with dimorphism variable information (col, label, description)
#' @param output_dir Directory to save plots
#' @param prefix Output filename prefix
#' @param save_plot Whether to save the plot to file
#' @param show_iterations How many iterations to show (default 20)
#' @return ggplot object
create_dimorphism_downsamples_violin_plot <- function(dfIn_phylo,
                                          dimorphism_downsample,
                                          dim_info,
                                          output_dir = NULL,
                                          prefix = NULL,
                                          save_plot = TRUE,
                                          show_iterations = 20) {
  
  library(ggplot2)
  library(dplyr)
  
  # Prepare data for violin plot
  # Original distributions
  original_data <- rbind(
    data.frame(
      group = "No FS Data",
      category = "Original",
      iteration = "Original",
      dimorphism = dfIn_phylo[[dim_info$col]][is.na(dfIn_phylo$FemaleSong_Agg01)]
    ),
    data.frame(
      group = "Has FS Data",
      category = "Original", 
      iteration = "Original",
      dimorphism = dfIn_phylo[[dim_info$col]][!is.na(dfIn_phylo$FemaleSong_Agg01)]
    )
  )
  
  # Remove NAs
  original_data <- original_data[!is.na(original_data$dimorphism), ]
  
  # Corrected iterations
  n_iter_to_show <- min(show_iterations, length(dimorphism_downsample$iterations))
  iterations_to_plot <- 1:n_iter_to_show
  
  iteration_data <- data.frame()
  for (i in iterations_to_plot) {
    iter_result <- dimorphism_downsample$iterations[[i]]$remaining_data
    temp_data <- data.frame(
      group = paste0("Iteration ", i),
      category = "Corrected",
      iteration = as.character(i),
      dimorphism = iter_result[[dim_info$col]][!is.na(iter_result$FemaleSong_Agg01)]
    )
    iteration_data <- rbind(iteration_data, temp_data[!is.na(temp_data$dimorphism), ])
  }
  
  # Combine all data
  plot_data <- rbind(original_data, iteration_data)
  
  # Set factor levels to control order
  group_levels <- c("No FS Data", "Has FS Data", 
                    paste0("Iteration ", 1:n_iter_to_show))
  plot_data$group <- factor(plot_data$group, levels = group_levels)
  
  # Calculate means for reference lines
  mean_no_data <- dimorphism_downsample$target_stats$target_mean
  mean_has_data <- dimorphism_downsample$current_stats$current_mean
  
  iteration_means <- sapply(1:n_iter_to_show, function(i) {
    iter_result <- dimorphism_downsample$iterations[[i]]$remaining_data
    mean(iter_result[[dim_info$col]][!is.na(iter_result$FemaleSong_Agg01)], na.rm = TRUE)
  })
  
  # Create violin plot
  pViolin <- ggplot(plot_data, aes(x = group, y = dimorphism, fill = category)) +
    geom_violin(trim = FALSE, scale = "width", alpha = 0.8) +
    geom_boxplot(width = 0.1, alpha = 0.6, outlier.shape = NA, 
                 position = position_dodge(width = 0.9)) +
    
    # Add vertical separator line
    geom_vline(xintercept = 2.5, linetype = "dashed", color = "gray50", size = 1) +
    
    # Add horizontal reference lines
    geom_hline(yintercept = mean_no_data, linetype = "dashed", color = "red", alpha = 0.7) +
    geom_hline(yintercept = mean_has_data, linetype = "dotted", color = "blue", alpha = 0.7) +
    
    # Add mean points with clear labeling
    geom_point(data = data.frame(
      x = c(1, 2, 3:(n_iter_to_show + 2)),
      y = c(mean_no_data, mean_has_data, iteration_means),
      category = c("Original", "Original", rep("Corrected", n_iter_to_show))
    ), aes(x = x, y = y), 
    color = "black", size = 3, shape = 18) +
    
    scale_fill_manual(values = c("Original" = "#377EB8", "Corrected" = "#4DAF4A"),
                      name = "Data Type") +
    
    labs(
      title = paste(dim_info$description, "Distribution: Bias Correction Across Iterations"),
      subtitle = paste("Removed", dimorphism_downsample$n_to_remove, 
                       "species per iteration using propensity weighting\n",
                       "Target mean (red dashed):", round(mean_no_data, 3),
                       " | Original mean (blue dotted):", round(mean_has_data, 3),
                       " | Mean of corrected means:", round(mean(iteration_means), 3),
                       "\nBlack diamonds show group means"),
      x = "",
      y = paste(dim_info$description, "(log-transformed)")
    ) +
    
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      legend.position = "bottom",
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10),
      panel.grid.major.x = element_blank()
    ) +
    
    # Add annotations
    annotate("text", x = 1.5, y = max(plot_data$dimorphism) * 0.95,
             label = "Original Data", fontface = "bold", size = 4) +
    annotate("text", x = (n_iter_to_show + 2) / 2 + 1.5, y = max(plot_data$dimorphism) * 0.95,
             label = "Bias-Corrected Iterations of Have-FS-Data", fontface = "bold", size = 4)
  
  # Save if requested
  if (save_plot && !is.null(output_dir) && !is.null(prefix)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    filename <- file.path(output_dir, paste0(prefix, "_violin_iterations.png"))
    ggsave(
      filename = filename,
      plot = pViolin,
      width = 14,  # Increased width to accommodate longer labels
      height = 7,
      dpi = 300
    )
    cat("Dimorphism violin plot saved to:", filename, "\n")
  }
  
  # Create summary statistics for iterations
  iteration_stats <- data.frame(
    Iteration = 1:n_iter_to_show,
    Mean = round(iteration_means, 4),
    Distance_from_target = round(abs(iteration_means - mean_no_data), 4)
  )
  
  if (save_plot && !is.null(output_dir) && !is.null(prefix)) {
    write.csv(
      iteration_stats,
      file = file.path(output_dir, paste0(prefix, "_iteration_means.csv")),
      row.names = FALSE
    )
  }
  
  return(list(
    plot = pViolin,
    iteration_stats = iteration_stats,
    overall_mean = mean(iteration_means),
    overall_sd = sd(iteration_means)
  ))
}

#### Wrapper Plotting functions from phylopath_plotting_comprehensive_fixed.R ----
#' Create all plots for non-downsampled phylopath runs
#'
#' @param phylopath_output Output from run_CB_FS_Terr_phylopath()
#' @param output_prefix Prefix for output files - should include all variables
#' @param save_png Save plots as PNG
#' @param save_pdf Save plots as PDF
create_nondownsampled_plots <- function(phylopath_output, 
                                        output_prefix = "phylopath",
                                        output_dir = "Outputs/PhylopathPlots",
                                        save_png = TRUE,
                                        save_pdf = FALSE) {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Extract components
  result <- phylopath_output$result
  summary_result <- summary(result)
  var_map <- phylopath_output$var_map
  nSpecies <- phylopath_output$nSpecies
  
  # Get variable names for file naming
  var_names <- names(var_map)
  var_string <- paste(unlist(var_map), collapse = " ")
  
  # Create file name base with all variables and species count
  file_base <- paste0("phylopath ", var_string, " ", nSpecies, "species")
  
  # Set up layout positions for DAGs with better spacing
  if (!is.null(var_map)) {
    # Create sensible positions based on variable names
    var_names <- names(var_map)
    n_vars <- length(var_names)
    
    if (n_vars == 4) {
      # Standard 4-variable layout with good spacing
      phylopath_map_positions <- data.frame(
        name = unlist(var_map),
        x = c(2, 8, 5, 5),  # Original positions
        y = c(9, 9, 1, 5),  # Original positions
        stringsAsFactors = FALSE
      )
    } else {
      # Circular layout for other numbers
      angles <- seq(0, 2*pi, length.out = n_vars + 1)[1:n_vars]
      phylopath_map_positions <- data.frame(
        name = unlist(var_map),
        x = 5 + 3 * cos(angles),  # Smaller radius to keep away from edges
        y = 5 + 3 * sin(angles),
        stringsAsFactors = FALSE
      )
    }
  } else {
    phylopath_map_positions <- NULL
  }
  
  # 1. Best model DAG
  cat("Creating best model DAG...\n")
  best_model <- best(result)
  p_best <- plot(best_model, 
                 manual_layout = phylopath_map_positions,
                 text_size = 3,  # Reduced from 4
                 box_x = 22, box_y = 16) +  # Further increased box size
    ggtitle("Best Model (Lowest CICc)") +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          plot.margin = margin(40, 40, 40, 40)) +  # Larger margins
    coord_cartesian(xlim = c(-1, 11), ylim = c(-1, 11), expand = TRUE)  # Expand plot area
  
  # 2. Fully-averaged best models DAG
  cat("Creating fully-averaged models DAG...\n")
  full_avg <- average(result, cut_off = 2, avg_method = "full")
  p_full_avg <- plot(full_avg,
                     manual_layout = phylopath_map_positions,
                     text_size = 3,  # Reduced from 4
                     box_x = 22, box_y = 16) +  # Further increased box size
    ggtitle("Fully-Averaged Models (Δ CICc < 2)") +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          plot.margin = margin(40, 40, 40, 40)) +  # Larger margins
    coord_cartesian(xlim = c(-1, 11), ylim = c(-1, 11), expand = TRUE)  # Expand plot area
  
  # 3. Conditional-averaged best models DAG
  cat("Creating conditional-averaged models DAG...\n")
  cond_avg <- average(result, cut_off = 2, avg_method = "conditional")
  p_cond_avg <- plot(cond_avg,
                     manual_layout = phylopath_map_positions,
                     text_size = 3,  # Reduced from 4
                     box_x = 22, box_y = 16) +  # Further increased box size
    ggtitle("Conditional-Averaged Models (Δ CICc < 2)") +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          plot.margin = margin(40, 40, 40, 40)) +  # Larger margins
    coord_cartesian(xlim = c(-1, 11), ylim = c(-1, 11), expand = TRUE)  # Expand plot area
  
  # 4. Bar plot of CICc values - KEEP ORIGINAL MODEL NAMES
  cat("Creating CICc bar plot...\n")
  # Prepare data - DON'T clean the model names
  model_data <- summary_result %>%
    mutate(
      is_top = delta_CICc < 2
    ) %>%
    arrange(CICc)
  
  # Limit to top 20 models for clarity
  n_models <- min(20, nrow(model_data))
  model_data_plot <- model_data[1:n_models, ]
  
  p_cicbar <- ggplot(model_data_plot, 
                     aes(x = reorder(model, -CICc),  # Use original model names
                         y = CICc,
                         fill = is_top)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_text(aes(label = round(delta_CICc, 2)), 
              hjust = -0.2, size = 3) +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "#2E86AB", "FALSE" = "#A7A7A7"),
                      labels = c("TRUE" = "Δ CICc < 2", "FALSE" = "Δ CICc ≥ 2"),
                      name = "Model Set") +
    labs(
      title = "Model Comparison by CICc",
      subtitle = paste("Showing top", n_models, "models. Values show Δ CICc"),
      x = NULL,
      y = "CICc"
    ) +
    theme_cowplot(12) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(size = 9),
      legend.position = "bottom"
    )
  
  # 5. Plot of all tested models
  cat("Creating summary plot of all models...\n")
  p_summary <- plot(summary_result) +
    ggtitle(paste("All Tested Models (n =", nrow(summary_result), ")")) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.y = element_text(size = 8)
    )
  
  # Combine DAGs into one figure
  dag_combined <- plot_grid(
    p_best, p_full_avg, p_cond_avg,
    ncol = 3,
    labels = c("A", "B", "C"),
    label_size = 16
  )
  
  # Save plots with proper file names
  if (save_png) {
    # Individual plots - increased sizes to prevent squishing
    ggsave(file.path(output_dir, paste0(file_base, " best_model.png")), 
           p_best, width = 15, height = 10, dpi = 300)
    ggsave(file.path(output_dir, paste0(file_base, " full_avg.png")), 
           p_full_avg, width = 15, height = 10, dpi = 300)
    ggsave(file.path(output_dir, paste0(file_base, " cond_avg.png")), 
           p_cond_avg, width = 15, height = 10, dpi = 300)
    ggsave(file.path(output_dir, paste0(file_base, " cicbar.png")), 
           p_cicbar, width = 10, height = 8, dpi = 300)
    ggsave(file.path(output_dir, paste0(file_base, " summary.png")), 
           p_summary, width = 12, height = 10, dpi = 300)
    
    # Combined DAG figure - increased height
    ggsave(file.path(output_dir, paste0(file_base, " dags_combined.png")), 
           dag_combined, width = 34, height = 10, dpi = 300)
    
    cat("PNG files saved to:", output_dir, "\n")
  }
  
  if (save_pdf) {
    pdf(file.path(output_dir, paste0(file_base, " all_plots.pdf")), 
        width = 15, height = 10)
    print(p_best)
    print(p_full_avg)
    print(p_cond_avg)
    print(p_cicbar)
    print(p_summary)
    dev.off()
    cat("PDF saved to:", file.path(output_dir, paste0(file_base, " all_plots.pdf")), "\n")
  }
  
  return(list(
    best_model = p_best,
    full_avg = p_full_avg,
    cond_avg = p_cond_avg,
    cicbar = p_cicbar,
    summary = p_summary,
    dags_combined = dag_combined
  ))
}

#' Create all plots for downsampled phylopath runs
#' Formerly create_downsampled_plots() - replacing with below function
#'
#' @param downsampling_results Output from run_multiple_phylopath()
#' @param downsampling_info Information about the downsampling
#' @param output_prefix Prefix including all variables and Remove[n][group] format
#' @param save_png Save plots as PNG
#' @param save_pdf Save plots as PDF
create_downsampled_plots_inflexible <- function(downsampling_results,
                                     detailed_models_input = NULL,
                                     downsampling_info = NULL,
                                     output_prefix = "phylopath_downsampled",
                                     output_dir = "Outputs/PhylopathPlots",
                                     save_png = TRUE,
                                     save_pdf = FALSE) {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 1. Species counts and fractions plot
  if (!is.null(downsampling_info)) {
    cat("Creating downsampling information plot...\n")
    
    # Prepare data for plotting
    if ("proportions" %in% names(downsampling_info)) {
      prop_data <- downsampling_info$proportions
      
      # Handle different column names (territoriality vs group)
      if ("territoriality" %in% colnames(prop_data)) {
        prop_data$group <- factor(prop_data$territoriality, levels = c("Low", "High"))
      } else if ("group" %in% colnames(prop_data)) {
        # Already has group column, just ensure it's a factor
        prop_data$group <- factor(prop_data$group)
      }
      
      # Create stacked bar plot
      prop_data_long <- prop_data %>%
        mutate(
          without_data = total_species - species_with_data
        ) %>%
        pivot_longer(cols = c(species_with_data, without_data),
                     names_to = "data_status",
                     values_to = "count")
      
      p_downsample_info <- ggplot(prop_data_long, 
                                  aes(x = group, y = count, fill = data_status)) +
        geom_bar(stat = "identity", position = "stack", alpha = 0.8) +
        geom_text(data = prop_data,
                  aes(x = group, y = total_species + 50, 
                      label = paste0("n = ", total_species, "\n",
                                     round(proportion * 100, 1), "% with data")),
                  inherit.aes = FALSE, size = 4) +
        scale_fill_manual(values = c("species_with_data" = "#2E86AB", 
                                     "without_data" = "#F4A261"),  # Changed from gray
                          labels = c("species_with_data" = "With FS & CB data",
                                     "without_data" = "Missing data"),
                          name = "") +
        labs(
          title = "Data Availability by Territoriality",
          subtitle = paste("Target proportion:", round(downsampling_info$target_proportion, 3),
                           "\nRemove", downsampling_info$n_to_remove, "species from",
                           downsampling_info$downsample_info$group_to_downsample),
          x = "Territoriality",
          y = "Number of Species"
        ) +
        theme_cowplot(12) +
        theme(legend.position = "bottom")
    } else {
      p_downsample_info <- NULL
    }
  } else {
    p_downsample_info <- NULL
  }
  
  # 2. Process detailed models if available
  # Fix: Look for files with pattern that matches the output_prefix
  detailed_models_files <- list.files(pattern = paste0("detailed_models_.*", 
                                                       gsub(" ", "_", output_prefix), 
                                                       ".*\\.csv$"))
  
  # If no files found, try a more general pattern
  if (length(detailed_models_files) == 0) {
    # Extract the key part of the prefix (e.g., "Remove66HighTerr")
    prefix_parts <- strsplit(output_prefix, " ")[[1]]
    key_pattern <- grep("Remove", prefix_parts, value = TRUE)
    if (length(key_pattern) > 0) {
      detailed_models_files <- list.files(pattern = paste0("detailed_models_.*", key_pattern, ".*\\.csv$"))
    }
  }
  
  detailed_models_file <- if (length(detailed_models_files) > 0) detailed_models_files[1] else NA
  
  if (!is.na(detailed_models_file)) {
    cat("Loading detailed models from:", detailed_models_file, "\n")
    detailed_models_df <- read.csv(detailed_models_file)
    
    # Extract edge information
    edge_cols <- grep("_est$", names(detailed_models_df), value = TRUE)
    
    # Aggregate by seed using the aggregate_by_seed function
    agg_results <- aggregate_by_seed(detailed_models_df, cutoff = 2)
    
    # Convert to dataframe format for plotting
    seed_dataframes <- convert_all_seed_results_to_dataframes(agg_results)
    seed_level_conditional <- seed_dataframes$conditionalAverage_coefficient_perSeed
    
    # 3. Violin plots of path coefficients
    cat("Creating violin plots of path coefficients...\n")
    
    # Clean path names for display
    seed_level_conditional$path_clean <- gsub("_to_", " → ", seed_level_conditional$path)
    seed_level_conditional$path_clean <- gsub("FemaleSong_Agg01", "Female Song", seed_level_conditional$path_clean)
    seed_level_conditional$path_clean <- gsub("HighConfidence_Coop", "Cooperation", seed_level_conditional$path_clean)
    seed_level_conditional$path_clean <- gsub("TerritorialityWeakVsStrong", "Territoriality", seed_level_conditional$path_clean)
    seed_level_conditional$path_clean <- gsub("Territory_12vs3", "Territory Type", seed_level_conditional$path_clean)
    seed_level_conditional$path_clean <- gsub("logMass_AVONET", "Body Mass", seed_level_conditional$path_clean)
    
    # Filter to paths that appear in multiple seeds
    path_counts <- seed_level_conditional %>%
      group_by(path) %>%
      summarise(n = n()) %>%
      filter(n > 10)
    
    seed_level_plot <- seed_level_conditional %>%
      filter(path %in% path_counts$path)
    
    p_violin <- ggplot(seed_level_plot, 
                       aes(x = path_clean, y = mean_coefficient)) +
      geom_violin(fill = "lightblue", alpha = 0.7) +
      geom_boxplot(width = 0.1, alpha = 0.8, outlier.shape = NA) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      coord_flip() +
      labs(
        title = "Distribution of Path Coefficients Across Iterations",
        subtitle = paste("Based on", length(unique(seed_level_plot$seed)), 
                         "downsampling iterations (conditional averaging)"),
        x = NULL,
        y = "Mean Coefficient"
      ) +
      theme_cowplot(12) +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 10)
      )
    
    # 4. Heatmap of mean coefficients
    cat("Creating heatmap of mean coefficients...\n")
    
    # Get the across-seeds summary
    path_summary <- agg_results$across_seeds_coefficient_summary$conditionalAverage_coefficient_acrossSeeds
    
    # Create from-to matrix for heatmap
    if (nrow(path_summary) > 0) {
      # Create a matrix of coefficients
      unique_from <- unique(path_summary$from)
      unique_to <- unique(path_summary$to)
      
      coef_matrix <- matrix(NA, 
                            nrow = length(unique_from), 
                            ncol = length(unique_to),
                            dimnames = list(unique_from, unique_to))
      
      for (i in 1:nrow(path_summary)) {
        coef_matrix[path_summary$from[i], path_summary$to[i]] <- path_summary$mean_coef_conditional[i]
      }
      
      # Convert to long format for ggplot
      coef_long <- as.data.frame(as.table(coef_matrix))
      names(coef_long) <- c("From", "To", "Coefficient")
      coef_long <- coef_long[!is.na(coef_long$Coefficient), ]
      
      # Clean variable names
      clean_names <- function(x) {
        x <- gsub("FemaleSong_Agg01", "Female\nSong", x)
        x <- gsub("HighConfidence_Coop", "Cooperation", x)
        x <- gsub("TerritorialityWeakVsStrong", "Territoriality", x)
        x <- gsub("Territory_12vs3", "Territory\nType", x)
        x <- gsub("logMass_AVONET", "Body Mass", x)
        return(x)
      }
      
      coef_long$From <- clean_names(coef_long$From)
      coef_long$To <- clean_names(coef_long$To)
      
      # Add significance info if available
      sig_info <- path_summary %>%
        mutate(
          From_clean = clean_names(from),
          To_clean = clean_names(to),
          sig_label = ifelse(mean_prop_sig_conditional > 0.95, "***",
                             ifelse(mean_prop_sig_conditional > 0.8, "**",
                                    ifelse(mean_prop_sig_conditional > 0.5, "*", "")))
        )
      
      coef_long <- coef_long %>%
        left_join(sig_info, by = c("From" = "From_clean", "To" = "To_clean"))
      
      p_heatmap <- ggplot(coef_long, aes(x = From, y = To, fill = Coefficient)) +
        geom_tile(color = "white", size = 0.5) +
        geom_text(aes(label = paste0(round(Coefficient, 3), sig_label)),
                  color = ifelse(abs(coef_long$Coefficient) > 0.3, "white", "black"),
                  size = 4) +
        scale_fill_gradient2(low = "#E63946", mid = "white", high = "#2E86AB",
                             midpoint = 0, name = "Mean\nCoefficient",
                             limits = c(-max(abs(coef_long$Coefficient)), 
                                        max(abs(coef_long$Coefficient)))) +
        labs(
          title = "Mean Path Coefficients Across Iterations",
          subtitle = "Significance: *** >95%, ** >80%, * >50% of iterations",
          x = "From Variable",
          y = "To Variable"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1)
        )
    } else {
      p_heatmap <- NULL
    }
    
    # 5. Model frequency bar plot - KEEP ORIGINAL MODEL NAMES
    cat("Creating model frequency plot...\n")
    if (!is.null(downsampling_results$model_frequencies)) {
      model_freq <- downsampling_results$model_frequencies %>%
        arrange(desc(frequency_in_sub2)) %>%
        head(15) %>%
        mutate(
          prop_sub2 = frequency_in_sub2 / max(downsampling_results$results_df$seed)
        )
      
      p_model_freq <- ggplot(model_freq, 
                             aes(x = reorder(model_name, frequency_in_sub2),  # Use original model names
                                 y = frequency_in_sub2)) +
        geom_bar(stat = "identity", fill = "#2E86AB", alpha = 0.8) +
        geom_text(aes(label = paste0(frequency_in_sub2, " (", 
                                     round(prop_sub2 * 100, 1), "%)")),
                  hjust = -0.1, size = 3) +
        coord_flip() +
        labs(
          title = "Top Models Across Downsampling Iterations",
          subtitle = paste("Models appearing in Δ CICc < 2 set across",
                           length(unique(downsampling_results$results_df$seed)), 
                           "iterations"),
          x = NULL,
          y = "Frequency in Top Model Set"
        ) +
        theme_cowplot(12) +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          axis.text.y = element_text(size = 9)
        ) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
    } else {
      p_model_freq <- NULL
    }
    
  } else {
    cat("No detailed models file found\n")
    p_violin <- NULL
    p_heatmap <- NULL
    p_model_freq <- NULL
  }
  
  # Save all plots with proper file names
  if (save_png) {
    if (!is.null(p_downsample_info)) {
      ggsave(file.path(output_dir, paste0(output_prefix, " downsample_info.png")),
             p_downsample_info, width = 8, height = 6, dpi = 300)
    }
    if (!is.null(p_violin)) {
      ggsave(file.path(output_dir, paste0(output_prefix, " violin.png")),
             p_violin, width = 10, height = 8, dpi = 300)
    }
    if (!is.null(p_heatmap)) {
      ggsave(file.path(output_dir, paste0(output_prefix, " heatmap.png")),
             p_heatmap, width = 8, height = 8, dpi = 300)
    }
    if (!is.null(p_model_freq)) {
      ggsave(file.path(output_dir, paste0(output_prefix, " model_freq.png")),
             p_model_freq, width = 10, height = 8, dpi = 300)
    }
    
    cat("PNG files saved to:", output_dir, "\n")
  }
  
  return(list(
    downsample_info = p_downsample_info,
    violin = p_violin,
    heatmap = p_heatmap,
    model_freq = p_model_freq
  ))
}


#' Updated create_downsampled_plots with flexible input
#' Formerly create_downsampled_plots_flexible()
#' Further developed in create_downsampled_plots_forDev2025-06-11.R to plot conditional average coefficients from full dataset on violin-box plots
#'
#' This version accepts detailed_models_df as:
#' A) A dataframe object
#' B) A file path to the detailed_models csv
#' C) NULL (searches for csv file using pattern matching)
#'
#' @param downsampling_results Output from run_multiple_phylopath or similar
#' @param detailed_models_input Either a dataframe, file path, or NULL
#' @param model_frequencies_input Either a dataframe, file path, or NULL
#' @param full_data_phylopath_input Either a phylopath result object or NULL
#' @param full_dataset Full dataset for running phylopath (if full_data_phylopath_input is NULL)
#' @param tree Phylogenetic tree for running phylopath (if full_data_phylopath_input is NULL)
#' @param downsampling_info Information about the downsampling
#' @param output_prefix Prefix for output files
#' @param output_dir Output directory
#' @param save_png Save plots as PNG
#' @param save_pdf Save plots as PDF
create_downsampled_plots <- function(downsampling_results,
                                     detailed_models_input = NULL,
                                     model_frequencies_input = NULL,
                                     full_data_phylopath_input = NULL,
                                     full_dataset = NULL,
                                     tree = NULL,
                                     downsampling_info = NULL,
                                     output_prefix = "phylopath_downsampled",
                                     output_dir = "Outputs/PhylopathPlots",
                                     save_png = TRUE,
                                     save_pdf = TRUE) {
  
  # Load required libraries
  require(ggplot2)
  require(dplyr)
  require(cowplot)
  require(phylopath)
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Skip downsampling info plot (not useful, prone to breaking)
  p_downsample_info <- NULL
  
  # Handle detailed models input
  detailed_models_df <- NULL
  
  if (is.data.frame(detailed_models_input)) {
    # Case A: Already a dataframe
    cat("Using provided detailed models dataframe\n")
    detailed_models_df <- detailed_models_input
    
  } else if (is.character(detailed_models_input) && file.exists(detailed_models_input)) {
    # Case B: File path provided
    cat("Loading detailed models from:", detailed_models_input, "\n")
    detailed_models_df <- read.csv(detailed_models_input)
    
  } else if (is.null(detailed_models_input)) {
    # Case C: Search for file using pattern
    cat("Searching for detailed models file...\n")
    
    # First try the standard pattern
    detailed_models_files <- list.files(pattern = paste0("detailed_models_.*", 
                                                         gsub(" ", "_", output_prefix), 
                                                         ".*\\.csv$"))
    
    # If no files found, try a more general pattern
    if (length(detailed_models_files) == 0) {
      # Extract the key part of the prefix (e.g., "Remove66HighTerr")
      prefix_parts <- strsplit(output_prefix, " ")[[1]]
      key_pattern <- grep("Remove|High", prefix_parts, value = TRUE)
      if (length(key_pattern) > 0) {
        pattern_str <- paste(key_pattern, collapse = ".*")
        detailed_models_files <- list.files(pattern = paste0("detailed_models_.*", pattern_str, ".*\\.csv$"))
      }
    }
    
    # If still no files, try searching in output directory
    if (length(detailed_models_files) == 0 && output_dir != ".") {
      detailed_models_files <- list.files(output_dir, 
                                          pattern = paste0("detailed_models_.*", 
                                                           gsub(" ", "_", output_prefix), 
                                                           ".*\\.csv$"),
                                          full.names = TRUE)
    }
    
    if (length(detailed_models_files) > 0) {
      detailed_models_file <- detailed_models_files[1]
      cat("Loading detailed models from:", detailed_models_file, "\n")
      detailed_models_df <- read.csv(detailed_models_file)
    } else {
      cat("No detailed models file found\n")
    }
    
  } else if (!is.null(downsampling_results$detailed_models)) {
    # Check if detailed_models is included in the results object
    cat("Using detailed models from results object\n")
    detailed_models_df <- downsampling_results$detailed_models
  }
  
  # Handle model frequencies input
  model_frequencies_df <- NULL
  
  if (is.data.frame(model_frequencies_input)) {
    # Case A: Already a dataframe
    cat("Using provided model frequencies dataframe\n")
    model_frequencies_df <- model_frequencies_input
    
  } else if (is.character(model_frequencies_input) && file.exists(model_frequencies_input)) {
    # Case B: File path provided
    cat("Loading model frequencies from:", model_frequencies_input, "\n")
    model_frequencies_df <- read.csv(model_frequencies_input)
    
  } else if (is.null(model_frequencies_input)) {
    # Case C: Try to get from downsampling_results or search for file
    if (!is.null(downsampling_results$model_frequencies)) {
      cat("Using model frequencies from results object\n")
      model_frequencies_df <- downsampling_results$model_frequencies
    } else {
      # Search for file using pattern
      cat("Searching for model frequencies file...\n")
      
      model_freq_files <- list.files(pattern = paste0("model_frequencies_.*", 
                                                      gsub(" ", "_", output_prefix), 
                                                      ".*\\.csv$"))
      
      if (length(model_freq_files) == 0 && output_dir != ".") {
        model_freq_files <- list.files(output_dir, 
                                       pattern = paste0("model_frequencies_.*", 
                                                        gsub(" ", "_", output_prefix), 
                                                        ".*\\.csv$"),
                                       full.names = TRUE)
      }
      
      if (length(model_freq_files) > 0) {
        model_freq_file <- model_freq_files[1]
        cat("Loading model frequencies from:", model_freq_file, "\n")
        model_frequencies_df <- read.csv(model_freq_file)
      }
    }
  }
  
  # Process detailed models if available
  if (!is.null(detailed_models_df) && nrow(detailed_models_df) > 0) {
    # Extract edge information
    edge_cols <- grep("_est$", names(detailed_models_df), value = TRUE)
    
    # Aggregate by seed
    agg_results <- aggregate_by_seed(detailed_models_df, cutoff = 2)
    
    # Convert to dataframe format for plotting
    seed_dataframes <- convert_all_seed_results_to_dataframes(agg_results)
    seed_level_conditional <- seed_dataframes$conditionalAverage_coefficient_perSeed
    
    # Process full dataset phylopath results if provided
    full_data_coefficients <- NULL
    
    if (!is.null(full_data_phylopath_input)) {
      cat("Processing full dataset phylopath results...\n")
      
      if (is.list(full_data_phylopath_input) && "result" %in% names(full_data_phylopath_input)) {
        # Extract conditional average coefficients from the phylopath result
        full_result <- full_data_phylopath_input$result
        
        # Get best models (delta_CICc < 2)
        full_summary <- summary(full_result)
        best_models <- full_summary[full_summary$delta_CICc < 2, ]
        
        if (nrow(best_models) > 0) {
          # Extract coefficients for each best model
          all_paths <- list()
          
          for (i in 1:nrow(best_models)) {
            model_name <- best_models$model[i]
            chosen_model <- choice(full_result, model_name)
            
            if (!is.null(chosen_model$coef) && is.matrix(chosen_model$coef)) {
              coef_matrix <- chosen_model$coef
              
              # Convert to path format
              for (from_idx in 1:nrow(coef_matrix)) {
                for (to_idx in 1:ncol(coef_matrix)) {
                  coef_value <- coef_matrix[from_idx, to_idx]
                  if (coef_value != 0) {
                    from_name <- rownames(coef_matrix)[from_idx]
                    to_name <- colnames(coef_matrix)[to_idx]
                    path_name <- paste0(from_name, "_to_", to_name)
                    
                    all_paths[[path_name]] <- c(all_paths[[path_name]], coef_value)
                  }
                }
              }
            }
          }
          
          # Calculate conditional averages
          full_data_coefficients <- data.frame(
            path = names(all_paths),
            full_data_mean = sapply(all_paths, mean),
            stringsAsFactors = FALSE
          )
        }
      }
    } else {
      # Option B: Run phylopath on full dataset if dataset and tree are provided
      if (!is.null(full_dataset) && !is.null(tree)) {
        cat("Running phylopath on full dataset...\n")
        
        # Load dataset and tree if they are file paths
        if (is.character(full_dataset)) {
          cat("Loading dataset from:", full_dataset, "\n")
          full_dataset <- read.csv(full_dataset)
        }
        
        if (is.character(tree)) {
          cat("Loading tree from:", tree, "\n")
          require(ape)
          if (grepl("\\.nex", tree)) {
            tree <- read.nexus(tree)
          } else if (grepl("\\.nwk|\\.tre", tree)) {
            tree <- read.tree(tree)
          }
        }
        
        # Filter dataset to species in tree
        full_dataset <- full_dataset[full_dataset$species %in% tree$tip.label, ]
        rownames(full_dataset) <- full_dataset$species
        
        # Determine which variables are being used based on output_prefix
        # Look for common patterns in the prefix
        if (grepl("PlumageDimorphism", output_prefix)) {
          mass_var <- "logMaleFemalePlumageDiffAbs"
        } else if (grepl("WingDimorphism", output_prefix)) {
          mass_var <- "PercentAbsLogWingDimorphism"
        } else {
          mass_var <- "logMass_AVONET"  # Default
        }
        
        # Determine territoriality variable
        if (grepl("Territory_12vs3|Terr3", output_prefix)) {
          territoriality_var <- "Territory_12vs3"
        } else {
          territoriality_var <- "TerritorialityWeakVsStrong"  # Default
        }
        
        # Run phylopath
        full_result <- run_CB_FS_Terr_phylopath(
          dfIn = full_dataset,
          tree = tree,
          female_song_var = "FemaleSong_Agg01",
          coop_breeding_var = "HighConfidence_Coop",
          territoriality_var = territoriality_var,
          mass_var = mass_var,
          plots2pdf = FALSE
        )
        
        # Now process the result
        if (!is.null(full_result) && !is.null(full_result$result)) {
          full_summary <- summary(full_result$result)
          best_models <- full_summary[full_summary$delta_CICc < 2, ]
          
          if (nrow(best_models) > 0) {
            all_paths <- list()
            
            for (i in 1:nrow(best_models)) {
              model_name <- best_models$model[i]
              chosen_model <- choice(full_result$result, model_name)
              
              if (!is.null(chosen_model$coef) && is.matrix(chosen_model$coef)) {
                coef_matrix <- chosen_model$coef
                
                for (from_idx in 1:nrow(coef_matrix)) {
                  for (to_idx in 1:ncol(coef_matrix)) {
                    coef_value <- coef_matrix[from_idx, to_idx]
                    if (coef_value != 0) {
                      from_name <- rownames(coef_matrix)[from_idx]
                      to_name <- colnames(coef_matrix)[to_idx]
                      path_name <- paste0(from_name, "_to_", to_name)
                      
                      all_paths[[path_name]] <- c(all_paths[[path_name]], coef_value)
                    }
                  }
                }
              }
            }
            
            # Calculate conditional averages
            full_data_coefficients <- data.frame(
              path = names(all_paths),
              full_data_mean = sapply(all_paths, mean),
              stringsAsFactors = FALSE
            )
          }
        }
      } else {
        cat("full_data_phylopath_input is NULL and dataset/tree not provided\n")
        cat("To enable full dataset comparison, provide either:\n")
        cat("  1. full_data_phylopath_input = result_bodymass (or similar)\n")
        cat("  2. full_dataset = dfIn_phylo and tree = tree_phylo\n")
      }
    }
    
    # Create VIOLIN-BOX PLOTS
    cat("Creating violin-box plots of path coefficients...\n")
    
    # Clean path names for display
    seed_level_conditional$path_clean <- gsub("_to_", " → ", seed_level_conditional$path)
    seed_level_conditional$path_clean <- gsub("FemaleSong_Agg01", "Female Song", seed_level_conditional$path_clean)
    seed_level_conditional$path_clean <- gsub("HighConfidence_Coop", "Cooperation", seed_level_conditional$path_clean)
    seed_level_conditional$path_clean <- gsub("TerritorialityWeakVsStrong", "Strong Territoriality", seed_level_conditional$path_clean)
    seed_level_conditional$path_clean <- gsub("Territory_12vs3", "Year-round Territory", seed_level_conditional$path_clean)
    seed_level_conditional$path_clean <- gsub("logMass_AVONET", "Body Mass", seed_level_conditional$path_clean)
    seed_level_conditional$path_clean <- gsub("logMaleFemalePlumageDiffAbs", "Plumage Dimorphism", seed_level_conditional$path_clean)
    seed_level_conditional$path_clean <- gsub("PercentAbsLogWingDimorphism", "Wing Dimorphism", seed_level_conditional$path_clean)
    
    # Filter to paths that appear in multiple seeds
    path_counts <- seed_level_conditional %>%
      group_by(path) %>%
      summarise(n = n()) %>%
      filter(n > 1)
    
    seed_level_plot <- seed_level_conditional %>%
      filter(path %in% path_counts$path)
    
    # Prepare full dataset coefficients for plotting if available
    full_data_plot <- NULL
    if (!is.null(full_data_coefficients)) {
      # Clean path names to match
      full_data_coefficients$path_clean <- gsub("_to_", " → ", full_data_coefficients$path)
      full_data_coefficients$path_clean <- gsub("FemaleSong_Agg01", "Female Song", full_data_coefficients$path_clean)
      full_data_coefficients$path_clean <- gsub("HighConfidence_Coop", "Cooperation", full_data_coefficients$path_clean)
      full_data_coefficients$path_clean <- gsub("TerritorialityWeakVsStrong", "Strong Territoriality", full_data_coefficients$path_clean)
      full_data_coefficients$path_clean <- gsub("Territory_12vs3", "Year-round Territory", full_data_coefficients$path_clean)
      full_data_coefficients$path_clean <- gsub("logMass_AVONET", "Body Mass", full_data_coefficients$path_clean)
      full_data_coefficients$path_clean <- gsub("logMaleFemalePlumageDiffAbs", "Plumage Dimorphism", full_data_coefficients$path_clean)
      full_data_coefficients$path_clean <- gsub("PercentAbsLogWingDimorphism", "Wing Dimorphism", full_data_coefficients$path_clean)
      
      # Filter to paths that are in the downsampled data
      full_data_plot <- full_data_coefficients[full_data_coefficients$path_clean %in% seed_level_plot$path_clean, ]
    }
    
    p_violin <- ggplot(seed_level_plot, 
                       aes(x = path_clean, y = mean_coefficient)) +
      geom_violin(fill = "lightblue", alpha = 0.7) +
      geom_boxplot(width = 0.1, alpha = 0.8, outlier.shape = NA) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      coord_flip() +
      labs(
        title = "Distribution of Path Coefficients Across Iterations",
        subtitle = paste("Based on", length(unique(seed_level_plot$seed)), 
                         "downsampling iterations (conditional averaging)"),
        x = NULL,
        y = "Mean Coefficient"
      ) +
      theme_cowplot(12) +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 10)
      )
    
    # Add full dataset coefficients if available
    if (!is.null(full_data_plot) && nrow(full_data_plot) > 0) {
      # Use geom_crossbar for a clean horizontal line that appears vertical after coord_flip
      p_violin <- p_violin +
        geom_crossbar(data = full_data_plot,
                      aes(x = path_clean, y = full_data_mean, 
                          ymin = full_data_mean, ymax = full_data_mean),
                      width = 0.4, color = "red", size = 0.5)
      
      # Update subtitle to indicate full dataset overlay
      p_violin <- p_violin +
        labs(subtitle = paste("Based on", length(unique(seed_level_plot$seed)), 
                              "downsampling iterations (conditional averaging)",
                              "\nRed lines show full dataset conditional averages"))
    }
    
    # Create COEFFICIENT HEATMAP
    cat("Creating heatmap of mean coefficients...\n")
    
    # Get the across-seeds summary
    path_summary <- agg_results$across_seeds_coefficient_summary$conditionalAverage_coefficient_acrossSeeds
    
    # Create from-to matrix for heatmap
    if (nrow(path_summary) > 0) {
      # Create a matrix of coefficients
      unique_from <- unique(path_summary$from)
      unique_to <- unique(path_summary$to)
      
      coef_matrix <- matrix(NA, 
                            nrow = length(unique_from), 
                            ncol = length(unique_to),
                            dimnames = list(unique_from, unique_to))
      
      for (i in 1:nrow(path_summary)) {
        coef_matrix[path_summary$from[i], path_summary$to[i]] <- path_summary$mean_coef_conditional[i]
      }
      
      # Convert to long format for ggplot
      coef_long <- as.data.frame(as.table(coef_matrix))
      names(coef_long) <- c("From", "To", "Coefficient")
      coef_long <- coef_long[!is.na(coef_long$Coefficient), ]
      
      # Clean variable names
      clean_names <- function(x) {
        x <- gsub("FemaleSong_Agg01", "Female\nSong", x)
        x <- gsub("HighConfidence_Coop", "Cooperation", x)
        x <- gsub("TerritorialityWeakVsStrong", "Strong\nTerritoriality", x)
        x <- gsub("Territory_12vs3", "Year-Round\nTerritory\nType", x)
        x <- gsub("logMass_AVONET", "Body Mass", x)
        return(x)
      }
      
      coef_long$From <- clean_names(coef_long$From)
      coef_long$To <- clean_names(coef_long$To)
      
      # Add significance info if available
      sig_info <- path_summary %>%
        mutate(
          From_clean = clean_names(from),
          To_clean = clean_names(to),
          sig_label = ifelse(mean_prop_sig_conditional > 0.95, "***",
                             ifelse(mean_prop_sig_conditional > 0.8, "**",
                                    ifelse(mean_prop_sig_conditional > 0.5, "*", "")))
        )
      
      coef_long <- coef_long %>%
        left_join(sig_info, by = c("From" = "From_clean", "To" = "To_clean"))
      
      p_heatmap <- ggplot(coef_long, aes(x = From, y = To, fill = Coefficient)) +
        geom_tile(color = "white", size = 0.5) +
        geom_text(aes(label = paste0(round(Coefficient, 3), sig_label)),
                  color = ifelse(abs(coef_long$Coefficient) > 0.3, "white", "black"),
                  size = 4) +
        scale_fill_gradient2(low = "#E63946", mid = "white", high = "#2E86AB",
                             midpoint = 0, name = "Mean\nCoefficient",
                             limits = c(-max(abs(coef_long$Coefficient)), 
                                        max(abs(coef_long$Coefficient)))) +
        labs(
          title = "Mean Path Coefficients Across Iterations",
          subtitle = "Significance: *** >95%, ** >80%, * >50% of iterations",
          x = "From Variable",
          y = "To Variable"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1)
        )
    } else {
      p_heatmap <- NULL
    }
    
    # MODEL FREQUENCY BAR PLOT
    cat("Creating model frequency plot...\n")
    if (!is.null(model_frequencies_df)) {
      model_freq <- model_frequencies_df %>%
        arrange(desc(frequency_in_sub2)) %>%
        head(15) %>%
        mutate(
          prop_sub2 = frequency_in_sub2 / length(unique(detailed_models_df$seed))
        )
      
      p_model_freq <- ggplot(model_freq, 
                             aes(x = reorder(model_name, frequency_in_sub2),
                                 y = frequency_in_sub2)) +
        geom_bar(stat = "identity", fill = "#2E86AB", alpha = 0.8) +
        geom_text(aes(label = paste0(frequency_in_sub2, " (", 
                                     round(prop_sub2 * 100, 1), "%)")),
                  hjust = -0.1, size = 3) +
        coord_flip() +
        labs(
          title = "Top Models Across Downsampling Iterations",
          subtitle = paste("Models appearing in Δ CICc < 2 set across",
                           length(unique(detailed_models_df$seed)), 
                           "iterations"),
          x = NULL,
          y = "Frequency in Top Model Set"
        ) +
        theme_cowplot(12) +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          axis.text.y = element_text(size = 9)
        ) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
    } else {
      p_model_freq <- NULL
    }
    
  } else {
    cat("No detailed models data available for plotting\n")
    p_violin <- NULL
    p_heatmap <- NULL
    p_model_freq <- NULL
  }
  
  # Save all plots with proper file names
  if (save_png) {
    if (!is.null(p_downsample_info)) {
      ggsave(file.path(output_dir, paste0(output_prefix, " downsample_info.png")),
             p_downsample_info, width = 8, height = 6, dpi = 300)
    }
    if (!is.null(p_violin)) {
      ggsave(file.path(output_dir, paste0(output_prefix, " violin-box_coefficients.png")),
             p_violin, width = 10, height = 8, dpi = 300)
    }
    if (!is.null(p_heatmap)) {
      ggsave(file.path(output_dir, paste0(output_prefix, " heatmap_coefficients.png")),
             p_heatmap, width = 8, height = 8, dpi = 300)
    }
    if (!is.null(p_model_freq)) {
      ggsave(file.path(output_dir, paste0(output_prefix, " model_freq.png")),
             p_model_freq, width = 10, height = 8, dpi = 300)
    }
    
    cat("PNG files saved to:", output_dir, "\n")
  }
  
  return(list(
    downsample_info = p_downsample_info,
    violin = p_violin,
    heatmap = p_heatmap,
    model_freq = p_model_freq
  ))
}


#' Wrapper function to create all phylopath plots
#'
#' @param analysis_type Either "nondownsampled" or "downsampled"
#' @param phylopath_output Output from phylopath analysis
#' @param downsampling_info Optional downsampling information
#' @param output_prefix Prefix for output files
#' @param output_dir Output directory
#' @param detailed_models_input Either a dataframe, file path, or NULL
#' @param model_frequencies_input Either a dataframe, file path, or NULL
#' @param full_data_phylopath_input Either a phylopath result object from non-downsampled analysis or NULL
#' @param full_dataset Full dataset for running phylopath (if full_data_phylopath_input is NULL and analysis_type is "downsampled")
#' @param tree Phylogenetic tree for running phylopath (if full_data_phylopath_input is NULL and analysis_type is "downsampled")
#' @param save_png Save as PNG
#' @param save_pdf Save as PDF
create_all_phylopath_plots <- function(analysis_type = c("nondownsampled", "downsampled"),
                                       phylopath_output,
                                       detailed_models_input = NULL,
                                       downsampling_info = NULL,
                                       output_prefix = "phylopath",
                                       output_dir = "Outputs/PhylopathPlots",
                                       model_frequencies_input = NULL,
                                       full_data_phylopath_input = NULL,
                                       full_dataset = NULL,
                                       tree = NULL,
                                       save_png = TRUE,
                                       save_pdf = FALSE) {
  
  require(cowplot)
  
  analysis_type <- match.arg(analysis_type)
  
  if (analysis_type == "nondownsampled") {
    plots <- create_nondownsampled_plots(
      phylopath_output = phylopath_output,
      output_prefix = output_prefix,
      output_dir = output_dir,
      save_png = save_png,
      save_pdf = save_pdf
    )
  } else {
    plots <- create_downsampled_plots(
      downsampling_results = phylopath_output,
      detailed_models_input = detailed_models_input,  
      downsampling_info = downsampling_info,
      output_prefix = output_prefix,
      output_dir = output_dir,
      model_frequencies_input = model_frequencies_input,
      full_data_phylopath_input = full_data_phylopath_input,
      full_dataset = full_dataset,
      tree = tree,
      save_png = save_png,
      save_pdf = save_pdf
    )
  }
  
  return(plots)
}


