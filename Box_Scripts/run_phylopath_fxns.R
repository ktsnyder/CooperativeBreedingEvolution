### Repeatedly downsample species to account for biases in data
# Kate Snyder
# May 15, 2025

library(phylopath)
library(phytools)


#### Functions ----

# Function to run multiple iterations and aggregate results
run_multiple_phylopath <- function(dfIn, tree, downsample_columns, downsample_values, numToRemove,
                                   n_iterations = 500, 
                                   female_song_var = "FemaleSong_Agg01", 
                                   coop_breeding_var = "HighConfidence_Coop", 
                                   territoriality_var = "TerritorialityWeakVsStrong", 
                                   mass_var = "logMass_AVONET",
                                   save_conditional_plots = TRUE,
                                   save_path_coefficients = TRUE,
                                   plotlabel = "") {
  
  require(dplyr)
  require(ggplot2)
  require(gridExtra)
  
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
    write.csv(detailed_models_df, paste0("detailed_models_", plotlabel, "_", n_iterations, "_", Sys.Date(), ".csv"), row.names = FALSE)
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
    pdf_name <- paste0("conditional_average_plots_", plotlabel, "_", n_iterations, "_", Sys.Date(), ".pdf")
    
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
  write.csv(model_freq_df, paste0("model_frequencies_", plotlabel, "_", n_iterations, "_", Sys.Date(), ".csv"), row.names = FALSE)
  
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
    pdf(paste0("top_models_", plotlabel, "_", n_iterations, "_", Sys.Date(), ".pdf"), width = 10, height = 8)
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
    write.csv(edge_summary, paste0("edge_summary_", plotlabel, "_", n_iterations, "_", Sys.Date(), ".csv"), row.names = FALSE)
    
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

run_CB_FS_Terr_phylopath <- function(dfIn, tree, female_song_var, coop_breeding_var, territoriality_var, mass_var = "none", include_terr_response = FALSE, plots2pdf = FALSE) {
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
  full_average_plot <- plot(average(result, cut_off = 2, avg_method = "full"), text_size = 3, manual_layout = phylopath_map_positions)
  conditional_average_plot <- plot(average(result, cut_off = 2, avg_method = "conditional"), text_size = 3, manual_layout = phylopath_map_positions)
  
  CICsub2_models = s$model[which(s$delta_CICc<2)]
  print(paste("There are", length(CICsub2_models), "models with delta CICc < 2"))
  
  
  topModelPlotList <- list()
  topModelPlotTitles <- list()
  for (i in 1:length(CICsub2_models)) {
    tempmodel = s$model[which(s$delta_CICc<2)][i]
    tempCIC = s$CICc[which(s$delta_CICc<2)][i]
    temp_delta_CIC = s$delta_CICc[which(s$delta_CICc<2)][i]
    tempmodel_clean <- gsub("→", "to", tempmodel)
    tempmodel_title <- paste(tempmodel_clean, "| CICc:", round(tempCIC, 2), "| delta_CICc:", round(temp_delta_CIC, 2))
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
    pdf(pdfname, width = 10, height = 8)
    
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
        theme(plot.title = element_text(size = 6, face = "bold"))
      
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
            theme(plot.title = element_text(size = 12, face = "bold")))
    
    # Print conditional average model with title
    print(conditional_average_plot + 
            ggtitle("Conditional averaged model") + 
            theme(plot.title = element_text(size = 12, face = "bold")))
    
    # Close the PDF device
    dev.off()
  }
  
  
  outlist <- list(
    var_map = var_map,
    result = result,
    summary_plot = s_plot,
    best_model_plot = best_model_plot,
    full_average_plot = full_average_plot,
    conditional_average_plot = conditional_average_plot,
    topModelPlotList = topModelPlotList,
    topModelPlotTitles = topModelPlotTitles,
    CICsub2_models = CICsub2_models,
    nSpecies = nSpecies
  )
  
  return(outlist)
  
} # end function run_CB_FS_Terr_phylopath


## Function to downsample dataset and run run_CB_FS_Terr_phylopath()
downsample_run_phylopath <- function(dfIn, tree, downsample_columns, downsample_values, numToRemove, seed, female_song_var = "FemaleSong_Agg01", coop_breeding_var = "HighConfidence_Coop", territoriality_var = "TerritorialityWeakVsStrong", mass_var = "none") {
  
  dfIn$GeographicRegion_Jetz = NA
  dfIn$GeographicRegion_Jetz[which(dfIn$Realm_Jetz2011 %in% c("PA", "NeA"))] <- "Holarctic"
  dfIn$GeographicRegion_Jetz[which(dfIn$Realm_Jetz2011 %in% c("AT", "NT", "IM", "AA", "OC"))] <- "Tropical"
  
  dfIn$GeographicRegion_Cockburn = NA
  dfIn$GeographicRegion_Cockburn[which(dfIn$Region_Cockburn2006 %in% c("Nearctic", "Holarctic", "Palearctic"))] <- "Holarctic"
  dfIn$GeographicRegion_Cockburn[which(dfIn$Region_Cockburn2006 %in% c("Indomalayan", "Australia", "Neotropical", "Africa"))] <- "Tropical"
  
  # Subset to just Oscine species
  dfOs = dfIn[which(dfIn$species %in% tree$tip.label),]
  dfOs$HaveFSData = !is.na(dfOs$FemaleSong_Agg01)
  dfOs$HaveCBData = !is.na(dfOs$HighConfidence_Coop)
  
  # We want to remove 83 birds that have GeographicRegion_Jetz == "Holarctic", HighConfidence_Coop == 0, HaveFSData == TRUE from the dataset (see "bias tests.R")
  # First, identify the rows that meet all criteria
  if (length(downsample_columns) == 2) {
    matching_rows <- which(dfOs[downsample_columns[1]] == downsample_values[1] & 
                             dfOs[downsample_columns[2]] == downsample_values[2])
  } else {
    matching_rows <- which(dfOs[downsample_columns[1]] == downsample_values[1] & 
                             dfOs[downsample_columns[2]] == downsample_values[2] & 
                             dfOs[downsample_columns[3]] == downsample_values[3])
  }

  # Check how many matching rows we have
  num_matching <- length(matching_rows)
  print(paste("Number of rows matching criteria:", num_matching))
  
  # Make sure there are at least 83 rows to remove
  if (num_matching < numToRemove) {
    stop("There are fewer than numToRemove rows matching the criteria")
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


# Helper function to build formulas from patterns
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


#### Functions for analyzing detailed models output from phylopath ----
# Author: Claude (based on user requirements)

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

# Helper function for null coalescing
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) y else x
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



