# Jackknife species analysis function for phylopath
# This function iteratively removes each species one at a time and runs phylopath analysis
# 7/11/2025 - moved output directory generation to in fxn

library(phylopath)
library(phytools)
library(dplyr)

#' Run jackknife species analysis for phylopath
#'
#' This function iteratively removes each species from the dataset and runs phylopath analysis,
#' tracking which species was removed in each iteration.
#' 
#' @param dfIn Data frame containing all species data
#' @param tree Phylogenetic tree
#' @param female_song_var Variable name for female song
#' @param coop_breeding_var Variable name for cooperative breeding
#' @param territoriality_var Variable name for territoriality
#' @param mass_var Variable name for mass (or "none")
#' @param save_conditional_plots Whether to save conditional average plots
#' @param save_path_coefficients Whether to save path coefficients
#' @param output_dir Output directory for results
#' @return List containing results data frames and summaries
run_jackknife_species_phylopath <- function(dfIn, 
                                          tree,
                                          female_song_var = "FemaleSong_Agg01", 
                                          coop_breeding_var = "HighConfidence_Coop", 
                                          territoriality_var = "TerritorialityWeakVsStrong", 
                                          mass_var = "logMass_AVONET",
                                          save_conditional_plots = FALSE,
                                          save_path_coefficients = TRUE,
                                          output_dir = file.path("Outputs", "PhylopathJackknife")) {
  
  require(dplyr)
  require(ggplot2)
  
  all_traits_phylopath_label <- paste(female_song_var, coop_breeding_var, territoriality_var, mass_var)
  
  output_dir = file.path("Outputs", "PhylopathJackknife", paste(all_traits_phylopath_label, "models"))
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Get list of species that have complete data for all variables
  required_vars <- c(female_song_var, coop_breeding_var, territoriality_var)
  if (mass_var != "none") {
    required_vars <- c(required_vars, mass_var)
  }
  
  # Filter to complete cases and species in tree
  dfIn_complete <- dfIn[complete.cases(dfIn[, required_vars]), ]
  dfIn_complete <- dfIn_complete[dfIn_complete$species %in% tree$tip.label, ]
  
  # Get the list of species to iterate through
  species_list <- dfIn_complete$species
  n_species <- length(species_list)
  
  cat(sprintf("Starting jackknife analysis for %d species...\n", n_species))
  
  # Create a data frame to store results
  results_df <- data.frame(
    seed = integer(n_species),
    RemovedSpecies = character(n_species),
    nSub2dCIC_Models = integer(n_species),
    CICsub2_models_string = character(n_species),
    nSpecies = integer(n_species),
    best_model = character(n_species),
    best_CICc = numeric(n_species),
    stringsAsFactors = FALSE
  )
  
  # Lists to store additional results
  conditional_average_plots <- list()
  path_coefficients <- list()
  detailed_models_list <- list()
  
  # Start time for progress tracking
  start_time <- Sys.time()
  
  # Run iterations
  for (i in 1:n_species) {
    # Use iteration number as seed for reproducibility
    current_seed <- i
    species_to_remove <- species_list[i]
    
    cat(sprintf("Iteration %d/%d: Removing species '%s'...\n", i, n_species, species_to_remove))
    
    # Create dataset without current species
    dfIn_jackknife <- dfIn_complete[dfIn_complete$species != species_to_remove, ]
    
    # Create tree without current species
    tree_jackknife <- drop.tip(tree, species_to_remove)
    
    # Run phylopath analysis
    tryCatch({
      current_output <- run_CB_FS_Terr_phylopath(
        dfIn = dfIn_jackknife,
        tree = tree_jackknife,
        female_song_var = female_song_var,
        coop_breeding_var = coop_breeding_var,
        territoriality_var = territoriality_var,
        mass_var = mass_var,
        plots2pdf = FALSE,
        plots2png = FALSE
      )
      
      # Extract summary from the result
      s <- summary(current_output$result)
      
      # Store basic results
      results_df$seed[i] <- current_seed
      results_df$RemovedSpecies[i] <- species_to_remove
      results_df$nSub2dCIC_Models[i] <- length(current_output$CICsub2_models)
      results_df$CICsub2_models_string[i] <- paste(current_output$CICsub2_models, collapse = ", ")
      results_df$nSpecies[i] <- current_output$nSpecies
      results_df$best_model[i] <- s$model[1]  # First model is the best
      results_df$best_CICc[i] <- s$CICc[1]    # CICc of the best model
      
      # Store the conditional average plot if requested
      if (save_conditional_plots) {
        conditional_average_plots[[i]] <- current_output$conditional_average_plot + 
          ggtitle(paste("Removed:", species_to_remove, "| nSpecies:", current_output$nSpecies)) +
          theme(plot.title = element_text(size = 10, face = "bold"))
      }
      
      # Create detailed model information for each model with delta_CICc < 2
      for (model_idx in which(s$delta_CICc < 2)) {
        model_name <- s$model[model_idx]
        model_CICc <- s$CICc[model_idx]
        model_delta_CICc <- s$delta_CICc[model_idx]
        
        # Get the chosen model to extract path coefficients
        chosen_model <- choice(current_output$result, model_name)
        
        # Create a row for this model
        model_row <- data.frame(
          seed = current_seed,
          RemovedSpecies = species_to_remove,
          model = model_name,
          CICc = model_CICc,
          delta_CICc = model_delta_CICc,
          nSpecies = current_output$nSpecies,
          stringsAsFactors = FALSE
        )
        
        # Extract edge information from coefficient matrix
        edges_coef <- chosen_model$coef
        edges_se <- chosen_model$se
        
        # Process coefficient matrix
        if (!is.null(edges_coef) && is.matrix(edges_coef)) {
          for (from_idx in 1:nrow(edges_coef)) {
            for (to_idx in 1:ncol(edges_coef)) {
              coef_value <- edges_coef[from_idx, to_idx]
              
              # Only process non-zero coefficients
              if (coef_value != 0) {
                from_var <- rownames(edges_coef)[from_idx]
                to_var <- colnames(edges_coef)[to_idx]
                
                # Create column names
                edge_name <- paste0(from_var, "_to_", to_var)
                est_col <- paste0(edge_name, "_est")
                se_col <- paste0(edge_name, "_se")
                
                # Add values to model row
                model_row[[est_col]] <- coef_value
                
                if (!is.null(edges_se)) {
                  se_value <- edges_se[from_idx, to_idx]
                  model_row[[se_col]] <- se_value
                  
                  # Calculate p-value (assuming normal distribution)
                  if (se_value > 0) {
                    z_score <- abs(coef_value / se_value)
                    p_value <- 2 * (1 - pnorm(z_score))
                    p_col <- paste0(edge_name, "_p")
                    model_row[[p_col]] <- p_value
                  }
                }
              }
            }
          }
        }
        
        # Add to detailed models list
        detailed_models_list[[length(detailed_models_list) + 1]] <- model_row
      }
      
      # Print progress
      if (i %% 10 == 0) {
        elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
        estimated_total <- elapsed_time * n_species / i
        cat(sprintf("Progress: %d/%d (%.1f%%) - Elapsed: %.1f min - Est. total: %.1f min\n",
                   i, n_species, 100*i/n_species, elapsed_time, estimated_total))
      }
      
    }, error = function(e) {
      cat(sprintf("Error processing species '%s': %s\n", species_to_remove, e$message))
      results_df$seed[i] <- current_seed
      results_df$RemovedSpecies[i] <- species_to_remove
      results_df$nSpecies[i] <- NA
    })
  }
  
  # Combine detailed models into a single data frame
  detailed_models_df <- bind_rows(detailed_models_list)
  
  # Save the detailed models data frame
  if (!is.null(detailed_models_df) && nrow(detailed_models_df) > 0) {
    filename <- paste0("detailed_models_JackknifeSpecies_n", n_species, "_", Sys.Date(), ".csv")
    write.csv(detailed_models_df, file.path(output_dir, filename), row.names = FALSE)
    cat(sprintf("Detailed model information saved to %s\n", filename))
  }
  
  # Calculate model frequency
  all_models <- unlist(strsplit(results_df$CICsub2_models_string, ", "))
  model_freq_df <- as.data.frame(table(all_models))
  names(model_freq_df) <- c("model_name", "frequency")
  model_freq_df$proportion <- model_freq_df$frequency / n_species
  model_freq_df <- model_freq_df %>% arrange(desc(frequency))
  
  # Save model frequencies
  freq_filename <- paste0("model_frequencies_JackknifeSpecies_n", n_species, "_", Sys.Date(), ".csv")
  write.csv(model_freq_df, file.path(output_dir, freq_filename), row.names = FALSE)
  
  # Save main results
  results_filename <- paste0("results_JackknifeSpecies_n", n_species, "_", Sys.Date(), ".csv")
  write.csv(results_df, file.path(output_dir, results_filename), row.names = FALSE)
  
  # Save conditional average plots if requested
  if (save_conditional_plots && length(conditional_average_plots) > 0) {
    pdf_name <- file.path(output_dir, paste0("conditional_average_plots_JackknifeSpecies_n", n_species, "_", Sys.Date(), ".pdf"))
    
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
      do.call(gridExtra::grid.arrange, c(page_plots, ncol = 2))
    }
    
    dev.off()
    cat(sprintf("Conditional average plots saved to %s\n", pdf_name))
  }
  
  # Create edge summary if path coefficients were saved
  if (save_path_coefficients && !is.null(detailed_models_df)) {
    est_cols <- grep("_est$", names(detailed_models_df), value = TRUE)
    
    edge_summary <- data.frame(
      edge = gsub("_est$", "", est_cols),
      stringsAsFactors = FALSE
    )
    
    for (col in est_cols) {
      p_col <- gsub("_est$", "_p", col)
      
      if (p_col %in% names(detailed_models_df)) {
        valid_data <- detailed_models_df[!is.na(detailed_models_df[[col]]), ]
        
        if (nrow(valid_data) > 0) {
          edge_index <- which(edge_summary$edge == gsub("_est$", "", col))
          
          edge_summary$n_occurrences[edge_index] <- nrow(valid_data)
          edge_summary$mean_estimate[edge_index] <- mean(valid_data[[col]], na.rm = TRUE)
          edge_summary$median_estimate[edge_index] <- median(valid_data[[col]], na.rm = TRUE)
          edge_summary$n_significant[edge_index] <- sum(valid_data[[p_col]] < 0.05, na.rm = TRUE)
          edge_summary$prop_significant[edge_index] <- edge_summary$n_significant[edge_index] / edge_summary$n_occurrences[edge_index]
          edge_summary$mean_p[edge_index] <- mean(valid_data[[p_col]], na.rm = TRUE)
        }
      }
    }
    
    edge_summary <- edge_summary %>% arrange(desc(prop_significant))
    
    # Save edge summary
    edge_filename <- paste0("edge_summary_JackknifeSpecies_n", n_species, "_", Sys.Date(), ".csv")
    write.csv(edge_summary, file.path(output_dir, edge_filename), row.names = FALSE)
  }
  
  # Print final summary
  total_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  cat(sprintf("\nJackknife analysis complete!\n"))
  cat(sprintf("Total time: %.1f minutes\n", total_time))
  cat(sprintf("Species analyzed: %d\n", n_species))
  cat(sprintf("Results saved to: %s\n", output_dir))
  
  # Return results
  return(list(
    results_df = results_df,
    model_frequencies = model_freq_df,
    detailed_models = detailed_models_df,
    edge_summary = edge_summary
  ))
}