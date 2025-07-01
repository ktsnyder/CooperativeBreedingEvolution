## Helper functions for simplifying simmap overlap runner scripts
## Kate Snyder / Claude
## 2025-06-28

# Load transition count files in a standardized way
loadTransitionCounts <- function(trait1, trait2, path = "Simmap Overlap Outputs", 
                                real_pattern = "REAL", dummy_pattern = "DUMMY") {
  require(stringr)
  
  # List all CSV files in the directory
  filelist <- list.files(path = path, pattern = "overlap_counts.*\\.csv$", full.names = TRUE)
  
  # Filter for the specific trait pair
  filelist <- filelist[str_detect(filelist, trait1) & str_detect(filelist, trait2)]
  
  # Separate real and dummy files
  real_file <- filelist[str_detect(filelist, real_pattern) & !str_detect(filelist, dummy_pattern)]
  dummy_file <- filelist[str_detect(filelist, dummy_pattern)]
  
  if (length(real_file) == 0) {
    stop(paste("No REAL data file found for", trait1, "and", trait2))
  }
  if (length(dummy_file) == 0) {
    stop(paste("No DUMMY data file found for", trait1, "and", trait2))
  }
  
  # Read the files
  df_real <- read.csv(real_file[1])
  df_dummy <- read.csv(dummy_file[1])
  
  # Get nsims
  nsims_real <- nrow(df_real)
  nsims_dummy <- nrow(df_dummy)
  
  cat("Loaded", trait1, "vs", trait2, "data:\n")
  cat("  Real file:", basename(real_file[1]), "(", nsims_real, "sims)\n")
  cat("  Dummy file:", basename(dummy_file[1]), "(", nsims_dummy, "sims)\n")
  
  return(list(
    df_real = df_real,
    df_dummy = df_dummy,
    nsims_real = nsims_real,
    nsims_dummy = nsims_dummy,
    files = list(real = real_file[1], dummy = dummy_file[1])
  ))
}

# Process transition statistics in a standardized way
processTransitionStats <- function(df_real, df_dummy = NULL, trait1, trait2) {
  require(stringr)
  require(emmeans)
  
  # If we have transition columns, process them
  transitioncols <- grep("^trait[12]_[01]to[01]_in_trait[12]_[01]$", names(df_real), value = TRUE)
  
  # If no generic columns, look for legacy columns
  if (length(transitioncols) == 0) {
    transitioncols <- c("Coop0to1inFS0", "Coop1to0inFS0", "Coop0to1inFS1", "Coop1to0inFS1",
                       "FS0to1inCoop0", "FS1to0inCoop0", "FS0to1inCoop1", "FS1to0inCoop1")
    transitioncols <- intersect(transitioncols, names(df_real))
  }
  
  if (length(transitioncols) == 0) {
    warning("No transition columns found in data")
    return(NULL)
  }
  
  ExpectedCols <- paste0(transitioncols, "Expected")
  
  # Calculate chi-square statistics
  ChiSqStat <- rowSums((df_real[,transitioncols] - df_real[,ExpectedCols])^2 / 
                      df_real[,ExpectedCols])
  ChiSqPvals <- pchisq(ChiSqStat, df = 7)
  df_real <- cbind(df_real, ChiSqStat, ChiSqPvals)
  
  # Calculate statistics using calcHuel if both real and dummy provided
  if (!is.null(df_dummy)) {
    source("calcHuel_corrected.R")
    calcHuelout <- calcHuel_corrected(df_real, df_dummy, trait1_name = trait1, trait2_name = trait2)
    
    if (!is.null(calcHuelout$TransitionStats)) {
      return(calcHuelout$TransitionStats)
    }
  }
  
  # Otherwise do basic transition stats
  hist(df_real$ChiSqPvals, main = paste(trait1, "vs", trait2), breaks = 20)
  abline(v = 0.05, col = "red")
  
  nsims <- nrow(df_real)
  
  # Count times actual > expected
  timesActualGreaterThanExpected <- df_real[, transitioncols] > df_real[, ExpectedCols]
  NtimesActualGreaterThanExpected <- colSums(timesActualGreaterThanExpected, na.rm = TRUE)
  
  return(list(
    df_with_chisq = df_real,
    NtimesActualGreaterThanExpected = NtimesActualGreaterThanExpected,
    nsims = nsims
  ))
}

# Prepare rate p-values for transition plots
prepareRatePvals <- function(TransitionStats, nsims) {
  require(dplyr)
  
  if (is.null(TransitionStats$logPairwisePostHoc)) {
    warning("No pairwise post-hoc results found")
    return(NULL)
  }
  
  # Extract transitions and p-values
  Transitions <- TransitionStats$logPairwisePostHoc$Transition
  pvals <- TransitionStats$logPairwisePostHoc$p.value
  
  pvaldf <- data.frame(
    Transition = Transitions,
    p.value = pvals,
    stringsAsFactors = FALSE
  )
  
  # Add significance labels
  pvaldf$SignificanceLabel <- "n.s."
  pvaldf$SignificanceLabel[pvaldf$p.value < 0.05] <- "p < 0.05"
  pvaldf$SignificanceLabel[pvaldf$p.value < 0.01] <- "p < 0.01"
  pvaldf$SignificanceLabel[pvaldf$p.value < 0.001] <- "p < 0.001"
  pvaldf$SignificanceLabel[pvaldf$p.value < 0.0001] <- "p < 0.0001"
  
  # Add N times actual > expected
  NtimesActualGreaterThanExpected <- TransitionStats$NtimesActualGreaterThanExpected
  NtimesActualGreaterThanExpecteddf <- data.frame(
    Transition = names(NtimesActualGreaterThanExpected),
    CountNActualGreaterThanExpected = as.integer(NtimesActualGreaterThanExpected),
    stringsAsFactors = FALSE
  )
  
  pvaldf <- merge(pvaldf, NtimesActualGreaterThanExpecteddf, by = "Transition")
  pvaldf$FractionActualGreaterThanExpected <- pvaldf$CountNActualGreaterThanExpected / nsims
  
  # Add percent trending labels
  pvaldf$PercentTrendingLabel <- "<80%"
  pvaldf$PercentTrendingLabel[pvaldf$FractionActualGreaterThanExpected > .8 | 
                              pvaldf$FractionActualGreaterThanExpected < .2] <- ">80%"
  pvaldf$PercentTrendingLabel[pvaldf$FractionActualGreaterThanExpected > .9 | 
                              pvaldf$FractionActualGreaterThanExpected < .1] <- ">90%"
  pvaldf$PercentTrendingLabel[pvaldf$FractionActualGreaterThanExpected > .95 | 
                              pvaldf$FractionActualGreaterThanExpected < .05] <- ">95%"
  pvaldf$PercentTrendingLabel[pvaldf$FractionActualGreaterThanExpected > .99 | 
                              pvaldf$FractionActualGreaterThanExpected < .01] <- ">99%"
  
  # Map to q rates if we have generic names
  if (grepl("trait[12]_", pvaldf$Transition[1])) {
    source("transition_count_helpers.R")
    pvaldf$qRate <- mapGenericToArrowNames(pvaldf$Transition)
    
    # Filter out any NA mappings
    ratePvals <- pvaldf[!is.na(pvaldf$qRate), ]
    names(ratePvals)[names(ratePvals) == "Transition"] <- "Transitions"
  } else {
    # Legacy mapping
    RateRef <- data.frame(
      Transitions = c("FS0to1inCoop0", "Coop0to1inFS0", "FS1to0inCoop0", "Coop0to1inFS1",
                     "Coop1to0inFS0", "FS0to1inCoop1", "Coop1to0inFS1", "FS1to0inCoop1"),
      qRate = c("q12", "q13", "q21", "q24", "q31", "q34", "q42", "q43"),
      stringsAsFactors = FALSE
    )
    ratePvals <- merge(RateRef, pvaldf, by.x = "Transitions", by.y = "Transition")
  }
  
  return(ratePvals)
}

# Run full simmap overlap analysis
runSimmapOverlapAnalysis <- function(trait1, trait2, 
                                   nsims_real = 500, nsims_dummy = 500,
                                   tree_file = NULL, data_file = NULL,
                                   calculate_transitions = TRUE,
                                   plot_transitions = TRUE,
                                   save_outputs = TRUE,
                                   output_dir = "Simmap Overlap Outputs",
                                   use_pregenerated_simmaps = FALSE,
                                   trait1_real_simmaps_file = NULL,
                                   trait2_real_simmaps_file = NULL,
                                   trait1_dummy_simmaps_file = NULL,
                                   trait2_dummy_simmaps_file = NULL) {
  
  require(phytools)
  source("CharacterSimmaps_modified.R")
  source("calcHuel_corrected.R")
  
  # Load tree and data if provided
  if (!is.null(tree_file)) {
    if (is.character(tree_file)) {
      tree <- read.nexus(tree_file)
    } else {
      tree <- tree_file
    }
  } else {
    stop("Tree file must be provided")
  }
  
  if (!is.null(data_file)) {
    if (is.character(data_file)) {
      df <- read.csv(data_file)
    } else {
      df <- data_file
    }
  } else {
    stop("Data file must be provided")
  }
  
  columns <- c(trait1, trait2)
  
  cat("Running simmap overlap analysis for:", trait1, "vs", trait2, "\n")
  
  # Run real data simmaps
  if (use_pregenerated_simmaps) {
    dfout_result <- CharacterSimmaps_modified(
      columns = columns, df = df, tree = tree,
      dummy = FALSE, nsims = nsims_real,
      treelabel = "tree", datalabel = NULL,
      output_dir = output_dir,
      return_simmaps = calculate_transitions,
      trait1_real_multisimmapRDS = trait1_real_simmaps_file,
      trait2_real_multisimmapRDS = trait2_real_simmaps_file
    )
  } else {
    dfout_result <- CharacterSimmaps_modified(
      columns = columns, df = df, tree = tree,
      dummy = FALSE, nsims = nsims_real,
      treelabel = "tree", datalabel = NULL,
      output_dir = output_dir,
      return_simmaps = calculate_transitions
    )
  }
  
  # Extract dataframe and simmaps
  if (is.list(dfout_result) && "dfout" %in% names(dfout_result)) {
    dfout <- dfout_result$dfout
    trait1_real_simmaps <- dfout_result$trait1_simmaps
    trait2_real_simmaps <- dfout_result$trait2_simmaps
  } else {
    dfout <- dfout_result
    trait1_real_simmaps <- NULL
    trait2_real_simmaps <- NULL
  }
  
  # Run dummy data simmaps
  if (use_pregenerated_simmaps) {
    dfDummy_result <- CharacterSimmaps_modified(
      columns = columns, df = df, tree = tree,
      dummy = TRUE, nsims = nsims_dummy,
      treelabel = "tree", datalabel = NULL,
      dummyMethod = "makeSimmap",
      output_dir = output_dir,
      return_simmaps = FALSE,  # Don't need dummy simmaps for transitions
      trait1_dummy_multisimmapRDS = trait1_dummy_simmaps_file,
      trait2_dummy_multisimmapRDS = trait2_dummy_simmaps_file
    )
  } else {
    dfDummy_result <- CharacterSimmaps_modified(
      columns = columns, df = df, tree = tree,
      dummy = TRUE, nsims = nsims_dummy,
      treelabel = "tree", datalabel = NULL,
      dummyMethod = "makeSimmap",
      output_dir = output_dir,
      return_simmaps = FALSE
    )
  }
  
  # Extract dummy dataframe
  if (is.list(dfDummy_result) && "dfout" %in% names(dfDummy_result)) {
    dfDummy <- dfDummy_result$dfout
  } else {
    dfDummy <- dfDummy_result
  }
  
  # Run calcHuel analysis
  calcHuelout <- calcHuel_corrected(
    dfout = dfout,
    dfDummy = dfDummy,
    trait1_name = trait1,
    trait2_name = trait2,
    calculate_transitions = calculate_transitions,
    plot_transitions = plot_transitions,
    trait1_simmaps = trait1_real_simmaps,
    trait2_simmaps = trait2_real_simmaps,
    plot_ggplots_pdf = save_outputs
  )
  
  # Save combined plots if requested
  if (save_outputs && !is.null(calcHuelout$plots)) {
    require(gridExtra)
    require(grid)
    
    # Determine filename based on whether transitions were calculated
    if (plot_transitions && !is.null(calcHuelout$TransitionStats)) {
      plotname_base <- file.path(output_dir, 
                                paste(trait1, trait2, nsims_real, nsims_dummy, "withTransCounts"))
    } else {
      plotname_base <- file.path(output_dir, 
                                paste(trait1, trait2, nsims_real, nsims_dummy))
    }
    
    # Get all plots
    all_plots <- calcHuelout$plots
    nPlots <- length(all_plots)
    
    if (nPlots > 0) {
      # Method 1: Use marrangeGrob without the arrow plot to avoid issues
      plots_for_pdf <- all_plots
      
      # Check if p4 or p4_gray exists
      p4_plot <- NULL
      p4_gray_plot <- NULL
      if ("p4" %in% names(plots_for_pdf)) {
        p4_plot <- plots_for_pdf$p4
        plots_for_pdf$p4 <- NULL
      }
      if ("p4_gray" %in% names(plots_for_pdf)) {
        p4_gray_plot <- plots_for_pdf$p4_gray
        plots_for_pdf$p4_gray <- NULL
      }
      
      # Save combined PDF without blank pages
      pdf_filename <- paste0(plotname_base, ".pdf")
      
      # Filter out NULL plots
      valid_plots <- list()
      for (pname in names(plots_for_pdf)) {
        p <- plots_for_pdf[[pname]]
        if (!is.null(p) && inherits(p, "ggplot")) {
          valid_plots[[pname]] <- p
        }
      }
      
      if (length(valid_plots) > 0) {
        # Use marrangeGrob to arrange plots without blank pages
        ml <- marrangeGrob(grobs = valid_plots, 
                          nrow = 1, ncol = 1,  # One plot per page
                          top = NULL)
        
        ggsave(pdf_filename, ml, width = 7.5, height = 3.8)
        cat("Combined plots saved to:", pdf_filename, " (", length(valid_plots), " plots)\n")
      }
      
      # Save as JPEG for viewing
      if (length(plots_for_pdf) > 0) {
        jpeg_filename <- paste0(plotname_base, ".jpg")
        jpeg(jpeg_filename, width = 7.5 * 100, height = 3.8 * length(plots_for_pdf) * 100, 
             quality = 95, units = "px")
        
        # Create a grid arrangement of the plots
        if (length(plots_for_pdf) == 1) {
          grid.draw(plots_for_pdf[[1]])
        } else {
          m3 <- marrangeGrob(grobs = plots_for_pdf, nrow = length(plots_for_pdf), ncol = 1, top = NULL)
          grid.draw(m3)
        }
        
        dev.off()
        cat("JPEG version saved to:", jpeg_filename, "\n")
      }
      
      # Save arrow plots separately if they exist
      if (!is.null(p4_plot)) {
        p4_filename <- paste0(plotname_base, "_arrow_plot.pdf")
        pdf(p4_filename, width = 10, height = 8)
        
        # Try different methods to print p4
        tryCatch({
          print(p4_plot)
        }, error = function(e) {
          cat("Error printing p4 as ggplot:", e$message, "\n")
          # Try as grob
          if (inherits(p4_plot, "grob")) {
            grid.draw(p4_plot)
          }
        })
        
        dev.off()
        cat("Arrow plot saved separately to:", p4_filename, "\n")
        
        # Also save as JPEG
        p4_jpeg <- paste0(plotname_base, "_arrow_plot.jpg")
        jpeg(p4_jpeg, width = 1000, height = 800, quality = 95)
        try(print(p4_plot))
        dev.off()
      }
      
      # Save gray arrow plot if it exists
      if (!is.null(p4_gray_plot)) {
        p4_gray_filename <- paste0(plotname_base, "_arrow_plot_gray_ns.pdf")
        pdf(p4_gray_filename, width = 10, height = 8)
        
        tryCatch({
          print(p4_gray_plot)
        }, error = function(e) {
          cat("Error printing p4_gray as ggplot:", e$message, "\n")
          if (inherits(p4_gray_plot, "grob")) {
            grid.draw(p4_gray_plot)
          }
        })
        
        dev.off()
        cat("Gray non-significant arrow plot saved to:", p4_gray_filename, "\n")
        
        # Also save as JPEG
        p4_gray_jpeg <- paste0(plotname_base, "_arrow_plot_gray_ns.jpg")
        jpeg(p4_gray_jpeg, width = 1000, height = 800, quality = 95)
        try(print(p4_gray_plot))
        dev.off()
      }
    }
  }
  
  return(calcHuelout)
}

# Helper to extract trait labels using getLabels
extractTraitLabels <- function(trait1, trait2) {
  # Try to source getLabels function
  if (file.exists("test_trait_overlap_simmaps.R")) {
    source("test_trait_overlap_simmaps.R")
    trait1_labels <- getLabels(trait1)
    trait2_labels <- getLabels(trait2)
  } else {
    # Use local version
    source("transition_count_helpers.R")
    trait1_labels <- getLabels_local(trait1)
    trait2_labels <- getLabels_local(trait2)
  }
  
  return(list(
    trait1_labels = trait1_labels,
    trait2_labels = trait2_labels
  ))
}