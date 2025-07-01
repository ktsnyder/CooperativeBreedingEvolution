## Corrected calcHuel function with proper state indexing
## Kate Snyder
## 2025-06-26
## Fixes the set.seed bug and uses generic column names

calcHuel_corrected <- function(dfout, dfDummy, nsims_real = NULL, nsims_dummy = NULL, 
                              otherlabel = NULL, newplot = TRUE, plot_ggplots_pdf = FALSE,
                              trait1_name = "Trait1", trait2_name = "Trait2",
                              calculate_transitions = FALSE, plot_transitions = FALSE,
                              trait1_simmaps = NULL, trait2_simmaps = NULL,
                              dirs = NULL, debug = FALSE) {
  require(tidyverse)
  require(ggplot2)
  
  # Source debug control if not already loaded
  if (!exists("debugPrint")) {
    if (file.exists("debug_control.R")) {
      source("debug_control.R")
    } else {
      # Fallback debug functions if debug_control.R not available
      debugPrint <- function(...) { if (debug) cat("DEBUG:", ..., "\n") }
      validationPrint <- function(...) { cat("VALIDATION:", ..., "\n") }
    }
  }
  
  # Set debug mode for this run
  if (exists("setDebugMode")) {
    setDebugMode(debug)
  }
  
  # Source file output helpers if not already loaded
  if (!exists("createStandardizedFilename")) {
    source("file_output_helpers.R")
  }
  
  # Use provided dirs or create simple list for backward compatibility
  if (is.null(dirs)) {
    dirs <- list(
      main = ".",
      simmaps = ".",
      data = ".",
      plots = ".",
      summaries = "."
    )
  }

 # if (is.null(nsims_real)) {
    nsims_real = length(dfout[,1])
#  }
 # if (is.null(nsims_dummy)) {
    nsims_dummy = length(dfDummy[,1])
  #}
  
  Nspecies = dfout[1,"Nspecies"]
  
  # Use provided trait names or extract from dataframe
  if ("trait1" %in% names(dfout)) {
    trait1 = dfout[1,"trait1"]
    trait2 = dfout[1,"trait2"]
  } else {
    trait1 = trait1_name
    trait2 = trait2_name
  }
  
  # Create standardized PDF name
  PDFname <- createStandardizedFilename(trait1, trait2, "combined_analysis",
                                       nsims_real, nsims_dummy,
                                       other_label = otherlabel,
                                       extension = "pdf")
  PDFpath <- file.path(dirs$plots, PDFname)
  
  # Calculate expected proportions assuming independence
  # Using generic column names
  if ("prop_trait1_state0" %in% names(dfout)) {
    # New format with generic names
    # Note: Expected proportions follow original order (trait2 * trait1) to match Map.Overlap output # this is ok here, but otherwise will be fixing this
    ExpProp_0_0 <- as.numeric(dfout$prop_trait2_state0) * as.numeric(dfout$prop_trait1_state0)
    ExpProp_1_0 <- as.numeric(dfout$prop_trait2_state0) * as.numeric(dfout$prop_trait1_state1)
    ExpProp_0_1 <- as.numeric(dfout$prop_trait2_state1) * as.numeric(dfout$prop_trait1_state0)
    ExpProp_1_1 <- as.numeric(dfout$prop_trait2_state1) * as.numeric(dfout$prop_trait1_state1)
    
    # Observed proportions
    ObsProp_0_0 <- as.numeric(dfout$Overlap_0_0)
    ObsProp_0_1 <- as.numeric(dfout$Overlap_0_1)
    ObsProp_1_0 <- as.numeric(dfout$Overlap_1_0)
    ObsProp_1_1 <- as.numeric(dfout$Overlap_1_1)
  } else {
    # Old format for compatibility
    ExpProp_0_0 <- as.numeric(dfout$propFSabsent) * as.numeric(dfout$propNoncoop)
    ExpProp_0_1 <- as.numeric(dfout$propFSabsent) * as.numeric(dfout$propCoop)
    ExpProp_1_0 <- as.numeric(dfout$propFSpresent) * as.numeric(dfout$propNoncoop)
    ExpProp_1_1 <- as.numeric(dfout$propFSpresent) * as.numeric(dfout$propCoop)
    
    ObsProp_0_0 <- as.numeric(dfout$ObsProp0Absent)
    ObsProp_0_1 <- as.numeric(dfout$ObsProp1Absent)
    ObsProp_1_0 <- as.numeric(dfout$ObsProp0Present)
    ObsProp_1_1 <- as.numeric(dfout$ObsProp1Present)
  }
  
  # Calculate D_real: sum of absolute deviations across all simmaps
  d_0_0_sum <- sum(abs(ObsProp_0_0 - ExpProp_0_0))
  d_0_1_sum <- sum(abs(ObsProp_0_1 - ExpProp_0_1))
  d_1_0_sum <- sum(abs(ObsProp_1_0 - ExpProp_1_0))
  d_1_1_sum <- sum(abs(ObsProp_1_1 - ExpProp_1_1))
  
  D_real <- (d_0_0_sum + d_0_1_sum + d_1_0_sum + d_1_1_sum) / nsims_real
  cat("D_real:", D_real, "\n")
  
  # Calculate distribution of D values for real data (per simmap)
  d_0_0 <- abs(ObsProp_0_0 - ExpProp_0_0)
  d_0_1 <- abs(ObsProp_0_1 - ExpProp_0_1)
  d_1_0 <- abs(ObsProp_1_0 - ExpProp_1_0)
  d_1_1 <- abs(ObsProp_1_1 - ExpProp_1_1)
  Real_dsims <- rowSums(cbind(d_0_0, d_0_1, d_1_0, d_1_1))
  
  # Process dummy data
  if ("prop_trait1_state0" %in% names(dfDummy)) {
    # New format
    # Note: Expected proportions follow original order (trait2 * trait1) to match Map.Overlap output
    ExpProp_0_0 <- as.numeric(dfDummy$prop_trait2_state0) * as.numeric(dfDummy$prop_trait1_state0)
    ExpProp_1_0 <- as.numeric(dfDummy$prop_trait2_state0) * as.numeric(dfDummy$prop_trait1_state1)
    ExpProp_0_1 <- as.numeric(dfDummy$prop_trait2_state1) * as.numeric(dfDummy$prop_trait1_state0)
    ExpProp_1_1 <- as.numeric(dfDummy$prop_trait2_state1) * as.numeric(dfDummy$prop_trait1_state1)
    
    ObsProp_0_0 <- as.numeric(dfDummy$Overlap_0_0)
    ObsProp_0_1 <- as.numeric(dfDummy$Overlap_0_1)
    ObsProp_1_0 <- as.numeric(dfDummy$Overlap_1_0)
    ObsProp_1_1 <- as.numeric(dfDummy$Overlap_1_1)
  } else {
    # Old format
    ExpProp_0_0 <- as.numeric(dfDummy$propFSabsent) * as.numeric(dfDummy$propNoncoop)
    ExpProp_0_1 <- as.numeric(dfDummy$propFSabsent) * as.numeric(dfDummy$propCoop)
    ExpProp_1_0 <- as.numeric(dfDummy$propFSpresent) * as.numeric(dfDummy$propNoncoop)
    ExpProp_1_1 <- as.numeric(dfDummy$propFSpresent) * as.numeric(dfDummy$propCoop)
    
    ObsProp_0_0 <- as.numeric(dfDummy$ObsProp0Absent)
    ObsProp_0_1 <- as.numeric(dfDummy$ObsProp1Absent)
    ObsProp_1_0 <- as.numeric(dfDummy$ObsProp0Present)
    ObsProp_1_1 <- as.numeric(dfDummy$ObsProp1Present)
  }
  
  d_0_0 <- abs(ObsProp_0_0 - ExpProp_0_0)
  d_0_1 <- abs(ObsProp_0_1 - ExpProp_0_1)
  d_1_0 <- abs(ObsProp_1_0 - ExpProp_1_0)
  d_1_1 <- abs(ObsProp_1_1 - ExpProp_1_1)
  Dummy_dsums <- rowSums(cbind(d_0_0, d_0_1, d_1_0, d_1_1))
  
  # Calculate p-value
  pval <- sum(Dummy_dsums > D_real) / nsims_dummy
  cat("Number of Dummy D values > D_real:", sum(Dummy_dsums > D_real), "\n")
  cat("P-value:", pval, "\n")
  
  # Calculate fraction of dummy simulations less than median real for each state
  if ("prop_trait1_state0" %in% names(dfout)) {
    # Using new column names
    median_real_0_0 <- median(as.numeric(dfout$Overlap_0_0))
    median_real_0_1 <- median(as.numeric(dfout$Overlap_0_1))
    median_real_1_0 <- median(as.numeric(dfout$Overlap_1_0))
    median_real_1_1 <- median(as.numeric(dfout$Overlap_1_1))
    
    Overlap_0_0_FractionDummyLessThanMedianReal <- sum(as.numeric(dfDummy$Overlap_0_0) <= median_real_0_0) / nsims_dummy
    Overlap_0_1_FractionDummyLessThanMedianReal <- sum(as.numeric(dfDummy$Overlap_0_1) <= median_real_0_1) / nsims_dummy
    Overlap_1_0_FractionDummyLessThanMedianReal <- sum(as.numeric(dfDummy$Overlap_1_0) <= median_real_1_0) / nsims_dummy
    Overlap_1_1_FractionDummyLessThanMedianReal <- sum(as.numeric(dfDummy$Overlap_1_1) <= median_real_1_1) / nsims_dummy
    
    # Diagnostic output
    cat("\nFraction calculation diagnostics:\n")
    cat("Median real overlaps: 0_0=", median_real_0_0, " 0_1=", median_real_0_1, 
        " 1_0=", median_real_1_0, " 1_1=", median_real_1_1, "\n")
    cat("Sum of median real overlaps:", median_real_0_0 + median_real_0_1 + median_real_1_0 + median_real_1_1, "\n")
    
    # Check dummy overlap ranges
    cat("\nDummy overlap ranges:\n")
    cat("  0_0: min=", min(as.numeric(dfDummy$Overlap_0_0)), " max=", max(as.numeric(dfDummy$Overlap_0_0)), "\n")
    cat("  0_1: min=", min(as.numeric(dfDummy$Overlap_0_1)), " max=", max(as.numeric(dfDummy$Overlap_0_1)), "\n")
    cat("  1_0: min=", min(as.numeric(dfDummy$Overlap_1_0)), " max=", max(as.numeric(dfDummy$Overlap_1_0)), "\n")
    cat("  1_1: min=", min(as.numeric(dfDummy$Overlap_1_1)), " max=", max(as.numeric(dfDummy$Overlap_1_1)), "\n")
    
    cat("\nFractions dummy <= median real: 0_0=", Overlap_0_0_FractionDummyLessThanMedianReal,
        " 0_1=", Overlap_0_1_FractionDummyLessThanMedianReal,
        " 1_0=", Overlap_1_0_FractionDummyLessThanMedianReal,
        " 1_1=", Overlap_1_1_FractionDummyLessThanMedianReal, "\n")
    cat("Sum of fractions (should not equal 1 or 2):", 
        Overlap_0_0_FractionDummyLessThanMedianReal + Overlap_0_1_FractionDummyLessThanMedianReal +
        Overlap_1_0_FractionDummyLessThanMedianReal + Overlap_1_1_FractionDummyLessThanMedianReal, "\n")
  } else {
    # Using old column names for compatibility
    median_real_0_0 <- median(as.numeric(dfout$ObsProp0Absent))
    median_real_0_1 <- median(as.numeric(dfout$ObsProp1Absent))
    median_real_1_0 <- median(as.numeric(dfout$ObsProp0Present))
    median_real_1_1 <- median(as.numeric(dfout$ObsProp1Present))
    
    Overlap_0_0_FractionDummyLessThanMedianReal <- sum(as.numeric(dfDummy$ObsProp0Absent) <= median_real_0_0) / nsims_dummy
    Overlap_0_1_FractionDummyLessThanMedianReal <- sum(as.numeric(dfDummy$ObsProp1Absent) <= median_real_0_1) / nsims_dummy
    Overlap_1_0_FractionDummyLessThanMedianReal <- sum(as.numeric(dfDummy$ObsProp0Present) <= median_real_1_0) / nsims_dummy
    Overlap_1_1_FractionDummyLessThanMedianReal <- sum(as.numeric(dfDummy$ObsProp1Present) <= median_real_1_1) / nsims_dummy
  }
  
  # Diagnostic output
  cat("\nDiagnostic statistics:\n")
  cat("Real D values - mean:", mean(Real_dsims), "sd:", sd(Real_dsims), "CV:", sd(Real_dsims)/mean(Real_dsims), "\n")
  cat("Dummy D values - mean:", mean(Dummy_dsums), "sd:", sd(Dummy_dsums), "CV:", sd(Dummy_dsums)/mean(Dummy_dsums), "\n")
  cat("Real D values range:", range(Real_dsims), "\n")
  cat("Dummy D values range:", range(Dummy_dsums), "\n")
  
  # Plot setup
  plotlabel = paste(trait1, "vs", trait2, "\nN species =", Nspecies, otherlabel)
  dummytitle = paste("Nsims =", nsims_dummy, otherlabel, "\nnum Dummy D > D_real:", sum(Dummy_dsums > D_real), ", pval =", pval)
  
  xmax = max(c(Real_dsims, Dummy_dsums)) * 1.1
  
  # Plot histograms
  if (newplot == TRUE) {
    font_size <- 1.0
    par(mfrow=c(2,1), mai=c(0.8,1.0,0.5,0.3), oma=c(2,3,2,3), 
        font.main=1, cex.main=1.25, cex.lab=font_size, cex.axis=1)
  }
  
  # Histogram for Real_dsims
  hist(Real_dsims, xlim=c(0,xmax), breaks=seq(0, xmax, length.out=21), 
       col=rgb(0.2,0.5,0.7,0.5),
       main=plotlabel, xlab="D statistic from real data simmaps", ylab="Frequency",
       border="white")
  abline(v = D_real, col="red", lwd=2.5)
  
  # Histogram for Dummy_dsums
  hist(Dummy_dsums, xlim=c(0,xmax), breaks=seq(0, xmax, length.out=21), 
       col=rgb(0.7,0.5,0.2,0.5),
       main=dummytitle, xlab="D statistic from simulated independent data simmaps", ylab="Frequency",
       border="white")
  
  # Create ggplot versions
  real_data <- data.frame(value = Real_dsims)
  dummy_data <- data.frame(value = Dummy_dsums)
  
  p1 <- ggplot(real_data, aes(x = value)) +
    geom_histogram(bins = 20, fill = rgb(0.2, 0.5, 0.7, 0.5), color = "white") +
    geom_vline(xintercept = D_real, color = "red", size = 1.5) +
    xlim(c(0, xmax)) +
    labs(title = plotlabel, x = "D statistic from real data simmaps", y = "Frequency") +
    theme_minimal(base_size = 10)

  p2 <- ggplot(dummy_data, aes(x = value)) +
    geom_histogram(bins = 20, fill = rgb(0.7, 0.5, 0.2, 0.5), color = "white") +
    xlim(c(0, xmax)) +
    labs(title = dummytitle, x = "D statistic from simulated independent data simmaps", y = "Frequency") +
    theme_minimal(base_size = 10)
  
  # Create boxplot for observed state proportions
  if ("trait1" %in% names(dfout) && "trait2" %in% names(dfout)) {
    # Create melted dataframes for boxplot
    library(tidyr)
    library(dplyr)
    
    # dfRealMelt <- dfout %>%
    #   select(any_of(c("treenum", "ObsProp0Absent", "ObsProp0Present", "ObsProp1Absent", "ObsProp1Present"))) %>%
    #   gather("ObservedState", "ObservedState.prop", ObsProp0Absent:ObsProp1Present) %>%
    #   mutate(Which = "Real")
    
    dfRealMelt <- dfout %>%
      select(any_of(c("treenum", "Overlap_0_0", "Overlap_0_1", "Overlap_1_0", "Overlap_1_1"))) %>%
      gather("ObservedState", "ObservedState.prop", Overlap_0_0:Overlap_1_1) %>%
      mutate(Which = "Real")
    
    # dfDummyMelt <- dfDummy %>%
    #   select(any_of(c("treenum", "ObsProp0Absent", "ObsProp0Present", "ObsProp1Absent", "ObsProp1Present"))) %>%
    #   gather("ObservedState", "ObservedState.prop", ObsProp0Absent:ObsProp1Present) %>%
    #   mutate(Which = "Dummy")
    dfDummyMelt <- dfDummy %>%
      select(any_of(c("treenum", "Overlap_0_0", "Overlap_0_1", "Overlap_1_0", "Overlap_1_1"))) %>%
      gather("ObservedState", "ObservedState.prop", Overlap_0_0:Overlap_1_1) %>%
      mutate(Which = "Dummy")
    
    
    dfCombined <- rbind(dfRealMelt, dfDummyMelt)
    dfCombined$ObservedState.prop <- as.numeric(dfCombined$ObservedState.prop)
    
    # Get trait labels
    source("getLabels.R")
    trait1StateLabels <- getLabels(trait1)
    trait2StateLabels <- getLabels(trait2)
    
    # Assign labels
    dfCombined$trait1StateLabel <- NA
    dfCombined$trait2StateLabel <- NA
    
    # dfCombined$trait1StateLabel[dfCombined$ObservedState %in% c("ObsProp0Absent", "ObsProp0Present")] <- trait1StateLabels[1]
    # dfCombined$trait1StateLabel[dfCombined$ObservedState %in% c("ObsProp1Absent", "ObsProp1Present")] <- trait1StateLabels[2]
    # dfCombined$trait2StateLabel[dfCombined$ObservedState %in% c("ObsProp0Absent", "ObsProp1Absent")] <- trait2StateLabels[1]
    # dfCombined$trait2StateLabel[dfCombined$ObservedState %in% c("ObsProp0Present", "ObsProp1Present")] <- trait2StateLabels[2]
    
    dfCombined$trait1StateLabel[dfCombined$ObservedState %in% c("Overlap_0_0", "Overlap_0_1")] <- trait1StateLabels[1]
    dfCombined$trait1StateLabel[dfCombined$ObservedState %in% c("Overlap_1_0", "Overlap_1_1")] <- trait1StateLabels[2]
    dfCombined$trait2StateLabel[dfCombined$ObservedState %in% c("Overlap_0_0", "Overlap_1_0")] <- trait2StateLabels[1]
    dfCombined$trait2StateLabel[dfCombined$ObservedState %in% c("Overlap_0_1", "Overlap_1_1")] <- trait2StateLabels[2]
    
    dfCombined <- dfCombined %>%
      mutate(Label = paste(trait1StateLabel, trait2StateLabel, sep = "\n"))
    
    # Order levels
    levels_ordered <- c(
      paste(trait1StateLabels[1], trait2StateLabels[1], sep = "\n"),
      paste(trait1StateLabels[1], trait2StateLabels[2], sep = "\n"),
      paste(trait1StateLabels[2], trait2StateLabels[1], sep = "\n"),
      paste(trait1StateLabels[2], trait2StateLabels[2], sep = "\n")
    )
    dfCombined$Label <- factor(dfCombined$Label, levels = levels_ordered, ordered = TRUE)
    
    # Create the boxplot
    p3 <- ggplot(dfCombined, aes(x = Label, y = ObservedState.prop, fill = Which)) +
      geom_boxplot(outlier.shape = NA) +
      theme_minimal() +
      labs(y = "Observed State Proportion", x = "", fill = "Simulation Data") +
      scale_fill_manual(values = c("Real" = "#762a83", "Dummy" = "#1b7837")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ggtitle(paste(trait1, trait2, "p =", round(pval, 3)))
  } else {
    # Create a placeholder if columns don't exist
    p3 <- ggplot() + 
      theme_void() + 
      ggtitle("Boxplot requires standard column names")
  }
  
  # Save plots if requested
  if (plot_ggplots_pdf) {
    library(gridExtra)
    pdf(PDFpath, width = 8, height = 14)
    grid.arrange(p1, p2, p3, ncol = 1, heights = c(1, 1, 1.2))
    dev.off()
    cat("Plots saved to:", PDFpath, "\n")
  }
  
  # Calculate median values for all metrics
  if ("prop_trait1_state0" %in% names(dfout)) {
    # Calculate medians for real data
    medians_real <- data.frame(
      prop_trait1_state0_MedianReal = median(as.numeric(dfout$prop_trait1_state0)),
      prop_trait1_state1_MedianReal = median(as.numeric(dfout$prop_trait1_state1)),
      prop_trait2_state0_MedianReal = median(as.numeric(dfout$prop_trait2_state0)),
      prop_trait2_state1_MedianReal = median(as.numeric(dfout$prop_trait2_state1)),
      Overlap_0_0_MedianReal = median_real_0_0,
      Overlap_0_1_MedianReal = median_real_0_1,
      Overlap_1_0_MedianReal = median_real_1_0,
      Overlap_1_1_MedianReal = median_real_1_1
    )
    
    # Calculate medians for dummy data
    medians_dummy <- data.frame(
      prop_trait1_state0_MedianDummy = median(as.numeric(dfDummy$prop_trait1_state0)),
      prop_trait1_state1_MedianDummy = median(as.numeric(dfDummy$prop_trait1_state1)),
      prop_trait2_state0_MedianDummy = median(as.numeric(dfDummy$prop_trait2_state0)),
      prop_trait2_state1_MedianDummy = median(as.numeric(dfDummy$prop_trait2_state1)),
      Overlap_0_0_MedianDummy = median(as.numeric(dfDummy$Overlap_0_0)),
      Overlap_0_1_MedianDummy = median(as.numeric(dfDummy$Overlap_0_1)),
      Overlap_1_0_MedianDummy = median(as.numeric(dfDummy$Overlap_1_0)),
      Overlap_1_1_MedianDummy = median(as.numeric(dfDummy$Overlap_1_1))
    )
  } else {
    # For old format compatibility
    medians_real <- data.frame(
      propNoncoop_MedianReal = median(as.numeric(dfout$propNoncoop)),
      propCoop_MedianReal = median(as.numeric(dfout$propCoop)),
      propFSabsent_MedianReal = median(as.numeric(dfout$propFSabsent)),
      propFSpresent_MedianReal = median(as.numeric(dfout$propFSpresent)),
      ObsProp0Absent_MedianReal = median_real_0_0,
      ObsProp1Absent_MedianReal = median_real_0_1,
      ObsProp0Present_MedianReal = median_real_1_0,
      ObsProp1Present_MedianReal = median_real_1_1
    )
    
    medians_dummy <- data.frame(
      propNoncoop_MedianDummy = median(as.numeric(dfDummy$propNoncoop)),
      propCoop_MedianDummy = median(as.numeric(dfDummy$propCoop)),
      propFSabsent_MedianDummy = median(as.numeric(dfDummy$propFSabsent)),
      propFSpresent_MedianDummy = median(as.numeric(dfDummy$propFSpresent)),
      ObsProp0Absent_MedianDummy = median(as.numeric(dfDummy$ObsProp0Absent)),
      ObsProp1Absent_MedianDummy = median(as.numeric(dfDummy$ObsProp1Absent)),
      ObsProp0Present_MedianDummy = median(as.numeric(dfDummy$ObsProp0Present)),
      ObsProp1Present_MedianDummy = median(as.numeric(dfDummy$ObsProp1Present))
    )
  }
  
  # Perform pairwise comparisons if we have the combined data
  if (exists("dfCombined") && nrow(dfCombined) > 0) {
    # Linear model for pairwise comparisons
    lmStates <- lm(ObservedState.prop ~ Which * ObservedState, data = dfCombined)
    
    # Load emmeans package for pairwise comparisons
    if (requireNamespace("emmeans", quietly = TRUE)) {
      emm <- emmeans::emmeans(lmStates, pairwise ~ Which | ObservedState)
      contrast <- emm$contrasts
      PairwisePostHoc <- summary(contrast, adjust = "tukey")
      
      # Extract p-values and format them
      # The PairwisePostHoc data frame contains one row per state comparison
      # We need to extract the correct p-values
      pvals <- PairwisePostHoc$p.value
      states <- unique(PairwisePostHoc$ObservedState)
      
      # Initialize with NAs
      PairwisePvals <- data.frame(
        Overlap_0_0_DummyVsRealPval = NA,
        Overlap_0_1_DummyVsRealPval = NA,
        Overlap_1_0_DummyVsRealPval = NA,
        Overlap_1_1_DummyVsRealPval = NA
      )
      
      # Fill in the p-values based on what states are present
      for (i in seq_along(states)) {
        state <- as.character(states[i])
        if (state == "Overlap_0_0") PairwisePvals$Overlap_0_0_DummyVsRealPval <- pvals[i]
        else if (state == "Overlap_0_1") PairwisePvals$Overlap_0_1_DummyVsRealPval <- pvals[i]
        else if (state == "Overlap_1_0") PairwisePvals$Overlap_1_0_DummyVsRealPval <- pvals[i]
        else if (state == "Overlap_1_1") PairwisePvals$Overlap_1_1_DummyVsRealPval <- pvals[i]
      }
    } else {
      # If emmeans not available, create NA placeholders
      PairwisePvals <- data.frame(
        Overlap_0_0_DummyVsRealPval = NA,
        Overlap_0_1_DummyVsRealPval = NA,
        Overlap_1_0_DummyVsRealPval = NA,
        Overlap_1_1_DummyVsRealPval = NA
      )
    }
  } else {
    # Create NA placeholders if dfCombined doesn't exist
    PairwisePvals <- data.frame(
      Overlap_0_0_DummyVsRealPval = NA,
      Overlap_0_1_DummyVsRealPval = NA,
      Overlap_1_0_DummyVsRealPval = NA,
      Overlap_1_1_DummyVsRealPval = NA
    )
  }
  
  # Create summary statistics
  mediansRow <- data.frame(
    trait1 = trait1,
    trait2 = trait2,
    Nspecies = Nspecies,
    nsims_real = nsims_real,
    nsims_dummy = nsims_dummy,
    pval = pval,
    D_real = D_real,
    D_real_mean = mean(Real_dsims),
    D_real_sd = sd(Real_dsims),
    D_dummy_mean = mean(Dummy_dsums),
    D_dummy_sd = sd(Dummy_dsums)
  )
  
  # Add all the additional statistics
  mediansRow <- cbind(mediansRow, medians_real, medians_dummy, PairwisePvals)
  
  # Add fraction dummy less than median real
  mediansRow$Overlap_0_0_FractionDummyLessThanMedianReal <- Overlap_0_0_FractionDummyLessThanMedianReal
  mediansRow$Overlap_0_1_FractionDummyLessThanMedianReal <- Overlap_0_1_FractionDummyLessThanMedianReal
  mediansRow$Overlap_1_0_FractionDummyLessThanMedianReal <- Overlap_1_0_FractionDummyLessThanMedianReal
  mediansRow$Overlap_1_1_FractionDummyLessThanMedianReal <- Overlap_1_1_FractionDummyLessThanMedianReal
  
  # Save mediansRow as CSV
  csv_filename <- getStandardizedPath(dirs, "summaries",
                                     trait1, trait2, 
                                     file_type = "analysis_summary",
                                     nsims_real = nsims_real,
                                     nsims_dummy = nsims_dummy,
                                     other_label = otherlabel,
                                     extension = "csv")
  write.csv(mediansRow, file = csv_filename, row.names = FALSE)
  cat("Summary statistics saved to:", csv_filename, "\n")
  
  # Calculate transition counts if requested
  TransitionStats <- NULL
  p4 <- p5 <- p6 <- p4_gray <- NULL
  
  if (calculate_transitions) {
    if (is.null(trait1_simmaps) || is.null(trait2_simmaps)) {
      warning("Transition counts requested but simmaps not provided. Skipping transition analysis.")
    } else {
      # Source required functions
      # Try to use optimized version if available
      if (file.exists("getTransitionStateCounts_optimized.R")) {
        source("getTransitionStateCounts_optimized.R")
        use_optimized <- TRUE
      } else {
        source("getTransitionStateCounts_generic.R")
        use_optimized <- FALSE
      }
      source("transition_count_helpers.R")
      
      cat("Calculating transition counts...\n")
      
      # Calculate transition counts using appropriate function
      if (use_optimized) {
        transStateCounts <- getTransitionStateCounts_optimized(
          trait1_simtrees = trait1_simmaps,
          trait2_simtrees = trait2_simmaps,
          trait1_name = trait1,
          trait2_name = trait2,
          show_progress = TRUE
        )
      } else {
        transStateCounts <- getTransitionStateCounts_generic(
          trait1_simtrees = trait1_simmaps,
          trait2_simtrees = trait2_simmaps,
          trait1_name = trait1,
          trait2_name = trait2
        )
      }
      
      # Validate the counts
      if (!validateTransitionCounts(transStateCounts)) {
        warning("Transition count validation failed. Check the results carefully.")
      }
      
      # Calculate chi-square statistics
      transitioncols <- grep("^trait[12]_[01]to[01]_in_trait[12]_[01]$", names(transStateCounts), value = TRUE)
      debugPrint("Transition columns found in transStateCounts:")
      debugPrint("Number of transition columns:", length(transitioncols))
      debugPrint("Columns:", paste(transitioncols, collapse=", "))
      ExpectedCols <- paste0(transitioncols, "Expected")
      
      ChiSqStat <- rowSums((transStateCounts[,transitioncols] - transStateCounts[,ExpectedCols])^2 / 
                          transStateCounts[,ExpectedCols])
      ChiSqPvals <- pchisq(ChiSqStat, df = 7)
      transStateCounts <- cbind(transStateCounts, ChiSqStat, ChiSqPvals)
      
      # Process transition statistics
      require(emmeans)
      
      # Melt the transition counts for analysis
      transition_cols_all <- c(transitioncols, ExpectedCols)
      dfMeltCounts <- transStateCounts %>%
        gather(key = "Transition", value = "Count", all_of(transition_cols_all))
      
      # Identify observed vs expected
      dfMeltCounts$ObservedVsExpected <- ifelse(grepl("Expected", dfMeltCounts$Transition), 
                                               "Expected", "Observed")
      
      # Clean transition names
      dfMeltCounts$Transition <- gsub("Expected", "", dfMeltCounts$Transition)
      
      # Convert counts and calculate log
      dfMeltCounts$Count <- as.numeric(dfMeltCounts$Count)
      dfMeltCounts$Count[dfMeltCounts$Count == 0] <- 0.001
      dfMeltCounts$logCount <- log(dfMeltCounts$Count)
      
      # Linear model on log-transformed counts
      lmLogMult <- lm(logCount ~ ObservedVsExpected * Transition, data = dfMeltCounts)
      logLmANOVA <- anova(lmLogMult)
      
      # Pairwise post-hoc tests
      emmLog <- emmeans(lmLogMult, pairwise ~ ObservedVsExpected | Transition)
      contrastLog <- emmLog$contrasts
      logPairwisePostHoc <- summary(contrastLog, adjust = "tukey")
      
      debugPrint("emmeans analysis results:")
      debugPrint("Number of rows in logPairwisePostHoc:", nrow(logPairwisePostHoc))
      if (nrow(logPairwisePostHoc) > 0) {
        debugPrint("Transition column exists:", "Transition" %in% names(logPairwisePostHoc))
        if ("Transition" %in% names(logPairwisePostHoc)) {
          debugPrint("Unique transitions:", paste(unique(logPairwisePostHoc$Transition), collapse=", "))
        }
      }
      
      # Count times actual > expected
      observed_keys <- transitioncols
      expected_keys <- ExpectedCols
      timesActualGreaterThanExpected <- transStateCounts[, observed_keys] > transStateCounts[, expected_keys]
      NtimesActualGreaterThanExpected <- colSums(timesActualGreaterThanExpected, na.rm = TRUE)
      
      # Create TransitionStats object
      TransitionStats <- list(
        transStateCounts = transStateCounts,
        logLmANOVA = logLmANOVA,
        logPairwisePostHoc = logPairwisePostHoc,
        NtimesActualGreaterThanExpected = NtimesActualGreaterThanExpected
      )
      
      # Create transition plots if requested
      if (plot_transitions) {
        # Use the fixed version if it exists, otherwise use original
        if (file.exists("transition_plot_fixed.R")) {
          source("transition_plot_fixed.R")
        } else {
          source("transition_plot.R")
        }
        
        # Merge transition counts with overlap data for plotting
        # The transition plot needs both transition counts and overlap proportions
        df_for_plot <- merge(dfout, transStateCounts, by.x = "treenum", by.y = "TreeNum")
        
        # Check if we have transition columns after merge
        # Also check for columns with .x suffix (from merge with duplicate names)
        trans_cols_in_merged <- grep("^trait[12]_[01]to[01]_in_trait[12]_[01](\\.x)?$", names(df_for_plot), value = TRUE)
        
        # Get trait labels
        trait1StateLabels <- getLabels(trait1)
        trait2StateLabels <- getLabels(trait2)
        
        # Create RateRef dataframe
        RateRef <- createRateRef(trait1, trait2)
        
        # Prepare ratePvals dataframe
        transitions_for_pvals <- logPairwisePostHoc$Transition
        pvals_for_transitions <- logPairwisePostHoc$p.value
        
        pvaldf <- data.frame(
          Transition = transitions_for_pvals,
          p.value = pvals_for_transitions,
          stringsAsFactors = FALSE
        )
        
        # Add significance labels
        pvaldf$SignificanceLabel <- "n.s."
        pvaldf$SignificanceLabel[pvaldf$p.value < 0.05] <- "p < 0.05"
        pvaldf$SignificanceLabel[pvaldf$p.value < 0.01] <- "p < 0.01"
        pvaldf$SignificanceLabel[pvaldf$p.value < 0.001] <- "p < 0.001"
        pvaldf$SignificanceLabel[pvaldf$p.value < 0.0001] <- "p < 0.0001"
        
        # Add fraction actual > expected
        fractions <- NtimesActualGreaterThanExpected / nsims_real
        fractiondf <- data.frame(
          Transition = names(fractions),
          FractionActualGreaterThanExpected = as.numeric(fractions),
          stringsAsFactors = FALSE
        )
        
        pvaldf <- merge(pvaldf, fractiondf, by = "Transition")
        
        # Add PercentTrendingLabel
        pvaldf$PercentTrendingLabel <- "<80%"
        pvaldf$PercentTrendingLabel[pvaldf$FractionActualGreaterThanExpected > .8 | 
                                    pvaldf$FractionActualGreaterThanExpected < .2] <- ">80%"
        pvaldf$PercentTrendingLabel[pvaldf$FractionActualGreaterThanExpected > .9 | 
                                    pvaldf$FractionActualGreaterThanExpected < .1] <- ">90%"
        pvaldf$PercentTrendingLabel[pvaldf$FractionActualGreaterThanExpected > .95 | 
                                    pvaldf$FractionActualGreaterThanExpected < .05] <- ">95%"
        pvaldf$PercentTrendingLabel[pvaldf$FractionActualGreaterThanExpected > .99 | 
                                    pvaldf$FractionActualGreaterThanExpected < .01] <- ">99%"
        
        # Map to arrow codes
        pvaldf$qRate <- mapGenericToArrowNames(pvaldf$Transition)
        
        ratePvals <- pvaldf[!is.na(pvaldf$qRate), ]
        
        names(ratePvals)[names(ratePvals) == "Transition"] <- "Transitions"
        
        # Calculate TransitionCountMedian BEFORE prepareTransitionDataForPlot modifies the columns
        if (length(trans_cols_in_merged) > 0) {
          # Map the generic transition names to the actual column names (which might have .x suffix)
          actual_col_mapping <- setNames(trans_cols_in_merged, gsub("\\.x$", "", trans_cols_in_merged))
          
          # Get the actual column names for the transitions we need
          cols_to_use <- actual_col_mapping[ratePvals$Transitions]
          
          if (any(!is.na(cols_to_use))) {
            # Ensure columns are numeric before calculating medians
            cols_valid <- cols_to_use[!is.na(cols_to_use)]
            data_subset <- df_for_plot[, cols_valid, drop = FALSE]
            
            # Convert to numeric if needed
            for (i in seq_along(data_subset)) {
              if (!is.numeric(data_subset[[i]])) {
                data_subset[[i]] <- as.numeric(data_subset[[i]])
              }
            }
            
            # Calculate medians from the transition count columns
            trans_medians <- apply(data_subset, 2, median, na.rm = TRUE)
            # Map back to the original transition names
            names(trans_medians) <- names(cols_valid)
            ratePvals$TransitionCountMedian <- trans_medians[ratePvals$Transitions]
          } else {
            warning("Could not map transition names to actual columns")
            ratePvals$TransitionCountMedian <- 5
          }
        } else {
          warning("No transition columns found in df_for_plot, using default TransitionCountMedian")
          ratePvals$TransitionCountMedian <- 5
        }
        
        # Prepare data for arrow plots (this converts transition columns to q columns)
        df_for_plot <- prepareTransitionDataForPlot(df_for_plot, trait1, trait2)
        
        # Create transition plot
        plottitle <- paste(trait1, trait2, "nsims:", nsims_real)
        
        # Check if we have the required data
        cat("\nChecking data for arrow plot:\n")
        cat("- df_for_plot rows:", nrow(df_for_plot), "\n")
        cat("- ratePvals rows:", nrow(ratePvals), "\n")
        cat("- Columns in df_for_plot:", paste(names(df_for_plot)[1:10], collapse=", "), "...\n")
        
        transition_plot_result <- try(transition_plot(
          df = df_for_plot,
          trait1StateLabels = trait1StateLabels,
          trait2StateLabels = trait2StateLabels,
          scale_area_by = 1,
          offset = 0.2,
          lengthen = 0.2,
          ratePvals = ratePvals,
          plottitle = plottitle,
          center = "median"
        ), silent = TRUE)
        
        if (inherits(transition_plot_result, "try-error")) {
          cat("\nError creating arrow plot:", as.character(transition_plot_result), "\n")
          p4 <- ggplot() + 
            ggtitle(paste("Arrow plot error for", trait1, "vs", trait2)) +
            theme_void()
          p4_gray <- NULL
        } else {
          # Extract the plots from the result
          p4 <- transition_plot_result$transition_plot
          p4_gray <- transition_plot_result$transitionplot_GrayNS
        }
        
        # Create transition count boxplots
        require(ggpubr)
        
        # Get interaction p-values for labels
        pvalInteraction <- logLmANOVA$`Pr(>F)`[3]
        pvalInteractionLabel <- ifelse(pvalInteraction < 0.0001, 
                                      "p < 0.0001", 
                                      paste("p =", round(pvalInteraction, digits = 5)))
        
        # Create group labels with N times actual > expected
        transition_levels <- unique(dfMeltCounts$Transition)
        groupLabels <- paste0(transition_levels, 
                            "\nNsims Actual >\nExpected: ", 
                            NtimesActualGreaterThanExpected[transition_levels], 
                            "/", nsims_real)
        names(groupLabels) <- transition_levels
        
        # Create transition count boxplot (non-log)
        p5 <- ggplot(dfMeltCounts, aes(x = Transition, y = Count, fill = ObservedVsExpected)) +
          geom_boxplot(outlier.shape = NA) +
          theme_minimal() +
          labs(y = "Transition Counts", x = "", fill = "Observed/Expected") +
          scale_fill_manual(values = c("Observed" = "#762a83", "Expected" = "#1b7837")) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1), 
                title = element_text(size = 8)) +
          ggtitle(paste(trait1, trait2, "\nnSimsObserved =", nsims_real, 
                       "    nSimsExpected =", nsims_real, 
                       "    Obs/Exp:TransCounts", pvalInteractionLabel)) +
          stat_compare_means(aes(label = after_stat(p.signif)), method = "t.test")
        
        # Create log transition count boxplot
        p6 <- ggplot(dfMeltCounts, aes(x = Transition, y = logCount, fill = ObservedVsExpected)) +
          geom_boxplot(outlier.shape = NA) +
          theme_minimal() +
          labs(y = "log(Transition Counts)", x = "", fill = "Observed/Expected") +
          scale_fill_manual(values = c("Observed" = "#762a83", "Expected" = "#1b7837")) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6), 
                title = element_text(size = 8)) +
          ggtitle(paste(trait1, trait2, "\nnSimsObserved =", nsims_real, 
                       "    nSimsExpected =", nsims_real, 
                       "    Obs/Exp:TransCounts", pvalInteractionLabel)) +
          stat_compare_means(aes(label = after_stat(p.signif)), method = "t.test") +
          scale_x_discrete(labels = groupLabels)
        
        # Extract arrow plots
        if (!inherits(transition_plot_result, "try-error") && is.list(transition_plot_result)) {
          p4 <- transition_plot_result$transition_plot
          # Debug info
          debugPrint("Arrow plot (p4) debug info:")
          debugPrint("- Class of p4:", class(p4))
          debugPrint("- Is NULL:", is.null(p4))
          if (inherits(p4, "ggplot")) {
            debugPrint("- Number of layers:", length(p4$layers))
            debugPrint("- Has data:", !is.null(p4$data))
          }
          # Could use gray NS plot as additional plot if needed
        } else {
          cat("\nArrow plot creation failed or returned non-list\n")
          p4 <- NULL
        }
        
        cat("Transition plots created\n")
      }
      
      # Save transition counts
      transition_csv <- getStandardizedPath(dirs, "data",
                                           trait1, trait2,
                                           file_type = "transition_counts",
                                           nsims_real = nsims_real,
                                           other_label = otherlabel,
                                           extension = "csv")
      write.csv(transStateCounts, file = transition_csv, row.names = FALSE)
      cat("Transition counts saved to:", transition_csv, "\n")
    }
  }
  
  # Return results
  result <- list(
    Real_dsims = Real_dsims,
    Dummy_dsums = Dummy_dsums,
    D_real = D_real,
    pval = pval,
    mediansRow = mediansRow
  )
  
  # Add plots
  if (exists("p3")) {
    result$plots <- list(p1 = p1, p2 = p2, p3 = p3)
  } else {
    result$plots <- list(p1 = p1, p2 = p2)
  }
  
  # Add transition results if calculated
  if (!is.null(TransitionStats)) {
    result$TransitionStats <- TransitionStats
    if (!is.null(p4)) result$plots$p4 <- p4
    if (!is.null(p4_gray)) result$plots$p4_gray <- p4_gray
    if (!is.null(p5)) result$plots$p5 <- p5
    if (!is.null(p6)) result$plots$p6 <- p6
  }
  
  return(result)
}
