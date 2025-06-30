## Corrected calcHuel function with proper state indexing
## Kate Snyder
## 2025-06-26
## Fixes the set.seed bug and uses generic column names

calcHuel_corrected <- function(dfout, dfDummy, nsims_real = NULL, nsims_dummy = NULL, 
                              otherlabel = NULL, newplot = TRUE, plot_ggplots_pdf = FALSE,
                              trait1_name = "Trait1", trait2_name = "Trait2") {
  require(tidyverse)
  require(ggplot2)

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
  
  PDFname = paste("Simmap Overlap Counts CORRECTED", trait1, trait2, nsims_real, nsims_dummy, otherlabel)
  
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
    
    # Get trait labels - define getLabels function locally
    getLabels <- function(trait) {
      if (grepl("coop", trait, ignore.case = TRUE)) {
        return(c("Non-Cooperative", "Cooperative"))
      } else if (grepl("FemaleSong", trait)) {
        return(c("Female Song Absent", "Female Song Present"))  
      } else if (grepl("Territory_12", trait)) {
        return(c("Non-Territorial", "Year-round Territorial"))
      } else if (grepl("WeakVsStrong", trait)) {
        return(c("Weak or No Territoriality", "Strong Territoriality"))
      } else {
        return(c(paste(trait, "0"), paste(trait, "1")))
      }
    }
    
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
    pdf_filename <- paste0("Simmap Overlap Outputs/", PDFname, ".pdf")
    pdf(pdf_filename, width = 8, height = 14)
    grid.arrange(p1, p2, p3, ncol = 1, heights = c(1, 1, 1.2))
    dev.off()
    cat("Plots saved to:", pdf_filename, "\n")
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
  if (!dir.exists("Simmap Overlap Outputs")) {
    dir.create("Simmap Overlap Outputs")
  }
  csv_filename <- paste0("Simmap Overlap Outputs/", 
                         gsub(" ", "_", PDFname), "_summary.csv")
  write.csv(mediansRow, file = csv_filename, row.names = FALSE)
  cat("Summary statistics saved to:", csv_filename, "\n")
  
  # Return results
  if (exists("p3")) {
    return(list(
      Real_dsims = Real_dsims,
      Dummy_dsums = Dummy_dsums,
      D_real = D_real,
      pval = pval,
      mediansRow = mediansRow,
      plots = list(p1 = p1, p2 = p2, p3 = p3)
    ))
  } else {
    return(list(
      Real_dsims = Real_dsims,
      Dummy_dsums = Dummy_dsums,
      D_real = D_real,
      pval = pval,
      mediansRow = mediansRow,
      plots = list(p1 = p1, p2 = p2)
    ))
  }
}
