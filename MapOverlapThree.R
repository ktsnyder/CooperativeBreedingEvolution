# Adapted from Map.Overlap from package 'phytools' version 1.9-16
# Using Claude 3.5 Sonnet
# By Kate Snyder
# 2024-12-11

Map.Overlap.Three <- function(tree1, tree2, tree3, tol = 1e-06, standardize = TRUE, ...) {
  # Check arguments
  if (hasArg(check.equal)) 
    check.equal <- list(...)$check.equal
  else check.equal <- TRUE
  
  # Verify tree classes
  if (!inherits(tree1, "phylo") || !inherits(tree2, "phylo") || !inherits(tree3, "phylo"))
    stop("All input trees should be objects of class \"phylo\".")
  
  # Check tree structures match if required
  if (check.equal) {
    if (!all.equal.phylo(tree1, tree2, tolerance = tol) || 
        !all.equal.phylo(tree2, tree3, tolerance = tol))
      stop("Mapped trees must have the same underlying structure.")
  }
  
  if (!inherits(tree1, "simmap") || !inherits(tree2, "simmap") || !inherits(tree3, "simmap"))
    stop("All input trees should be objects of class \"simmap\".")
  
  # Get unique states from each tree
  s1 <- mapped.states(tree1)
  s2 <- mapped.states(tree2)
  s3 <- mapped.states(tree3)
  
  # Create 3D array for results
  R <- array(0, 
             dim = c(length(s1), length(s2), length(s3)),
             dimnames = list(s1, s2, s3))
  
  # Iterate through edges
  for (i in 1:nrow(tree1$edge)) {
    # Create matrices for start/end points for each tree
    XX <- matrix(0, length(tree1$maps[[i]]), 2, 
                 dimnames = list(names(tree1$maps[[i]]), c("start", "end")))
    YY <- matrix(0, length(tree2$maps[[i]]), 2,
                 dimnames = list(names(tree2$maps[[i]]), c("start", "end")))
    ZZ <- matrix(0, length(tree3$maps[[i]]), 2,
                 dimnames = list(names(tree3$maps[[i]]), c("start", "end")))
    
    # Fill first row
    XX[1, 2] <- tree1$maps[[i]][1]
    YY[1, 2] <- tree2$maps[[i]][1]
    ZZ[1, 2] <- tree3$maps[[i]][1]
    
    # Fill remaining rows for tree1
    if (length(tree1$maps[[i]]) > 1) {
      for (j in 2:length(tree1$maps[[i]])) {
        XX[j, 1] <- XX[j - 1, 2]
        XX[j, 2] <- XX[j, 1] + tree1$maps[[i]][j]
      }
    }
    
    # Fill remaining rows for tree2
    if (length(tree2$maps[[i]]) > 1) {
      for (j in 2:length(tree2$maps[[i]])) {
        YY[j, 1] <- YY[j - 1, 2]
        YY[j, 2] <- YY[j, 1] + tree2$maps[[i]][j]
      }
    }
    
    # Fill remaining rows for tree3
    if (length(tree3$maps[[i]]) > 1) {
      for (j in 2:length(tree3$maps[[i]])) {
        ZZ[j, 1] <- ZZ[j - 1, 2]
        ZZ[j, 2] <- ZZ[j, 1] + tree3$maps[[i]][j]
      }
    }
    
    # Calculate overlaps
    for (j in 1:nrow(XX)) {
      lower_y <- max(which(YY[, 1] <= XX[j, 1]))
      upper_y <- which(YY[, 2] >= (XX[j, 2] - tol))[1]
      
      for (k in lower_y:upper_y) {
        # Find overlapping region between XX and YY
        overlap_start <- max(YY[k, 1], XX[j, 1])
        overlap_end <- min(YY[k, 2], XX[j, 2])
        
        # Find which segments in ZZ overlap with this region
        lower_z <- max(which(ZZ[, 1] <= overlap_start))
        upper_z <- which(ZZ[, 2] >= (overlap_end - tol))[1]
        
        for (l in lower_z:upper_z) {
          overlap_length <- min(min(ZZ[l, 2], overlap_end) - 
                                  max(ZZ[l, 1], overlap_start))
          
          if (overlap_length > tol) {
            R[rownames(XX)[j], 
              rownames(YY)[k],
              rownames(ZZ)[l]] <- R[rownames(XX)[j],
                                    rownames(YY)[k],
                                    rownames(ZZ)[l]] + overlap_length
          }
        }
      }
    }
  }
  
  if (standardize) 
    R/sum(R)
  else R
}


# First, let's create a function to handle this conversion
convert_overlap_3d_to_vector <- function(overlap_array, suffix) {
  # Convert 3D array to vector
  overlapVec <- as.vector(overlap_array)
  
  # Create names for each combination
  dim_names <- dimnames(overlap_array)
  
  # Use expand.grid to get all combinations of states
  name_combinations <- expand.grid(
    dim_names[[1]],  # states from first tree
    dim_names[[2]],  # states from second tree
    dim_names[[3]]   # states from third tree
  )
  
  # Create new names by combining states with underscores
  new_names <- paste(
    name_combinations$Var1,
    name_combinations$Var2,
    name_combinations$Var3,
    sep = "_"
  )
  
  # Add the "_REAL or _DUMMY" suffix
  names(overlapVec) <- paste0(new_names, suffix)
  
  return(overlapVec)
}


calcHuelflex_three = function(overlapdf) {
  require(dplyr)
  require(tidyr)
  require(stringr)
  require(ggplot2)
  
  if("X" %in% colnames(overlapdf)) {
    overlapdf$X = NULL
  }
  
  nsims = length(overlapdf[,1])
  numericColumnStart = max(which(str_detect(colnames(overlapdf), "trait"))) + 1
  
  # Convert numeric columns
  overlapdf[,numericColumnStart:length(colnames(overlapdf))] = 
    apply(overlapdf[,numericColumnStart:length(colnames(overlapdf))], MARGIN = 2, FUN = as.numeric)
  
  # Modified pivot_longer to handle three traits
  overlaplonger = overlapdf %>%
    pivot_longer(
      cols = !c(tree, trait1, trait2, trait3),
      names_to = "state",
      values_to = "proportion"
    )
  
  overlaplonger$Which = NA
  overlaplonger$Which[which(str_detect(overlaplonger$state, "DUMMY"))] = "Dummy"
  overlaplonger$Which[which(str_detect(overlaplonger$state, "REAL"))] = "Real"
  
  overlaplonger$state <- gsub("_DUMMY", "", overlaplonger$state)
  overlaplonger$state <- gsub("_REAL", "", overlaplonger$state)
  
  # Modified to handle three-way state combinations
  df_long <- overlaplonger %>%
    separate(state, into = c("trait_state", "second_state", "third_state"), sep = "_") %>%
    mutate(
      FS = paste(second_state, third_state, sep = "_")
    )
  
  total_times_trait <- df_long %>%
    group_by(tree, trait1, trait2, trait3, trait_state, Which) %>%
    summarize(total_trait = sum(proportion), .groups = "drop")
  
  total_times_FS <- df_long %>%
    group_by(tree, trait1, trait2, trait3, FS, Which) %>%
    summarize(total_FS = sum(proportion), .groups = "drop")
  
  df_long <- df_long %>%
    left_join(total_times_trait, by = c("tree", "trait1", "trait2", "trait3", "trait_state", "Which")) %>%
    left_join(total_times_FS, by = c("tree", "trait1", "trait2", "trait3", "FS", "Which"))
  
  df_long <- df_long %>%
    mutate(expected_proportion = total_trait * total_FS)
  
  df_real <- df_long %>%
    filter(Which == "Real")
  
  df_real <- df_real %>%
    mutate(abs_diff = abs(proportion - expected_proportion))
  
  D_real <- df_real %>%
    summarize(total_abs_diff = sum(abs_diff)) %>%
    pull(total_abs_diff) / nsims
  
  Real_dsims <- df_real %>%
    select(tree, trait1, trait2, trait3, trait_state, FS, abs_diff) %>%
    group_by(tree, trait1, trait2, trait3) %>%
    summarize(Real_dsim = sum(abs_diff), .groups = "drop")
  
  # Dummy data analysis
  df_dummy <- df_long %>%
    filter(Which == "Dummy")
  
  df_dummy <- df_dummy %>%
    mutate(abs_diff = abs(proportion - expected_proportion))
  
  Dummy_dsums <- df_dummy %>%
    group_by(tree, trait1, trait2, trait3) %>%
    summarize(Dummy_dsum = sum(abs_diff), .groups = "drop")
  
  numGreater = sum(Dummy_dsums$Dummy_dsum > D_real)
  pval = numGreater/nsims
  
  # Calculate medians for Real data
  medians_real <- df_long %>%
    filter(Which == "Real") %>%
    group_by(trait_state, FS) %>%
    summarize(median_real = median(proportion), .groups = "drop")
  
  # Join medians to Dummy data
  df_dummy <- df_dummy %>%
    left_join(medians_real, by = c("trait_state", "FS"))
  
  fraction_dummy_less_than_median_real <- df_dummy %>%
    group_by(trait_state, FS) %>%
    summarize(fraction = sum(proportion <= median_real) / n(), .groups = "drop")
  
  trait1 = overlaplonger$trait1[1]
  trait2 = overlaplonger$trait2[1]
  trait3 = overlaplonger$trait3[1]
  
  plotlabel = paste(trait1, trait2, trait3)
  dummytitle = paste("Nsims =", nsims, "\nnum Dummy dsums > D_real:", numGreater, ", pval =", pval)
  
  xmax = max(c(Real_dsims$Real_dsim, Dummy_dsums$Dummy_dsum))*1.1
  
  # Plotting
  p1 <- ggplot(Real_dsims, aes(x = Real_dsim)) +
    geom_histogram(binwidth = xmax/20, fill = rgb(0.2, 0.5, 0.7, 0.5), color = "white") +
    geom_vline(xintercept = D_real, color = "red") +
    xlim(c(0, xmax)) +
    labs(title = plotlabel, x = "D statistic from real data simmaps", y = "Frequency") +
    theme_minimal(base_size = 10)
  
  p2 <- ggplot(Dummy_dsums, aes(x = Dummy_dsum)) +
    geom_histogram(binwidth = xmax/20, fill = rgb(0.7, 0.5, 0.2, 0.5), color = "white") +
    xlim(c(0, xmax)) +
    labs(title = dummytitle, x = "D statistic from simulated independent data simmaps", y = "Frequency") +
    theme_minimal(base_size = 10)
  
  boxplotStates <- ggplot(overlaplonger, aes(x = state, y = proportion, fill = Which)) +
    geom_boxplot(outlier.shape = NA) +
    theme_minimal() +
    labs(y = "Observed State Proportion", x = "", fill = "Simulation Data") +
    scale_fill_manual(values = c("Real" = "blue", "Dummy" = "red")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(paste(trait1, trait2, trait3, "p =", pval))
  
  
  boxplotStatesLowestProportions <- ggplot(overlaplonger, aes(x = state, y = proportion, fill = Which)) +
    geom_boxplot(outlier.shape = NA) +
    theme_minimal() +
    labs(y = "Observed State Proportion - low range", x = "", fill = "Simulation Data") +
    scale_fill_manual(values = c("Real" = "blue", "Dummy" = "red")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    coord_cartesian(ylim = c(0, 0.10)) +
    ggtitle(paste(trait1, trait2, trait3, "p =", pval))
  
  # New faceted boxplot
  # First, extract the third trait state from the state column
  overlaplonger_split <- overlaplonger %>%
    mutate(
      third_trait_state = str_sub(state, -1, -1),  # Get last character of state
      other_states = substr(state, 1, nchar(state)-2)  # Get all but last two characters
    )
  
  boxplotStatesFaceted_third <- ggplot(overlaplonger_split, 
                                 aes(x = other_states, y = proportion, fill = Which)) +
    geom_boxplot(outlier.shape = NA) +
    facet_wrap(~third_trait_state, 
               labeller = labeller(third_trait_state = 
                                     function(x) paste(trait3, "=", x))) +
    theme_minimal() +
    labs(y = "Observed State Proportion", 
         x = "", 
         fill = "Simulation Data") +
    scale_fill_manual(values = c("Real" = "blue", "Dummy" = "red")) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 10),
      strip.background = element_rect(fill = "lightgray", color = NA)
    ) +
    ggtitle(paste(trait1, trait2, "split by", trait3))
  
  
  # New faceted plot by second trait
  overlaplonger_split_second <- overlaplonger %>%
    mutate(
      second_trait_state = str_sub(state, -3, -3),  # Get second-to-last number
      other_states = paste(substr(state, 1, nchar(state)-4),  # First part
                           str_sub(state, -1, -1), sep = "_")   # Last number
    )
  
  boxplotStatesFaceted_second <- ggplot(overlaplonger_split_second, 
                                        aes(x = other_states, y = proportion, fill = Which)) +
    geom_boxplot(outlier.shape = NA) +
    facet_wrap(~second_trait_state, 
               labeller = labeller(second_trait_state = 
                                     function(x) paste(trait2, "=", x))) +
    theme_minimal() +
    labs(y = "Observed State Proportion", 
         x = "", 
         fill = "Simulation Data") +
    scale_fill_manual(values = c("Real" = "blue", "Dummy" = "red")) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 10),
      strip.background = element_rect(fill = "lightgray", color = NA)
    ) +
    ggtitle(paste(trait1, trait3, "split by", trait2))
  
  return(list(
    p1 = p1, 
    p2 = p2, 
    boxplotStates = boxplotStates, 
    boxplotStatesYcutoff = boxplotStatesLowestProportions,
    boxplotStatesFaceted_second = boxplotStatesFaceted_second, 
    boxplotStatesFaceted_third = boxplotStatesFaceted_third, 
    fraction_dummy_less_than_median_real = fraction_dummy_less_than_median_real
  ))
}
