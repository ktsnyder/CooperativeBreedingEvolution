## Optimized generic transition state count function with progress tracking
## Based on getTransitionStateCounts_generic.R
## Includes performance improvements and progress indicators
## Kate Snyder with assistance from Claude (Anthropic)
## 2025-06-30

getTransitionStateCounts_optimized <- function(trait1_simtrees, trait2_simtrees, 
                                              trait1_name = "trait1", trait2_name = "trait2",
                                              show_progress = TRUE) {
  
  require(phytools)
  require(dplyr)
  
  nTrees = min(c(length(trait1_simtrees), length(trait2_simtrees)))
  
  if (show_progress) {
    cat("Calculating transition counts for", nTrees, "trees...\n")
    pb <- txtProgressBar(min = 0, max = nTrees, style = 3)
  }
  
  # Pre-allocate result matrix for efficiency
  result_matrix <- matrix(NA, nrow = nTrees, ncol = 12)
  colnames(result_matrix) <- c(
    "trait1_0to1_in_trait2_0", "trait1_1to0_in_trait2_0", 
    "trait1_0to1_in_trait2_1", "trait1_1to0_in_trait2_1", 
    "trait2_0to1_in_trait1_0", "trait2_1to0_in_trait1_0", 
    "trait2_0to1_in_trait1_1", "trait2_1to0_in_trait1_1", 
    "trait1_PropTime0", "trait1_PropTime1", 
    "trait2_PropTime0", "trait2_PropTime1"
  )
  
  # Track timing for performance monitoring
  timing_info <- list()
  start_time <- Sys.time()
  
  for (t in 1:nTrees) {
    tree_start <- Sys.time()
    
    trait1_simtree1 = trait1_simtrees[[t]]
    trait2_simtree1 = trait2_simtrees[[t]]
    
    # Time state proportions
    trait1_timePerState = describe.simmap(trait1_simtree1)
    trait1_PropTime0 = trait1_timePerState$times["prop","0"]
    trait1_PropTime1 = trait1_timePerState$times["prop","1"]
    
    trait2_timePerState = describe.simmap(trait2_simtree1)
    trait2_PropTime0 = trait2_timePerState$times["prop","0"]
    trait2_PropTime1 = trait2_timePerState$times["prop","1"]
    
    # Get segments using optimized function
    segment_start <- Sys.time()
    trait1_segments = getSimmapSegments_optimized(trait1_simtree1)
    trait2_segments = getSimmapSegments_optimized(trait2_simtree1)
    segment_time <- difftime(Sys.time(), segment_start, units = "secs")
    
    trait1_Transitions = trait1_segments$transitions_df
    trait1_Segments = trait1_segments$segments_df
    trait2_Transitions = trait2_segments$transitions_df
    trait2_Segments = trait2_segments$segments_df
    
    # Initialize counts
    counts <- numeric(8)
    names(counts) <- c(
      "trait1_0to1_in_trait2_0", "trait1_1to0_in_trait2_0",
      "trait1_0to1_in_trait2_1", "trait1_1to0_in_trait2_1",
      "trait2_0to1_in_trait1_0", "trait2_1to0_in_trait1_0",
      "trait2_0to1_in_trait1_1", "trait2_1to0_in_trait1_1"
    )
    
    # Process trait1 transitions with optimized lookup
    trans_start <- Sys.time()
    if (nrow(trait1_Transitions) > 0) {
      trait1_Transitions <- trait1_Transitions[!is.na(trait1_Transitions$Transition01or10), ]
      
      for (i in 1:nrow(trait1_Transitions)) {
        transition <- trait1_Transitions[i, ]
        
        # Optimized state lookup using vectorized operations
        trait2_state_at_transition <- find_other_trait_state_optimized(
          transition$EdgeNumber, 
          transition$DistanceFromRootwardNode, 
          trait2_Segments
        )
        
        if (!is.na(trait2_state_at_transition)) {
          # Build count key
          count_key <- paste0("trait1_", transition$Transition01or10, "_in_trait2_", trait2_state_at_transition)
          if (count_key %in% names(counts)) {
            counts[count_key] <- counts[count_key] + 1
          }
        }
      }
    }
    
    # Process trait2 transitions
    if (nrow(trait2_Transitions) > 0) {
      trait2_Transitions <- trait2_Transitions[!is.na(trait2_Transitions$Transition01or10), ]
      
      for (i in 1:nrow(trait2_Transitions)) {
        transition <- trait2_Transitions[i, ]
        
        trait1_state_at_transition <- find_other_trait_state_optimized(
          transition$EdgeNumber, 
          transition$DistanceFromRootwardNode, 
          trait1_Segments
        )
        
        if (!is.na(trait1_state_at_transition)) {
          count_key <- paste0("trait2_", transition$Transition01or10, "_in_trait1_", trait1_state_at_transition)
          if (count_key %in% names(counts)) {
            counts[count_key] <- counts[count_key] + 1
          }
        }
      }
    }
    trans_time <- difftime(Sys.time(), trans_start, units = "secs")
    
    # Store results in pre-allocated matrix
    result_matrix[t, ] <- c(counts, trait1_PropTime0, trait1_PropTime1, 
                           trait2_PropTime0, trait2_PropTime1)
    
    # Update progress
    if (show_progress) {
      setTxtProgressBar(pb, t)
    }
    
    # Store timing info for first few trees
    if (t <= 5) {
      timing_info[[t]] <- list(
        total = difftime(Sys.time(), tree_start, units = "secs"),
        segments = segment_time,
        transitions = trans_time
      )
    }
  }
  
  if (show_progress) {
    close(pb)
  }
  
  # Convert to dataframe
  countsdf <- as.data.frame(result_matrix)
  countsdf$TreeNum <- 1:nTrees
  
  # Calculate totals (these were missing from optimized version)
  countsdf$Total_trait1_0to1 <- countsdf$trait1_0to1_in_trait2_0 + countsdf$trait1_0to1_in_trait2_1
  countsdf$Total_trait1_1to0 <- countsdf$trait1_1to0_in_trait2_0 + countsdf$trait1_1to0_in_trait2_1
  countsdf$Total_trait2_0to1 <- countsdf$trait2_0to1_in_trait1_0 + countsdf$trait2_0to1_in_trait1_1
  countsdf$Total_trait2_1to0 <- countsdf$trait2_1to0_in_trait1_0 + countsdf$trait2_1to0_in_trait1_1
  
  # Calculate expected counts
  countsdf$trait1_0to1Expected <- (countsdf$trait1_0to1_in_trait2_0 + countsdf$trait1_0to1_in_trait2_1) * countsdf$trait2_PropTime0
  countsdf$trait1_1to0Expected <- (countsdf$trait1_1to0_in_trait2_0 + countsdf$trait1_1to0_in_trait2_1) * countsdf$trait2_PropTime0
  countsdf$trait2_0to1Expected <- (countsdf$trait2_0to1_in_trait1_0 + countsdf$trait2_0to1_in_trait1_1) * countsdf$trait1_PropTime0
  countsdf$trait2_1to0Expected <- (countsdf$trait2_1to0_in_trait1_0 + countsdf$trait2_1to0_in_trait1_1) * countsdf$trait1_PropTime0
  
  # Expected counts by state
  countsdf$trait1_0to1_in_trait2_0Expected <- countsdf$trait1_0to1Expected
  countsdf$trait1_0to1_in_trait2_1Expected <- (countsdf$trait1_0to1_in_trait2_0 + countsdf$trait1_0to1_in_trait2_1) * countsdf$trait2_PropTime1
  countsdf$trait1_1to0_in_trait2_0Expected <- countsdf$trait1_1to0Expected
  countsdf$trait1_1to0_in_trait2_1Expected <- (countsdf$trait1_1to0_in_trait2_0 + countsdf$trait1_1to0_in_trait2_1) * countsdf$trait2_PropTime1
  
  countsdf$trait2_0to1_in_trait1_0Expected <- countsdf$trait2_0to1Expected
  countsdf$trait2_0to1_in_trait1_1Expected <- (countsdf$trait2_0to1_in_trait1_0 + countsdf$trait2_0to1_in_trait1_1) * countsdf$trait1_PropTime1
  countsdf$trait2_1to0_in_trait1_0Expected <- countsdf$trait2_1to0Expected
  countsdf$trait2_1to0_in_trait1_1Expected <- (countsdf$trait2_1to0_in_trait1_0 + countsdf$trait2_1to0_in_trait1_1) * countsdf$trait1_PropTime1
  
  # Reorder columns to match original function
  col_order <- c("TreeNum", 
                 "trait1_0to1_in_trait2_0", "trait1_1to0_in_trait2_0", 
                 "trait1_0to1_in_trait2_1", "trait1_1to0_in_trait2_1", 
                 "trait2_0to1_in_trait1_0", "trait2_1to0_in_trait1_0", 
                 "trait2_0to1_in_trait1_1", "trait2_1to0_in_trait1_1", 
                 "trait1_PropTime0", "trait1_PropTime1", 
                 "trait2_PropTime0", "trait2_PropTime1",
                 "Total_trait1_0to1", "Total_trait1_1to0",
                 "Total_trait2_0to1", "Total_trait2_1to0",
                 grep("Expected", names(countsdf), value = TRUE))
  
  countsdf <- countsdf[, col_order]
  
  # Print performance summary
  total_time <- difftime(Sys.time(), start_time, units = "secs")
  if (show_progress) {
    cat("\nTransition count calculation complete.\n")
    cat("Total time:", round(total_time, 2), "seconds\n")
    cat("Average time per tree:", round(total_time / nTrees, 3), "seconds\n")
    
    if (length(timing_info) > 0) {
      avg_segment_time <- mean(sapply(timing_info, function(x) x$segments))
      avg_trans_time <- mean(sapply(timing_info, function(x) x$transitions))
      cat("\nTiming breakdown (first 5 trees):\n")
      cat("  Average segment extraction:", round(avg_segment_time, 3), "seconds\n")
      cat("  Average transition counting:", round(avg_trans_time, 3), "seconds\n")
    }
  }
  
  return(countsdf)
}

# Optimized version of getSimmapSegments using pre-allocation
getSimmapSegments_optimized <- function(tree) {
  outlist <- list()
  
  # Count total segments and transitions first
  n_transitions <- 0
  n_segments <- 0
  for (edge_data in tree$maps) {
    n_segments <- n_segments + length(edge_data)
    if (length(edge_data) > 1) {
      n_transitions <- n_transitions + length(edge_data) - 1
    }
  }
  
  # Pre-allocate data frames
  transitions_df <- data.frame(
    EdgeNumber = integer(n_transitions),
    Transition01or10 = character(n_transitions),
    DistanceFromRootwardNode = numeric(n_transitions),
    stringsAsFactors = FALSE
  )
  
  segments_df <- data.frame(
    EdgeNumber = integer(n_segments),
    State = numeric(n_segments),
    Start = numeric(n_segments),
    End = numeric(n_segments),
    stringsAsFactors = FALSE
  )
  
  trans_idx <- 1
  seg_idx <- 1
  
  # Process each edge
  for (edge in 1:length(tree$maps)) {
    edge_data <- tree$maps[[edge]]
    states <- names(edge_data)
    
    # Process transitions
    if (length(states) > 1) {
      distance <- 0
      for (i in seq_along(edge_data)) {
        distance <- distance + edge_data[i]
        if (i < length(edge_data)) {
          transitions_df$EdgeNumber[trans_idx] <- edge
          transitions_df$Transition01or10[trans_idx] <- paste0(states[i], "to", states[i+1])
          transitions_df$DistanceFromRootwardNode[trans_idx] <- distance
          trans_idx <- trans_idx + 1
        }
      }
    }
    
    # Process segments
    start_position <- 0
    for (i in seq_along(edge_data)) {
      end_position <- start_position + edge_data[i]
      segments_df$EdgeNumber[seg_idx] <- edge
      segments_df$State[seg_idx] <- as.numeric(states[i])
      segments_df$Start[seg_idx] <- start_position
      segments_df$End[seg_idx] <- end_position
      seg_idx <- seg_idx + 1
      start_position <- end_position
    }
  }
  
  # Remove any empty rows if we over-allocated
  if (trans_idx <= n_transitions) {
    transitions_df <- transitions_df[1:(trans_idx-1), ]
  }
  if (seg_idx <= n_segments) {
    segments_df <- segments_df[1:(seg_idx-1), ]
  }
  
  outlist$transitions_df <- transitions_df
  outlist$segments_df <- segments_df
  return(outlist)
}

# Optimized version using vectorized operations where possible
find_other_trait_state_optimized <- function(edge, transition_distance, other_segments) {
  # Filter to the specific edge first
  edge_segments <- other_segments[other_segments$EdgeNumber == edge, ]
  
  if (nrow(edge_segments) == 0) {
    return(NA)
  }
  
  # Find segment containing the transition point using vectorized comparison
  matching_idx <- which(edge_segments$Start < transition_distance & 
                       edge_segments$End >= transition_distance)
  
  if (length(matching_idx) > 0) {
    return(edge_segments$State[matching_idx[1]])
  } else {
    return(NA)
  }
}