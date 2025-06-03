
# After detecting states for each trait
fs_states <- detect_unique_states(FSsimtrees[1:10])
coop_states <- detect_unique_states(Coopsimtrees[1:10])
terr_states <- detect_unique_states(Terrsimtrees[1:10])

# Create a list of state vectors
states_per_trait <- list(
  fs = fs_states,
  coop = coop_states,
  terr = terr_states
)

# Create the multidimensional matrix
overlap_matrix <- create_multidimensional_matrix(states_per_trait)

generate_transition_keys(fs_states, coop_states, terr_states, trait_names = c("FS", "CB", ""))

stateCounts = getTransitionStateCounts3(trait1_simtrees = trait1_simtrees, trait2_simtrees = trait2_simtrees, trait3_simtrees = trait3_simtrees, trait_names = c("FS", "CB", ""))


#### Helper functions below here ----

getTransitionStateCounts3 <- function(trait1_simtrees, trait2_simtrees, trait3_simtrees, 
                                      trait_names = c("Trait1", "Trait2", "Trait3")) {
  # use the minimum number of trees available across all three traits
  nTrees <- min(length(trait1_simtrees), length(trait2_simtrees), length(trait3_simtrees))
  
  # Detect states for each trait from the first tree of each type
  trait1_states <- detect_unique_states(trait1_simtrees[1])
  trait2_states <- detect_unique_states(trait2_simtrees[1])
  trait3_states <- detect_unique_states(trait3_simtrees[1])
  
  # Check that the first two traits are binary (2 states)
  if (length(trait1_states) != 2 || length(trait2_states) != 2) {
    stop("The first two traits must be binary (have exactly 2 states)")
  }
  
  # Generate all transition keys
  keys <- generate_transition_keys(trait1_states, trait2_states, trait3_states, trait_names)
  
  # Extract the keys for each trait's transitions
  trait1_keys <- keys$trait1_keys
  trait2_keys <- keys$trait2_keys 
  trait3_keys <- keys$trait3_keys
  
  # Load helper function for map overlap
  source("MapOverlapThree.R")
  
  # Prepare to store one row of output per simmap triple:
  results_list <- list()
  
  # Continue with the rest of the function...
  # Loop over each simmap triple (each "simulation")
  for (t in 1:nTrees) {
    # The rest of the function would be updated to use the dynamic keys
    # Instead of hard-coded keys and assumptions about state values

    # Load the current simmap trees for this simulation
    trait1_sim <- trait1_simtrees[[t]]
    trait2_sim <- trait2_simtrees[[t]]
    trait3_sim <- trait3_simtrees[[t]]
    
    # Calculate overlap matrices for each pair of traits and all three traits
    OverlapMatTrait1Trait2 <- Map.Overlap(trait1_sim, trait2_sim)
    OverlapMatTrait1Trait3 <- Map.Overlap(trait1_sim, trait3_sim)
    OverlapMatTrait2Trait3 <- Map.Overlap(trait2_sim, trait3_sim)
    OverlapTrait1Trait2Trait3 <- Map.Overlap.Three(trait1_sim, trait2_sim, trait3_sim)
    
    # Create named vectors for the pairwise overlaps
    ObsPropsTrait1Trait2 <- c()
    for (i in 1:length(trait1_states)) {
      for (j in 1:length(trait2_states)) {
        key <- paste0("ObsProp", trait_names[1], trait1_states[i], trait_names[2], trait2_states[j])
        ObsPropsTrait1Trait2[key] <- OverlapMatTrait1Trait2[i, j]
      }
    }
    
    ObsPropsTrait1Trait3 <- c()
    for (i in 1:length(trait1_states)) {
      for (j in 1:length(trait3_states)) {
        key <- paste0("ObsProp", trait_names[1], trait1_states[i], trait_names[3], trait3_states[j])
        ObsPropsTrait1Trait3[key] <- OverlapMatTrait1Trait3[i, j]
      }
    }
    
    ObsPropsTrait2Trait3 <- c()
    for (i in 1:length(trait2_states)) {
      for (j in 1:length(trait3_states)) {
        key <- paste0("ObsProp", trait_names[2], trait2_states[i], trait_names[3], trait3_states[j])
        ObsPropsTrait2Trait3[key] <- OverlapMatTrait2Trait3[i, j]
      }
    }
    
    # Create named vector for the three-way overlaps
    ObsPropsAll <- c()
    for (i in 1:length(trait1_states)) {
      for (j in 1:length(trait2_states)) {
        for (k in 1:length(trait3_states)) {
          key <- paste0("ObsProp", trait_names[1], trait1_states[i], trait_names[2], 
                        trait2_states[j], trait_names[3], trait3_states[k])
          ObsPropsAll[key] <- OverlapTrait1Trait2Trait3[i, j, k]
        }
      }
    }
    
    # Get the overall proportion of branch lengths in each state for each trait
    trait1_desc <- describe.simmap(trait1_sim)
    trait1_props <- setNames(
      trait1_desc$times["prop", ], 
      paste0("ObsProp", trait_names[1], names(trait1_desc$times["prop", ]))
    )
    
    trait2_desc <- describe.simmap(trait2_sim)
    trait2_props <- setNames(
      trait2_desc$times["prop", ], 
      paste0("ObsProp", trait_names[2], names(trait2_desc$times["prop", ]))
    )
    
    trait3_desc <- describe.simmap(trait3_sim)
    trait3_props <- setNames(
      trait3_desc$times["prop", ], 
      paste0("ObsProp", trait_names[3], names(trait3_desc$times["prop", ]))
    )
    
    # Get segment and transition dataframes (using your helper function)
    trait1_segs <- getSimmapSegments(trait1_sim)
    trait1_trans <- trait1_segs$transitions_df
    trait1_seg <- trait1_segs$segments_df
    
    trait2_segs <- getSimmapSegments(trait2_sim)
    trait2_trans <- trait2_segs$transitions_df
    trait2_seg <- trait2_segs$segments_df
    
    trait3_segs <- getSimmapSegments(trait3_sim)
    trait3_trans <- trait3_segs$transitions_df
    trait3_seg <- trait3_segs$segments_df
    
    # Initialize counters for transitions
    trait1_counts <- setNames(rep(0, length(trait1_keys)), trait1_keys)
    trait2_counts <- setNames(rep(0, length(trait2_keys)), trait2_keys)
    trait3_counts <- setNames(rep(0, length(trait3_keys)), trait3_keys)
    
    ## (A) Count Trait1 transitions
    if(nrow(trait1_trans) > 0) {
      for (i in 1:nrow(trait1_trans)) {
        row_i <- trait1_trans[i, ]
        if (is.na(row_i$Transition01or10)) next  # skip if no change
        
        # Get the trait2 state at the precise location on the branch:
        trait2_state <- trait2_seg$State[
          trait2_seg$EdgeNumber == row_i$EdgeNumber &
            trait2_seg$Start < row_i$DistanceFromRootwardNode &
            trait2_seg$End >= row_i$DistanceFromRootwardNode
        ]
        if(length(trait2_state) == 0) next
        
        # Get the trait3 state at that same point:
        trait3_state <- trait3_seg$State[
          trait3_seg$EdgeNumber == row_i$EdgeNumber &
            trait3_seg$Start < row_i$DistanceFromRootwardNode &
            trait3_seg$End >= row_i$DistanceFromRootwardNode
        ]
        if(length(trait3_state) == 0) next
        
        # Create the key for this transition
        key <- paste0(trait_names[1], row_i$Transition01or10, "in", 
                      trait_names[2], trait2_state[1], 
                      trait_names[3], trait3_state[1])
        
        if(key %in% names(trait1_counts)) {
          trait1_counts[key] <- trait1_counts[key] + 1
        }
      }
    }
    
    ## (B) Count Trait2 transitions
    if(nrow(trait2_trans) > 0) {
      for (i in 1:nrow(trait2_trans)) {
        row_i <- trait2_trans[i, ]
        if (is.na(row_i$Transition01or10)) next
        
        # Look up trait1 state at the moment of the trait2 change:
        trait1_state <- trait1_seg$State[
          trait1_seg$EdgeNumber == row_i$EdgeNumber &
            trait1_seg$Start < row_i$DistanceFromRootwardNode &
            trait1_seg$End >= row_i$DistanceFromRootwardNode
        ]
        if(length(trait1_state) == 0) next
        
        # And trait3 state:
        trait3_state <- trait3_seg$State[
          trait3_seg$EdgeNumber == row_i$EdgeNumber &
            trait3_seg$Start < row_i$DistanceFromRootwardNode &
            trait3_seg$End >= row_i$DistanceFromRootwardNode
        ]
        if(length(trait3_state) == 0) next
        
        # Create the key for this transition
        key <- paste0(trait_names[2], row_i$Transition01or10, "in", 
                      trait_names[1], trait1_state[1], 
                      trait_names[3], trait3_state[1])
        
        if(key %in% names(trait2_counts)) {
          trait2_counts[key] <- trait2_counts[key] + 1
        }
      }
    }
    
    # ## (C) Count Trait3 transitions
    # if(nrow(trait3_trans) > 0) {
    #   for (i in 1:nrow(trait3_trans)) {
    #     row_i <- trait3_trans[i, ]
    #     if (is.na(row_i$Transition01or10)) next
    # 
    #     # Get the state of trait1 and trait2 at the transition:
    #     trait1_state <- trait1_seg$State[
    #       trait1_seg$EdgeNumber == row_i$EdgeNumber &
    #         trait1_seg$Start < row_i$DistanceFromRootwardNode &
    #         trait1_seg$End >= row_i$DistanceFromRootwardNode
    #     ]
    #     if(length(trait1_state) == 0) next
    # 
    #     trait2_state <- trait2_seg$State[
    #       trait2_seg$EdgeNumber == row_i$EdgeNumber &
    #         trait2_seg$Start < row_i$DistanceFromRootwardNode &
    #         trait2_seg$End >= row_i$DistanceFromRootwardNode
    #     ]
    #     if(length(trait2_state) == 0) next
    # 
    #     # Create the key for this transition
    #     key <- paste0(trait_names[3], row_i$Transition01or10, "in",
    #                   trait_names[1], trait1_state[1],
    #                   trait_names[2], trait2_state[1])
    # 
    #     if(key %in% names(trait3_counts)) {
    #       trait3_counts[key] <- trait3_counts[key] + 1
    #     }
    #   }
    }
    
    ## Compute overall totals for each focal transition type (for expected counts)
    # For trait1 (e.g., FS)
    transition_types_trait1 <- unique(gsub("in.*$", "", trait1_keys))
    trait1_totals <- list()
    for(tt in transition_types_trait1) {
      key_pattern <- paste0("^", tt)
      trait1_totals[[tt]] <- sum(trait1_counts[grep(key_pattern, names(trait1_counts))])
    }
    
    # For trait2 (e.g., Coop)
    transition_types_trait2 <- unique(gsub("in.*$", "", trait2_keys))
    trait2_totals <- list()
    for(tt in transition_types_trait2) {
      key_pattern <- paste0("^", tt)
      trait2_totals[[tt]] <- sum(trait2_counts[grep(key_pattern, names(trait2_counts))])
    }
    
    # # For trait3 (e.g., Terr)
    # # Extract the transition types (e.g., "Terr1toTerr2", "Terr2toTerr1", etc.)
    # transition_types_trait3 <- unique(gsub("in.*$", "", trait3_keys))
    # trait3_totals <- list()
    # for(tt in transition_types_trait3) {
    #   key_pattern <- paste0("^", tt)
    #   trait3_totals[[tt]] <- sum(trait3_counts[grep(key_pattern, names(trait3_counts))])
    # }
    
    ## Compute expected counts for each transition key
    # For trait1 transitions:
    trait1_expected <- trait1_counts
    for (key in names(trait1_counts)) {
      # Extract the transition type (e.g., "FS0to1")
      trans_type <- gsub("in.*$", "", key)
      total_dir <- trait1_totals[[trans_type]]
      
      # Extract trait2 and trait3 states from the key
      pattern <- paste0(".*in", trait_names[2], "(.*)", trait_names[3], "(.*)")
      matches <- regmatches(key, regexec(pattern, key))[[1]]
      if (length(matches) >= 3) {
        trait2_state_val <- matches[2]
        trait3_state_val <- matches[3]
        
        # Compute the expected count based on proportions
        overlap_key <- paste0("ObsProp", trait_names[2], trait2_state_val, 
                              trait_names[3], trait3_state_val)
        if (overlap_key %in% names(ObsPropsTrait2Trait3)) {
          trait1_expected[key] <- total_dir * ObsPropsTrait2Trait3[overlap_key]
        } else {
          trait1_expected[key] <- NA
        }
      } else {
        trait1_expected[key] <- NA
      }
    }
    
    # For trait2 transitions:
    trait2_expected <- trait2_counts
    for (key in names(trait2_counts)) {
      # Extract the transition type (e.g., "Coop0to1")
      trans_type <- gsub("in.*$", "", key)
      total_dir <- trait2_totals[[trans_type]]
      
      # Extract trait1 and trait3 states from the key
      pattern <- paste0(".*in", trait_names[1], "(.*)", trait_names[3], "(.*)")
      matches <- regmatches(key, regexec(pattern, key))[[1]]
      if (length(matches) >= 3) {
        trait1_state_val <- matches[2]
        trait3_state_val <- matches[3]
        
        # Compute the expected count based on proportions
        overlap_key <- paste0("ObsProp", trait_names[1], trait1_state_val, 
                              trait_names[3], trait3_state_val)
        if (overlap_key %in% names(ObsPropsTrait1Trait3)) {
          trait2_expected[key] <- total_dir * ObsPropsTrait1Trait3[overlap_key]
        } else {
          trait2_expected[key] <- NA
        }
      } else {
        trait2_expected[key] <- NA
      }
    }
    
    # # For trait3 transitions:
    # trait3_expected <- trait3_counts
    # for (key in names(trait3_counts)) {
    #   # Extract the transition type (e.g., "Terr1toTerr2")
    #   trans_type <- gsub("in.*$", "", key)
    #   total_dir <- trait3_totals[[trans_type]]
    #   
    #   # Extract trait1 and trait2 states from the key
    #   pattern <- paste0(".*in", trait_names[1], "(.*)", trait_names[2], "(.*)")
    #   matches <- regmatches(key, regexec(pattern, key))[[1]]
    #   if (length(matches) >= 3) {
    #     trait1_state_val <- matches[2]
    #     trait2_state_val <- matches[3]
    #     
    #     # Compute the expected count based on proportions
    #     overlap_key <- paste0("ObsProp", trait_names[1], trait1_state_val, 
    #                           trait_names[2], trait2_state_val)
    #     if (overlap_key %in% names(ObsPropsTrait1Trait2)) {
    #       trait3_expected[key] <- total_dir * ObsPropsTrait1Trait2[overlap_key]
    #     } else {
    #       trait3_expected[key] <- NA
    #     }
    #   } else {
    #     trait3_expected[key] <- NA
    #   }
    # }
    
    ## Combine all results for this simulation into one named vector.
    outvec <- c(TreeNum = t,
                # Observed counts:
                trait1_counts,
                trait2_counts,
                trait3_counts,
                # Expected counts:
                setNames(trait1_expected, paste0(names(trait1_expected), "Expected")),
                setNames(trait2_expected, paste0(names(trait2_expected), "Expected")),
              #  setNames(trait3_expected, paste0(names(trait3_expected), "Expected")),
                # Overall proportions:
                trait1_props, 
                trait2_props,
                trait3_props,
                # Pairwise overlaps:
                ObsPropsTrait1Trait2,
                ObsPropsTrait1Trait3,
                ObsPropsTrait2Trait3,
                # Three-way overlaps:
                ObsPropsAll
    )
    
    results_list[[t]] <- outvec
  }
  
  ## Combine all simulation results into one data.frame.
  results_df <- do.call(rbind, lapply(results_list, function(x) as.data.frame(t(x), stringsAsFactors = FALSE)))
  
  ## Convert any numeric columns (they may be character because of the rbind)
  results_df[] <- lapply(results_df, function(col) as.numeric(as.character(col)))
  
  # Add metadata about the traits and states
  attr(results_df, "trait1_states") <- trait1_states
  attr(results_df, "trait2_states") <- trait2_states
  attr(results_df, "trait3_states") <- trait3_states
  attr(results_df, "trait_names") <- trait_names
  
  return(results_df)
}
    


generate_transition_keys <- function(trait1_states, trait2_states, trait3_states, trait_names = NULL) {
  # Set default trait names if not provided
  if (is.null(trait_names)) {
    trait_names <- c("Trait1", "Trait2", "Trait3")
  }
  
  # Validate inputs
  if (length(trait_names) != 3) {
    stop("Must provide exactly 3 trait names")
  }
  
  # For the first two traits (expected to be binary), create transition keys
  # Check if the first two traits are binary
  if (length(trait1_states) != 2 || length(trait2_states) != 2) {
    warning("First two traits should ideally be binary (have exactly 2 states)")
  }
  
  # Create empty vectors to store keys
  trait1_keys <- c()
  trait2_keys <- c()
  terr_keys <- c()
  
  # Name the traits for readability
  trait1_name <- trait_names[1]
  trait2_name <- trait_names[2]
  trait3_name <- trait_names[3]
  
  # For transitions in trait1 (e.g., Coop)
  for (from_state in trait1_states) {
    for (to_state in trait1_states) {
      if (from_state != to_state) {  # Only consider actual transitions
        transition <- paste0(from_state, "to", to_state)
        
        # For each state of trait2 (e.g., FS)
        for (state2 in trait2_states) {
          # For each state of trait3 (e.g., Terr)
          for (state3 in trait3_states) {
            key <- paste0(trait1_name, transition, "in", trait2_name, state2, trait3_name, state3)
            trait1_keys <- c(trait1_keys, key)
          }
        }
      }
    }
  }
  
  # For transitions in trait2 (e.g., FS)
  for (from_state in trait2_states) {
    for (to_state in trait2_states) {
      if (from_state != to_state) {  # Only consider actual transitions
        transition <- paste0(from_state, "to", to_state)
        
        # For each state of trait1 (e.g., Coop)
        for (state1 in trait1_states) {
          # For each state of trait3 (e.g., Terr)
          for (state3 in trait3_states) {
            key <- paste0(trait2_name, transition, "in", trait1_name, state1, trait3_name, state3)
            trait2_keys <- c(trait2_keys, key)
          }
        }
      }
    }
  }
  
  # For transitions in trait3 (e.g., Terr)
  # Create all possible transition pairs for trait3
  if (length(trait3_states) > 1) {
    terr_transition_types <- c()
    for (i in 1:length(trait3_states)) {
      for (j in 1:length(trait3_states)) {
        if (i != j) {
          terr_transition_types <- c(terr_transition_types, 
                                     paste0(trait3_states[i], "to", trait3_states[j]))
        }
      }
    }
    
    # For each transition type in trait3
    for (tt in terr_transition_types) {
      # For each state of trait1
      for (state1 in trait1_states) {
        # For each state of trait2
        for (state2 in trait2_states) {
          key <- paste0(trait3_name, tt, "in", trait1_name, state1, trait2_name, state2)
          terr_keys <- c(terr_keys, key)
        }
      }
    }
  }
  
  # Return a list of all key types
  return(list(
    trait1_keys = trait1_keys,
    trait2_keys = trait2_keys,
    trait3_keys = terr_keys,
    all_keys = c(trait1_keys, trait2_keys, terr_keys)
  ))
}

create_multidimensional_matrix <- function(states_per_trait) {
  # Check input
  if (!is.list(states_per_trait) || length(states_per_trait) == 0) {
    stop("states_per_trait must be a non-empty list of character vectors")
  }
  
  # Get dimensions for the array
  dims <- sapply(states_per_trait, length)
  
  # Create an empty array with the right dimensions
  result <- array(0, dim = dims)
  
  # Set proper dimnames for each dimension
  dimnames_list <- states_per_trait
  names(dimnames_list) <- paste0("trait", 1:length(states_per_trait))
  dimnames(result) <- dimnames_list
  
  # Print information about the created matrix
  dim_str <- paste(dims, collapse=" × ")
  cat("Created a", dim_str, "matrix for", length(dims), "traits\n")
  
  # Show state configuration
  for (i in 1:length(states_per_trait)) {
    cat("Trait", i, "states:", paste(states_per_trait[[i]], collapse=", "), "\n")
  }
  
  return(result)
}

detect_unique_states <- function(simmap_collection) {
  # Check if we have a single simmap or a collection
  if (inherits(simmap_collection, "simmap")) {
    # Handle a single simmap
    simmap_list <- list(simmap_collection)
  } else {
    # Handle a collection (list, multiSimmap, multiPhylo)
    simmap_list <- simmap_collection
  }
  
  # Container for all unique states
  all_states <- character(0)
  
  # Process each simmap in the collection
  for (i in 1:length(simmap_list)) {
    simmap <- simmap_list[[i]]
    
    # Extract states from each edge in the maps
    for (edge_map in simmap$maps) {
      edge_states <- names(edge_map)
      all_states <- union(all_states, edge_states)
    }
    
    # Optionally check mapped.edge column names if they exist
    if (!is.null(simmap$mapped.edge)) {
      if (!is.null(colnames(simmap$mapped.edge))) {
        mapped_states <- colnames(simmap$mapped.edge)
        all_states <- union(all_states, mapped_states)
      }
    }
  }
  
  # Sort the states (preserving character/numeric distinctions)
  if (all(grepl("^[0-9]+$", all_states))) {
    # If all states are numeric strings, sort numerically
    all_states <- all_states[order(as.numeric(all_states))]
  } else {
    # Otherwise, sort alphabetically
    all_states <- sort(all_states)
  }
  
  # Print a verification message
  cat("Detected states:", paste(all_states, collapse=", "), "\n")
  
  return(all_states)
}
  
  #' Get labels for trait states
  #' 
  #' @param trait Trait name
  #' @return A character vector of labels for this trait's states
  getLabels <- function(trait) {
    require(stringr)
    if (trait == "Final.polygyny") {
      return(c("Monogamy", "Polygyny"))
    } else if (grepl("coop", trait, ignore.case = TRUE)) {
      return(c("Non-Cooperative", "Cooperative"))
    } else if (str_detect(trait, "Kin")) {
      return(c("Non-kin", "Kin"))
    } else if (str_detect(trait, "Familial")) {
      return(c("Non-Familial Living", "Familial Living"))
    } else if (str_detect(trait, "Colonial")) {
      return(c("Non-Colonial", "Colonial"))
    } else if (str_detect(trait, "GroupsLargerThanPair")) {
      return(c("Asocial or pair", "Small or large groups"))
    } else if (str_detect(trait, "LongSocialBonds")) {
      return(c("Bonds last one season or less", "Multi-year bonds"))
    } else if (str_detect(trait, "MoreThanTwoCaretakers")) {
      return(c("Two or fewer caretakers", "More than two caretakers"))
    } else if (str_detect(trait, "TwoOrMoreCaretakers")) {
      return(c("Fewer than two caretakers", "Two or more caretakers"))
    } else if (str_detect(trait, "Asocial")) {
      return(c("Asocial", "Pair or group sociality"))
    } else if (str_detect(trait, "SeasonOrLonger")) {
      return(c("Bonds last less than one season", "Season or longer social bonds"))
    } else if (str_detect(trait, "LargestGroupSizes")) {
      return(c("Asocial, pair, or small groups", "Large groups"))
    } else if (str_detect(trait, "FemaleSong")) {
      return(c("Female Song Absent", "Female Song Present"))
    } else if (str_detect(trait, "HighConfidence_Coop")) {
      return(c("Non-Cooperative", "Cooperative"))
    } else {
      return(c(paste(trait, "0"), paste(trait, "1")))
    }
  }