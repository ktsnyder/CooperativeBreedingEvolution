## Modified CharacterSimmaps function
## Based on test_trait_overlap_simmaps.R
## Key modifications:
## 1. Saves Real and Dummy simmaps as RDS files with trait names and method indicator
## 2. Uses general variable names (trait1/trait2 instead of FS/Coop)
## 3. Adds calculate_transitions arg to optionally skip getTransitionStateCounts
## 4. Compatible with calcHuel_corrected.R
## Kate Snyder with assistance from Claude (Anthropic)
## 2025-06-27

library(phytools)
library(stringr)
library(tidyverse)

CharacterSimmaps_modified <- function(columns, df, tree, dummy, nsims, treelabel, 
                                    datalabel = NULL, 
                                    dummyMethod = c("simHistory", "makeSimmap"), 
                                    plotSampleSimmaps = FALSE, 
                                    cladesubsetvalue = NULL, 
                                    setQratesTree = NULL, 
                                    setQratesData = NULL, 
                                    columnForGlobalQ = NULL, 
                                    columnGlobalQrates = NULL,
                                    calculate_transitions = TRUE,
                                    save_simmaps = TRUE,
                                    output_dir = "Simmap Overlap Outputs",
                                    trait1_real_multisimmapRDS = NULL,
                                    trait2_real_multisimmapRDS = NULL,
                                    trait1_dummy_multisimmapRDS = NULL,
                                    trait2_dummy_multisimmapRDS = NULL,
                                    return_simmaps = FALSE,
                                    other_label = NULL,
                                    dirs = NULL) {
  
  # Source file output helpers if not already loaded
  if (!exists("createStandardizedFilename")) {
    source("file_output_helpers.R")
  }
  
  # Use provided dirs or create output directory if not provided
  if (is.null(dirs)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir)
    }
    # Create a simple dirs list for backward compatibility
    dirs <- list(
      main = output_dir,
      simmaps = output_dir,
      data = output_dir,
      plots = output_dir,
      summaries = output_dir
    )
  }
  
  # Handle default dummyMethod selection
  if (length(dummyMethod) > 1) {
    dummyMethod <- dummyMethod[1]
  }
  
  # Set default datalabel if not provided
  if (is.null(datalabel)) {
    datalabel = paste(columns[1], columns[2])
  }
  
  # Handle Q tree parameter
  if (is.null(setQratesTree)) {
    Qtree = tree
  } else if (is.character(setQratesTree)) {
    Qtree = read.nexus(setQratesTree)
  } else {
    Qtree = setQratesTree
  }
  
  # Handle Q data parameter
  if (is.null(setQratesData)) {
    Qdata = df
  } else if (is.character(setQratesData)) {
    Qdata = read.csv(setQratesData)
  } else {
    Qdata = setQratesData
  }
  
  # Source required functions (assuming they're in parent directory)
  source("subsettreedata.R")
  source("findQrates.R")
  
  # Find Q rates for trait1 (previously cooprates)
  trait1_rates <- findQrates(columns = columns[1], newdata = Qdata, newtree = Qtree)
  trait1_Q <- trait1_rates$qrates
  trait1_Q01 <- trait1_Q[3]
  trait1_Q10 <- trait1_Q[2]
  trait1_Anc = str_remove(trait1_rates$ARDlikanc, "ARDlik.anc ")
  trait1_Anc = as.numeric(trait1_Anc)
  names(trait1_Anc) <- c("0","1")
  
  # Find Q rates for trait2 (previously FSrates)
  trait2_rates <- findQrates(columns = columns[2], newdata = Qdata, newtree = Qtree)
  trait2_Q <- trait2_rates$qrates
  trait2_Q01 <- trait2_Q[3]
  trait2_Q10 <- trait2_Q[2]
  trait2_Anc = str_remove(trait2_rates$ARDlikanc, "ARDlik.anc ")
  trait2_Anc = as.numeric(trait2_Anc)
  names(trait2_Anc) <- c("0","1")
  
  # Subset tree and data to common species
  subsets <- subsettreedata(columns = columns, newdata = df, newtree = tree)
  subsetdf <- subsets$subsetdf
  subsettree <- subsets$subsettree
  
  # Get trait vectors
  trait2_vec <- subsetdf[,columns[2]]
  names(trait2_vec) <- subsetdf$species
  trait1_vec <- subsetdf[,columns[1]]
  names(trait1_vec) <- subsetdf$species
  
  # Apply global Q rates if provided
  if (is.matrix(columnGlobalQrates)) {
    if (columnForGlobalQ == 1) {
      trait1_Q <- columnGlobalQrates
      trait1_Q01 <- trait1_Q[3]
      trait1_Q10 <- trait1_Q[2]
    } else {
      trait2_Q <- columnGlobalQrates
      trait2_Q01 <- trait2_Q[3]
      trait2_Q10 <- trait2_Q[2]
    }
  }
  
  # Check if we're using pre-generated simmaps
  use_pregen_simmaps <- FALSE
  simmap_generation_method <- "drop_tips_then_simmap"
  
  if (dummy == FALSE && !is.null(trait1_real_multisimmapRDS) && !is.null(trait2_real_multisimmapRDS)) {
    use_pregen_simmaps <- TRUE
    simmap_generation_method <- "full-tree-generated"
  } else if (dummy == TRUE && !is.null(trait1_dummy_multisimmapRDS) && !is.null(trait2_dummy_multisimmapRDS)) {
    use_pregen_simmaps <- TRUE
    simmap_generation_method <- "full-tree-generated"
  }
  
  # Generate or load simmaps
  if (dummy == FALSE) {
    # Real data simmaps
    if (use_pregen_simmaps) {
      cat("Loading pre-generated real data simmaps...\n")
      
      # Load trait1 simmaps
      cat("Loading", columns[1], "simmaps from:", trait1_real_multisimmapRDS, "\n")
      trait1_data <- readRDS(trait1_real_multisimmapRDS)
      if (is.list(trait1_data) && "multisimmap" %in% names(trait1_data)) {
        trait1_simtrees <- trait1_data$multisimmap
      } else {
        trait1_simtrees <- trait1_data
      }
      
      # Load trait2 simmaps
      cat("Loading", columns[2], "simmaps from:", trait2_real_multisimmapRDS, "\n")
      trait2_data <- readRDS(trait2_real_multisimmapRDS)
      if (is.list(trait2_data) && "multisimmap" %in% names(trait2_data)) {
        trait2_simtrees <- trait2_data$multisimmap
      } else {
        trait2_simtrees <- trait2_data
      }
      
      # Check if we have valid simmap objects
      if (!inherits(trait1_simtrees, "multiSimmap") && !inherits(trait1_simtrees[[1]], "simmap")) {
        stop("trait1_real_multisimmapRDS does not contain valid simmap objects")
      }
      if (!inherits(trait2_simtrees, "multiSimmap") && !inherits(trait2_simtrees[[1]], "simmap")) {
        stop("trait2_real_multisimmapRDS does not contain valid simmap objects")
      }
      
      # Drop tips to common species if needed
      common_tips <- intersect(trait1_simtrees[[1]]$tip.label, trait2_simtrees[[1]]$tip.label)
      if (length(common_tips) < length(trait1_simtrees[[1]]$tip.label)) {
        cat("Dropping tips to", length(common_tips), "common species...\n")
        trait1_tips_to_drop <- setdiff(trait1_simtrees[[1]]$tip.label, common_tips)
        trait2_tips_to_drop <- setdiff(trait2_simtrees[[1]]$tip.label, common_tips)
        
        if (length(trait1_tips_to_drop) > 0) {
          trait1_simtrees <- lapply(trait1_simtrees, drop.tip.simmap, tip = trait1_tips_to_drop)
          class(trait1_simtrees) <- c("multiSimmap", "multiPhylo")
        }
        if (length(trait2_tips_to_drop) > 0) {
          trait2_simtrees <- lapply(trait2_simtrees, drop.tip.simmap, tip = trait2_tips_to_drop)
          class(trait2_simtrees) <- c("multiSimmap", "multiPhylo")
        }
      }
      
      # Update nsims to match loaded simmaps
      nsims <- length(trait1_simtrees)
      cat("Using", nsims, "pre-generated simmaps\n")
      
    } else {
      cat("Generating real data simmaps for", columns[2], "...\n")
      trait2_simtrees <- make.simmap(tree = subsettree, x = trait2_vec, model = "ARD", nsim = nsims, Q = trait2_Q)
      cat("Generating real data simmaps for", columns[1], "...\n")
      trait1_simtrees <- make.simmap(tree = subsettree, x = trait1_vec, model = "ARD", nsim = nsims, Q = trait1_Q)
    }
    datalabel = paste(datalabel, "REAL")
    
    # Save simmaps as RDS if requested
    if (save_simmaps) {
      # Save trait1 simmaps
      trait1_filename <- getStandardizedPath(dirs, "simmaps", 
                                           trait1 = columns[1], 
                                           trait2 = columns[2], 
                                           file_type = "simmaps",
                                           real_dummy = "REAL",
                                           nsims_real = nsims,
                                           other_label = other_label,
                                           extension = "rds")
      trait1_simmap_data <- list(
        multisimmap = trait1_simtrees,
        trait = columns[1],
        paired_trait = columns[2],
        tree_file = treelabel,
        subsettree = subsettree,
        subsetdf = subsetdf,
        Q_rates = trait1_Q,
        Q01 = trait1_Q01,
        Q10 = trait1_Q10,
        ancestral_state = trait1_Anc,
        method = simmap_generation_method,
        nsims = nsims,
        timestamp = format(Sys.time(), "%Y-%m-%d_%H%M")
      )
      saveRDS(trait1_simmap_data, trait1_filename)
      cat("Saved", columns[1], "simmaps to:", trait1_filename, "\n")
      
      # Save trait2 simmaps
      trait2_filename <- getStandardizedPath(dirs, "simmaps", 
                                           trait1 = columns[2], 
                                           trait2 = columns[1], 
                                           file_type = "simmaps",
                                           real_dummy = "REAL",
                                           nsims_real = nsims,
                                           other_label = other_label,
                                           extension = "rds")
      trait2_simmap_data <- list(
        multisimmap = trait2_simtrees,
        trait = columns[2],
        paired_trait = columns[1],
        tree_file = treelabel,
        subsettree = subsettree,
        subsetdf = subsetdf,
        Q_rates = trait2_Q,
        Q01 = trait2_Q01,
        Q10 = trait2_Q10,
        ancestral_state = trait2_Anc,
        method = "drop_tips_then_simmap",
        nsims = nsims,
        timestamp = format(Sys.time(), "%Y-%m-%d_%H%M")
      )
      saveRDS(trait2_simmap_data, trait2_filename)
      cat("Saved", columns[2], "simmaps to:", trait2_filename, "\n")
    }
    
    # Plot sample simmaps if requested
    if (plotSampleSimmaps) {
      eg_filename <- getStandardizedPath(dirs, "plots",
                                       trait1 = columns[1],
                                       trait2 = columns[2],
                                       file_type = "egSimmaps",
                                       real_dummy = "REAL",
                                       nsims_real = nsims,
                                       other_label = other_label,
                                       extension = "pdf")
      
      pdf(file = eg_filename, height=12, width=6)
      layout(matrix(1:6, nrow = 2, ncol=3))
      
      for (i in 1:3) {
        # Plot trait1 simmap
        simmap1 <- trait1_simtrees[[i]]
        py = c("black","red")
        pynamed <- py
        names(pynamed) <- c(0,1)
        
        plotSimmap(simmap1, fsize=0.2, lwd = 0.8, colors = pynamed)
        treetiplabels <- simmap1$tip.label %in% names(trait1_vec[trait1_vec==1])
        tiplabels(pch=21, bg=py[as.numeric(treetiplabels)+1], col = py[as.numeric(treetiplabels)+1], cex=0.1)
        numrates1 <- lapply(trait1_Q, round, digits=6)
        title(paste(" ", "\nmake.simmap","Qrates:", numrates1[2], numrates1[3], columns[1]), cex.main = 0.5)
        
        # Plot trait2 simmap
        simmap2 <- trait2_simtrees[[i]]
        plotSimmap(simmap2, fsize=0.2, lwd=0.8, colors = pynamed)
        treetiplabels2 <- simmap2$tip.label %in% names(trait2_vec[trait2_vec==1])
        tiplabels(pch=21, bg=py[as.numeric(treetiplabels2)+1], col = py[as.numeric(treetiplabels2)+1], cex=0.1)
        numrates2 <- lapply(trait2_Q, round, digits=6)
        title(main=paste(" ","\nARDmodel","Qrates:", numrates2[2], numrates2[3], columns[2]), cex.main = 0.5)
      }
      dev.off()
    }
    
  } else {
    # Dummy simmaps
    if (use_pregen_simmaps) {
      cat("Loading pre-generated dummy data simmaps...\n")
      
      # Load trait1 dummy simmaps
      cat("Loading", columns[1], "dummy simmaps from:", trait1_dummy_multisimmapRDS, "\n")
      trait1_data <- readRDS(trait1_dummy_multisimmapRDS)
      if (is.list(trait1_data) && "multisimmap" %in% names(trait1_data)) {
        trait1_simtrees <- trait1_data$multisimmap
      } else {
        trait1_simtrees <- trait1_data
      }
      
      # Load trait2 dummy simmaps
      cat("Loading", columns[2], "dummy simmaps from:", trait2_dummy_multisimmapRDS, "\n")
      trait2_data <- readRDS(trait2_dummy_multisimmapRDS)
      if (is.list(trait2_data) && "multisimmap" %in% names(trait2_data)) {
        trait2_simtrees <- trait2_data$multisimmap
      } else {
        trait2_simtrees <- trait2_data
      }
      
      # Check if we have valid simmap objects
      if (!inherits(trait1_simtrees, "multiSimmap") && !inherits(trait1_simtrees[[1]], "simmap")) {
        stop("trait1_dummy_multisimmapRDS does not contain valid simmap objects")
      }
      if (!inherits(trait2_simtrees, "multiSimmap") && !inherits(trait2_simtrees[[1]], "simmap")) {
        stop("trait2_dummy_multisimmapRDS does not contain valid simmap objects")
      }
      
      # Drop tips to common species if needed
      common_tips <- intersect(trait1_simtrees[[1]]$tip.label, trait2_simtrees[[1]]$tip.label)
      if (length(common_tips) < length(trait1_simtrees[[1]]$tip.label)) {
        cat("Dropping tips to", length(common_tips), "common species...\n")
        trait1_tips_to_drop <- setdiff(trait1_simtrees[[1]]$tip.label, common_tips)
        trait2_tips_to_drop <- setdiff(trait2_simtrees[[1]]$tip.label, common_tips)
        
        if (length(trait1_tips_to_drop) > 0) {
          trait1_simtrees <- lapply(trait1_simtrees, drop.tip.simmap, tip = trait1_tips_to_drop)
          class(trait1_simtrees) <- c("multiSimmap", "multiPhylo")
        }
        if (length(trait2_tips_to_drop) > 0) {
          trait2_simtrees <- lapply(trait2_simtrees, drop.tip.simmap, tip = trait2_tips_to_drop)
          class(trait2_simtrees) <- c("multiSimmap", "multiPhylo")
        }
      }
      
      # Update nsims to match loaded simmaps
      nsims <- length(trait1_simtrees)
      cat("Using", nsims, "pre-generated dummy simmaps\n")
      datalabel <- paste(datalabel, "DUMMY")
      
    } else {
      if (dummyMethod == "simHistory") {
        cat("Generating dummy simmaps using sim.history for", columns[1], "...\n")
        trait1_simtrees <- sim.history(tree = subsettree, Q = trait1_Q, nsim = nsims, anc = trait1_Anc)
        cat("Generating dummy simmaps using sim.history for", columns[2], "...\n")
        trait2_simtrees <- sim.history(tree = subsettree, Q = trait2_Q, nsim = nsims, anc = trait2_Anc)
        datalabel <- paste(datalabel, "DUMMYSimHist")
        
      } else if (dummyMethod == "makeSimmap") {
        # Make randomized versions of trait1 simmaps
        trait1_simtrees_rand <- list()
        cat("Generating dummy simmaps using makeSimmap for", columns[1], "...\n")
        for (j in 1:nsims) {
          trait1_vec_temp <- subsetdf[,columns[1]]
          trait1_vec_random <- sample(trait1_vec_temp)
          names(trait1_vec_random) <- subsetdf$species
          trait1_simtree <- make.simmap(tree = subsettree, x = trait1_vec_random, model = "ARD", nsim = 1, Q = trait1_Q)
          trait1_simtrees_rand[[j]] <- trait1_simtree
          if (j %% 10 == 0) cat(j, "")
        }
        trait1_simtrees <- trait1_simtrees_rand
        cat("\n")
        
        # Make randomized versions of trait2 simmaps
        trait2_simtrees_rand <- list()
        cat("Generating dummy simmaps using makeSimmap for", columns[2], "...\n")
        for (j in 1:nsims) {
          trait2_vec_temp <- subsetdf[,columns[2]]
          trait2_vec_random <- sample(trait2_vec_temp)
          names(trait2_vec_random) <- subsetdf$species
          trait2_simtree <- make.simmap(tree = subsettree, x = trait2_vec_random, model = "ARD", nsim = 1, Q = trait2_Q)
          trait2_simtrees_rand[[j]] <- trait2_simtree
          if (j %% 10 == 0) cat(j, "")
        }
        trait2_simtrees <- trait2_simtrees_rand
        cat("\n")
        
        datalabel <- paste(datalabel, "DUMMYResampledMkSimmap")
      }
    }
    
    # Save dummy simmaps as RDS if requested
    if (save_simmaps) {
      # Save trait1 dummy simmaps
      trait1_filename <- getStandardizedPath(dirs, "simmaps", 
                                           trait1 = columns[1], 
                                           trait2 = columns[2], 
                                           file_type = "simmaps",
                                           real_dummy = "DUMMY",
                                           nsims_dummy = nsims,
                                           other_label = other_label,
                                           extension = "rds")
      trait1_simmap_data <- list(
        multisimmap = trait1_simtrees,
        trait = columns[1],
        paired_trait = columns[2],
        tree_file = treelabel,
        subsettree = subsettree,
        subsetdf = subsetdf,
        Q_rates = trait1_Q,
        Q01 = trait1_Q01,
        Q10 = trait1_Q10,
        ancestral_state = trait1_Anc,
        method = ifelse(use_pregen_simmaps, simmap_generation_method, paste0("drop_tips_then_", dummyMethod)),
        nsims = nsims,
        timestamp = format(Sys.time(), "%Y-%m-%d_%H%M")
      )
      saveRDS(trait1_simmap_data, trait1_filename)
      cat("Saved", columns[1], "dummy simmaps to:", trait1_filename, "\n")
      
      # Save trait2 dummy simmaps
      trait2_filename <- getStandardizedPath(dirs, "simmaps", 
                                           trait1 = columns[2], 
                                           trait2 = columns[1], 
                                           file_type = "simmaps",
                                           real_dummy = "DUMMY",
                                           nsims_dummy = nsims,
                                           other_label = other_label,
                                           extension = "rds")
      trait2_simmap_data <- list(
        multisimmap = trait2_simtrees,
        trait = columns[2],
        paired_trait = columns[1],
        tree_file = treelabel,
        subsettree = subsettree,
        subsetdf = subsetdf,
        Q_rates = trait2_Q,
        Q01 = trait2_Q01,
        Q10 = trait2_Q10,
        ancestral_state = trait2_Anc,
        method = ifelse(use_pregen_simmaps, simmap_generation_method, paste0("drop_tips_then_", dummyMethod)),
        nsims = nsims,
        timestamp = format(Sys.time(), "%Y-%m-%d_%H%M")
      )
      saveRDS(trait2_simmap_data, trait2_filename)
      cat("Saved", columns[2], "dummy simmaps to:", trait2_filename, "\n")
    }
  }
  
  # Set up for output dataframe
  column1 <- columns[1]
  column2 <- columns[2]
  
  # Update Nspecies based on actual simmap tips if using pre-generated simmaps
  if (use_pregen_simmaps) {
    Nspecies <- length(trait1_simtrees[[1]]$tip.label)
    # Update subsettree to match the simmaps
    subsettree <- trait1_simtrees[[1]]
    # Create a simple subsetdf with just species names
    subsetdf <- data.frame(species = trait1_simtrees[[1]]$tip.label)
  } else {
    Nspecies <- length(subsetdf$species)
  }
  
  # Calculate overlap statistics
  set.seed(6000)
  dfout <- data.frame()
  
  cat("Calculating overlap statistics...\n")
  for (i in 1:nsims) {
    trait2_simtree1 <- trait2_simtrees[[i]]
    trait1_simtree1 <- trait1_simtrees[[i]]
    
    # Get proportions for each trait
    trait2_described <- describe.simmap(trait2_simtree1)
    prop_trait2_state0 <- trait2_described$times[2,1]
    prop_trait2_state1 <- trait2_described$times[2,2]
    
    trait1_described <- describe.simmap(trait1_simtree1)
    prop_trait1_state0 <- trait1_described$times[2,1]
    prop_trait1_state1 <- trait1_described$times[2,2]
    
    totaltime <- trait2_described$times[1,3]
    
    # Calculate overlap
    OverlapMat <- Map.Overlap(trait2_simtree1, trait1_simtree1)
    # Note: Map.Overlap returns matrix with trait2 in rows, trait1 in columns
    # Overlap_0_0 <- OverlapMat[1,1]  # trait2=0, trait1=0
    # Overlap_0_1 <- OverlapMat[1,2]  # trait2=0, trait1=1
    # Overlap_1_0 <- OverlapMat[2,1]  # trait2=1, trait1=0
    # Overlap_1_1 <- OverlapMat[2,2]  # trait2=1, trait1=1
    Overlap_0_0 <- OverlapMat[1,1]  # trait1=0, trait2=0 
    Overlap_0_1 <- OverlapMat[2,1]  # trait1=0, trait2=1 # KTS ordered so that traits are now in order in Overlap_trait1_trait2
    Overlap_1_0 <- OverlapMat[1,2]  # trait1=1, trait2=0
    Overlap_1_1 <- OverlapMat[2,2]  # trait1=1, trait2=1
    
    treenum <- i
    
    # Create row with generic column names for compatibility with calcHuel_corrected
    # Order must match original: ObsProp0Absent, ObsProp0Present, ObsProp1Absent, ObsProp1Present
    # temprow <- c(treenum, column1, column2, trait1_Q01, trait1_Q10, trait2_Q01, trait2_Q10, 
    #              Nspecies, prop_trait2_state0, prop_trait2_state1, prop_trait1_state0, prop_trait1_state1,
    #              Overlap_0_0, Overlap_1_0, Overlap_0_1, Overlap_1_1, totaltime)
    temprow <- c(treenum, column1, column2, trait1_Q01, trait1_Q10, trait2_Q01, trait2_Q10, Nspecies, prop_trait1_state0, prop_trait1_state1, prop_trait2_state0, prop_trait2_state1, Overlap_0_0, Overlap_1_0, Overlap_0_1, Overlap_1_1, totaltime) # KTS reordered trait1 and trait2
    
    dfout <- rbind(dfout, temprow)
  }
  
  # Convert to dataframe and add column names
  dfout <- as.data.frame(dfout)
  
  # Use old column names for compatibility with existing code
  # colnames(dfout) <- c("treenum", "column1", "column2", "coopQ01", "coopQ10", "FSQAbsPres", "FSQPresAbs", 
  #                      "Nspecies", "propFSabsent", "propFSpresent", "propNoncoop", "propCoop", 
  #                      "ObsProp0Absent", "ObsProp0Present", "ObsProp1Absent", "ObsProp1Present", "totaltime")
  colnames(dfout) <- c("treenum", "trait1", "trait2", "trait1_Q01", "trait1_Q10", "trait2_Q01", "trait2_Q10",
                       "Nspecies", "prop_trait1_state0", "prop_trait1_state1", "prop_trait2_state0", "prop_trait2_state1",
                       "Overlap_0_0", "Overlap_1_0", "Overlap_0_1", "Overlap_1_1", "totaltime")
  
  # Convert numeric columns from character to numeric
  numeric_cols <- c("treenum", "trait1_Q01", "trait1_Q10", "trait2_Q01", "trait2_Q10",
                    "Nspecies", "prop_trait1_state0", "prop_trait1_state1",
                    "prop_trait2_state0", "prop_trait2_state1",
                    "Overlap_0_0", "Overlap_1_0", "Overlap_0_1", "Overlap_1_1",
                    "totaltime")
  
  for (col in numeric_cols) {
    dfout[[col]] <- as.numeric(as.character(dfout[[col]]))
  }
  
  # Also add new generic column names for future compatibility
  # dfout$trait1 <- column1
  # dfout$trait2 <- column2
  # dfout$trait1_prop_state0 <- dfout$propNoncoop
  # dfout$trait1_prop_state1 <- dfout$propCoop
  # dfout$trait2_prop_state0 <- dfout$propFSabsent
  # dfout$trait2_prop_state1 <- dfout$propFSpresent
  # dfout$Overlap_0_0 <- dfout$ObsProp0Absent
  # dfout$Overlap_0_1 <- dfout$ObsProp1Absent
  # dfout$Overlap_1_0 <- dfout$ObsProp0Present
  # dfout$Overlap_1_1 <- dfout$ObsProp1Present
  
  # Calculate transition counts if requested
  if (calculate_transitions) {
    cat("Calculating transition state counts...\n")
    
    # Try to use optimized version if available
    if (file.exists("getTransitionStateCounts_optimized.R")) {
      source("getTransitionStateCounts_optimized.R")
      transStateCounts <- getTransitionStateCounts_optimized(
        trait1_simtrees = trait1_simtrees,
        trait2_simtrees = trait2_simtrees,
        trait1_name = columns[1],
        trait2_name = columns[2],
        show_progress = TRUE
      )
    } else if (file.exists(file.path("Simmap_Overlap_functions","getTransitionStateCounts_optimized.R"))) {
      source(file.path("Simmap_Overlap_functions","getTransitionStateCounts_optimized.R"))
      transStateCounts <- getTransitionStateCounts_optimized(
        trait1_simtrees = trait1_simtrees,
        trait2_simtrees = trait2_simtrees,
        trait1_name = columns[1],
        trait2_name = columns[2],
        show_progress = TRUE
      )
    } else {
      # Fall back to original function
      if (file.exists("find transition counts by state for 2 Discrete traits.R")) {
        source("find transition counts by state for 2 Discrete traits.R")
        transStateCounts = getTransitionStateCounts(Coopsimtrees = trait1_simtrees, FSsimtrees = trait2_simtrees)
      } else {
        print("Could not calculate transition state counts")
      }
    }
    
    dfout = merge(dfout, transStateCounts, by.x = "treenum", by.y = "TreeNum")
  }
  
  # Save final output
  real_dummy_indicator <- ifelse(dummy, "DUMMY", "REAL")
  nsims_param <- if(dummy) "nsims_dummy" else "nsims_real"
  output_filename <- getStandardizedPath(dirs, "data",
                                       trait1 = columns[1],
                                       trait2 = columns[2],
                                       file_type = "overlap_counts",
                                       real_dummy = real_dummy_indicator,
                                       nsims_real = if(!dummy) nsims else NULL,
                                       nsims_dummy = if(dummy) nsims else NULL,
                                       other_label = other_label,
                                       extension = "csv")
  write.csv(dfout, output_filename, row.names = FALSE)
  cat("Final output saved to:", output_filename, "\n")
  
  # Plot sample simmaps for both traits (Real and Dummy)
  if (save_simmaps && nsims >= 3) {
    cat("Creating simmap visualization PDFs...\n")
    
    # Define colors for states
    state_colors <- c("0" = "black", "1" = "red")
    
    # Plot trait1 simmaps
    trait1_pdf_filename <- getStandardizedPath(dirs, "plots",
                                              trait1 = columns[1],
                                              trait2 = columns[2],
                                              file_type = "simmap_examples_trait1",
                                              real_dummy = real_dummy_indicator,
                                              nsims_real = if(!dummy) nsims else NULL,
                                              nsims_dummy = if(dummy) nsims else NULL,
                                              other_label = other_label,
                                              extension = "pdf")
    
    pdf(trait1_pdf_filename, width = 8, height = 20)  # Extra tall for large phylogeny
    layout(matrix(1:7, ncol = 1))  # 6 plots + 1 legend
    
    # Set margins for better visualization
    par(mar = c(2, 2, 3, 2))
    
    # Plot first 3 simmaps for trait1
    for (i in 1:3) {
      plotSimmap(trait1_simtrees[[i]], 
                 colors = state_colors,
                 fsize = 0.4,  # Slightly larger font for readability
                 lwd = 1.2,
                 mar = c(1, 1, 2, 1))
      
      # Add title indicating which simmap this is
      title(main = paste0(columns[1], " - ", ifelse(dummy, "Dummy", "Real"), " Simmap #", i),
            cex.main = 1.2)
      
      # Add tip labels colored by state
      if (use_pregen_simmaps) {
        # Extract tip states from the simmap
        tip_states <- getStates(trait1_simtrees[[i]], "tips")
        tip_colors <- state_colors[as.character(tip_states)]
      } else {
        tip_states <- trait1_vec[match(trait1_simtrees[[i]]$tip.label, names(trait1_vec))]
        tip_colors <- state_colors[as.character(tip_states)]
      }
      tiplabels(pch = 21, bg = tip_colors, col = tip_colors, cex = 0.3)
    }
    
    # Add legend
    par(mar = c(1, 1, 1, 1))
    plot.new()
    legend("center", 
           legend = c(paste(columns[1], "= 0"), paste(columns[1], "= 1")),
           fill = c("black", "red"),
           title = "States",
           cex = 1.5,
           bty = "n")
    
    dev.off()
    cat("Saved", columns[1], "simmap examples to:", trait1_pdf_filename, "\n")
    
    # Plot trait2 simmaps
    trait2_pdf_filename <- getStandardizedPath(dirs, "plots",
                                              trait1 = columns[2],
                                              trait2 = columns[1],
                                              file_type = "simmap_examples_trait2",
                                              real_dummy = real_dummy_indicator,
                                              nsims_real = if(!dummy) nsims else NULL,
                                              nsims_dummy = if(dummy) nsims else NULL,
                                              other_label = other_label,
                                              extension = "pdf")
    
    pdf(trait2_pdf_filename, width = 8, height = 20)  # Extra tall for large phylogeny
    layout(matrix(1:7, ncol = 1))  # 6 plots + 1 legend
    
    # Set margins for better visualization
    par(mar = c(2, 2, 3, 2))
    
    # Plot first 3 simmaps for trait2
    for (i in 1:3) {
      plotSimmap(trait2_simtrees[[i]], 
                 colors = state_colors,
                 fsize = 0.4,  # Slightly larger font for readability
                 lwd = 1.2,
                 mar = c(1, 1, 2, 1))
      
      # Add title indicating which simmap this is
      title(main = paste0(columns[2], " - ", ifelse(dummy, "Dummy", "Real"), " Simmap #", i),
            cex.main = 1.2)
      
      # Add tip labels colored by state
      if (use_pregen_simmaps) {
        # Extract tip states from the simmap
        tip_states <- getStates(trait2_simtrees[[i]], "tips")
        tip_colors <- state_colors[as.character(tip_states)]
      } else {
        tip_states <- trait2_vec[match(trait2_simtrees[[i]]$tip.label, names(trait2_vec))]
        tip_colors <- state_colors[as.character(tip_states)]
      }
      tiplabels(pch = 21, bg = tip_colors, col = tip_colors, cex = 0.3)
    }
    
    # Add legend
    par(mar = c(1, 1, 1, 1))
    plot.new()
    legend("center", 
           legend = c(paste(columns[2], "= 0"), paste(columns[2], "= 1")),
           fill = c("black", "red"),
           title = "States",
           cex = 1.5,
           bty = "n")
    
    dev.off()
    cat("Saved", columns[2], "simmap examples to:", trait2_pdf_filename, "\n")
  }
  
  # Return results
  if (return_simmaps) {
    return(list(
      dfout = dfout,
      trait1_simmaps = trait1_simtrees,
      trait2_simmaps = trait2_simtrees
    ))
  } else {
    return(dfout)
  }
}
