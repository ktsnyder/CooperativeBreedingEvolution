# File Output Helper Functions for Simmap Overlap Pipeline
# Provides standardized file naming and directory creation

#' Create standardized filename for simmap overlap outputs
#' 
#' @param trait1 First trait name
#' @param trait2 Second trait name
#' @param file_type Type of file (e.g., "overlap_counts", "transition_counts", "arrow_plot")
#' @param nsims_real Number of real simulations (optional)
#' @param nsims_dummy Number of dummy simulations (optional)
#' @param real_dummy "REAL" or "DUMMY" indicator (optional)
#' @param other_label Optional user-defined label (default NULL)
#' @param extension File extension (default "csv")
#' @return Standardized filename
createStandardizedFilename <- function(trait1, trait2, file_type, 
                                     nsims_real = NULL, nsims_dummy = NULL, 
                                     real_dummy = NULL, other_label = NULL, 
                                     extension = "csv") {
  
  # Start with trait names and file type
  filename_parts <- c(trait1, trait2, file_type)
  
  # Add REAL/DUMMY indicator if provided
  if (!is.null(real_dummy)) {
    filename_parts <- c(filename_parts, real_dummy)
  }
  
  # Add nsims - use both if provided, otherwise single value
  if (!is.null(nsims_real) && !is.null(nsims_dummy)) {
    filename_parts <- c(filename_parts, paste0(nsims_real, "_", nsims_dummy))
  } else if (!is.null(nsims_real)) {
    filename_parts <- c(filename_parts, nsims_real)
  } else if (!is.null(nsims_dummy)) {
    filename_parts <- c(filename_parts, nsims_dummy)
  }
  
  # Add other_label if provided
  if (!is.null(other_label) && other_label != "") {
    filename_parts <- c(filename_parts, other_label)
  }
  
  # Combine with underscores and add extension
  filename <- paste(filename_parts, collapse = "_")
  filename <- paste0(filename, ".", extension)
  
  return(filename)
}

#' Create output directory structure for simmap overlap analysis
#' 
#' @param base_dir Base directory (default "Simmap_Overlap_Outputs")
#' @param trait1 First trait name
#' @param trait2 Second trait name
#' @param nsims_real Number of real simulations
#' @param nsims_dummy Number of dummy simulations
#' @param other_label Optional user-defined label (default NULL)
#' @param create_subdirs Whether to create subdirectories (default TRUE)
#' @return List of directory paths
createOutputDirectory <- function(base_dir = "Simmap_Overlap_Outputs",
                                trait1, trait2, nsims_real, nsims_dummy,
                                other_label = NULL, create_subdirs = TRUE) {
  
  # Create timestamp
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  # Build main directory name
  dir_parts <- c(trait1, "vs", trait2, nsims_real, nsims_dummy)
  
  # Add other_label if provided
  if (!is.null(other_label) && other_label != "") {
    dir_parts <- c(dir_parts, other_label)
  }
  
  # Add timestamp
  dir_parts <- c(dir_parts, timestamp)
  
  # Create main directory name
  main_dir_name <- paste(dir_parts, collapse = "_")
  main_dir <- file.path(base_dir, main_dir_name)
  
  # Create main directory
  if (!dir.exists(main_dir)) {
    dir.create(main_dir, recursive = TRUE)
  }
  
  # Define subdirectories
  subdirs <- list(
    main = main_dir,
    simmaps = file.path(main_dir, "simmaps"),
    data = file.path(main_dir, "data"),
    plots = file.path(main_dir, "plots"),
    summaries = file.path(main_dir, "summaries")
  )
  
  # Create subdirectories if requested
  if (create_subdirs) {
    for (subdir in subdirs[-1]) {  # Skip main dir (already created)
      if (!dir.exists(subdir)) {
        dir.create(subdir)
      }
    }
  }
  
  return(subdirs)
}

#' Get file path using standardized naming
#' 
#' @param dirs Directory list from createOutputDirectory
#' @param subdir Which subdirectory to use ("simmaps", "data", "plots", "summaries")
#' @param ... Arguments to pass to createStandardizedFilename
#' @return Full file path
getStandardizedPath <- function(dirs, subdir, ...) {
  filename <- createStandardizedFilename(...)
  return(file.path(dirs[[subdir]], filename))
}

#' Extract other_label from existing filename or directory
#' 
#' @param path File or directory path
#' @param trait1 First trait name to remove
#' @param trait2 Second trait name to remove
#' @return Extracted other_label or NULL
extractOtherLabel <- function(path, trait1, trait2) {
  # Get just the filename/dirname without path
  basename_path <- basename(path)
  
  # Remove file extension if present
  name_parts <- tools::file_path_sans_ext(basename_path)
  
  # First remove the trait names (which might contain underscores)
  name_parts <- gsub(trait1, "TRAIT1PLACEHOLDER", name_parts)
  name_parts <- gsub(trait2, "TRAIT2PLACEHOLDER", name_parts)
  
  # Split by underscores
  parts <- strsplit(name_parts, "_")[[1]]
  
  # Remove known components including placeholders
  known_parts <- c("TRAIT1PLACEHOLDER", "TRAIT2PLACEHOLDER", "vs", "REAL", "DUMMY", 
                   "overlap", "counts", "transition", "arrow", "plot", "gray", "ns", 
                   "combined", "analysis", "simmap", "examples", "summary", "sims",
                   "with", "transitions")
  
  # Also remove numeric parts (nsims, timestamps)
  numeric_parts <- grepl("^[0-9]+$", parts)
  
  # Find parts that aren't known or numeric
  other_parts <- parts[!(parts %in% known_parts | numeric_parts)]
  
  # If we have remaining parts, join them
  if (length(other_parts) > 0) {
    return(paste(other_parts, collapse = "_"))
  } else {
    return(NULL)
  }
}