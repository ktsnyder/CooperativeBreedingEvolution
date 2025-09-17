# Function to check LFS files
check_lfs_files <- function() {
  required_files <- list(
    list(
      name = "BirdzillaHackett4_Stage2_1000trees.tre",
      min_size_mb = 400,  # Actual ~480MB
      lfs_text = "version https://git-lfs"
    ),
    list(
      name = "PasserineMultiphy1000Hackett4_nondicho.nex", 
      min_size_mb = 100,  
      lfs_text = "version https://git-lfs"
    )
  )
  
  all_ok <- TRUE
  
  for (file_info in required_files) {
    if (!file.exists(file_info$name)) {
      cat("❌ ERROR: Required file", file_info$name, "not found\n")
      all_ok <- FALSE
    } else {
      # Check file size
      size_mb <- file.info(file_info$name)$size / (1024^2)
      
      if (size_mb < file_info$min_size_mb) {
        # Check if it's an LFS pointer
        first_line <- readLines(file_info$name, n = 1, warn = FALSE)
        if (grepl(file_info$lfs_text, first_line[1])) {
          cat("❌ ERROR:", file_info$name, "is a Git LFS pointer file, not the actual data file.\n")
          cat("   Found pointer file (", round(size_mb * 1024, 2), "KB) instead of full file (>", 
              file_info$min_size_mb, "MB)\n")
          cat("   Please download the full file from GitHub or clone with Git LFS installed.\n")
          cat("Alternatively, skip the sections that require these files: 'Re-generate consensus trees' and 'Multitree tests' ")
        } else {
          cat("⚠️  WARNING:", file_info$name, "seems unusually small (", round(size_mb, 2), 
              "MB, expected >", file_info$min_size_mb, "MB)\n")
        }
        all_ok <- FALSE
      } else {
        cat("✅", file_info$name, "found (", round(size_mb, 1), "MB)\n")
      }
    }
  }
  
  if (!all_ok) {
    stop("Required data files are missing or incomplete. See messages above for details.")
  }
  
  invisible(TRUE)
}

# Run the check at script start
cat("Checking required phylogenetic data files...\n")
check_lfs_files()
cat("All required files present and valid.\n\n")