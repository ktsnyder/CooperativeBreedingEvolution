#!/usr/bin/env Rscript

# Package Installation Checker for Cooperative Breeding Evolution Analysis
# This script checks for required packages across all R scripts in the project
# and provides options to install missing packages

# Define required packages (extracted from all R scripts in the project)
required_packages <- c(
  "ape",
  "cowplot", 
  "dplyr",
  "emmeans",
  "flextable",
  "ggplot2",
  "ggpubr",
  "ggrepel", 
  "grid",
  "gridExtra",
  "mnormt",
  "patchwork",
  "phylolm",
  "phylopath",
  "phytools",
  "R.utils",
  "stringr",
  "tidyr",
  "tidyverse"
)

bioconductor_packages <- c("graph", "RBGL")

# Function to check if packages are installed and get their versions
check_packages <- function(packages) {
  installed_info <- list()
  missing_packages <- character(0)
  
  for (pkg in packages) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      # Package is installed, get version
      version <- tryCatch({
        packageVersion(pkg)
      }, error = function(e) {
        "Version unavailable"
      })
      installed_info[[pkg]] <- as.character(version)
    } else {
      # Package is missing
      missing_packages <- c(missing_packages, pkg)
    }
  }
  
  return(list(installed = installed_info, missing = missing_packages))
}

# Function to write installed package info to file
write_package_info <- function(installed_info, filename = "installed_packages_info.txt") {
  con <- file(filename, "w")
  
  writeLines("=== INSTALLED PACKAGES AND VERSIONS ===", con)
  writeLines(paste("Generated on:", Sys.time()), con)
  writeLines("", con)
  
  if (length(installed_info) > 0) {
    for (pkg_name in sort(names(installed_info))) {
      writeLines(paste(pkg_name, ":", installed_info[[pkg_name]]), con)
    }
  } else {
    writeLines("No packages from the required list are currently installed.", con)
  }
  
  close(con)
  cat("Package information written to:", filename, "\n")
}

# Function to install missing packages
install_missing_packages <- function(missing_packages) {
  if (length(missing_packages) == 0) {
    cat("No packages need to be installed!\n")
    return(TRUE)
  }
  
  cat("Attempting to install missing packages...\n")
  success_count <- 0
  failed_packages <- character(0)
  
  for (pkg in missing_packages) {
    cat("Installing", pkg, "...\n")
    
    # Check if this is a Bioconductor package
    if (pkg %in% bioconductor_packages) {
      # Install BiocManager if needed
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        cat("Installing BiocManager first...\n")
        install.packages("BiocManager", quiet = TRUE)
      }
      
      # Install the Bioconductor package
      tryCatch({
        BiocManager::install(pkg, update = FALSE, ask = FALSE, quiet = TRUE)
        if (requireNamespace(pkg, quietly = TRUE)) {
          cat("✓", pkg, "installed successfully from Bioconductor\n")
          success_count <- success_count + 1
        } else {
          cat("✗", pkg, "installation failed (package not loadable)\n")
          failed_packages <- c(failed_packages, pkg)
        }
      }, error = function(e) {
        cat("✗", pkg, "installation failed:", e$message, "\n")
        
        # OS-specific advice for Bioconductor packages
        if (Sys.info()["sysname"] == "Windows") {
          cat("  Note: On Windows, you may need Rtools installed.\n")
          cat("  Download from: https://cran.r-project.org/bin/windows/Rtools/\n")
        }
        failed_packages <- c(failed_packages, pkg)
      })
      
    } else {
      # Regular CRAN package installation
      tryCatch({
        install.packages(pkg, dependencies = TRUE, quiet = TRUE)
        if (requireNamespace(pkg, quietly = TRUE)) {
          cat("✓", pkg, "installed successfully\n")
          success_count <- success_count + 1
        } else {
          cat("✗", pkg, "installation failed (package not loadable)\n")
          failed_packages <- c(failed_packages, pkg)
        }
      }, error = function(e) {
        cat("✗", pkg, "installation failed:", e$message, "\n")
        failed_packages <- c(failed_packages, pkg)
      })
    }
  }
  
  cat("\nInstallation Summary:\n")
  cat("Successfully installed:", success_count, "packages\n")
  if (length(failed_packages) > 0) {
    cat("Failed to install:", length(failed_packages), "packages:", paste(failed_packages, collapse = ", "), "\n")
  }
  
  return(length(failed_packages) == 0)
}

# Main execution
cat("=== PACKAGE DEPENDENCY CHECKER ===\n")

# Combine all packages for checking (but they'll be installed differently)
all_packages <- c(required_packages, bioconductor_packages)

cat("Checking", length(all_packages), "required packages...\n")
cat("(Including", length(bioconductor_packages), "Bioconductor packages)\n\n")

# Check package status
result <- check_packages(all_packages)

# Display installed packages
if (length(result$installed) > 0) {
  cat("INSTALLED PACKAGES (", length(result$installed), "):\n")
  for (pkg_name in sort(names(result$installed))) {
    cat("✓", pkg_name, ":", result$installed[[pkg_name]], "\n")
  }
  cat("\n")
}

# Write installed package info to file
write_package_info(result$installed)

# Display missing packages
if (length(result$missing) > 0) {
  cat("MISSING PACKAGES (", length(result$missing), "):\n")
  for (pkg in sort(result$missing)) {
    cat("✗", pkg, "\n")
  }
  cat("\n")
  
  # Ask user if they want to install missing packages
  user_input <- readline("Would you like to install the missing packages? (y/n): ")
  
  if (tolower(trimws(user_input)) %in% c("y", "yes")) {
    print("User responded yes, attempting installation.")
    install_missing_packages(result$missing)
    
    # Re-check packages after installation
    cat("\nRe-checking package status after installation...\n")
    final_result <- check_packages(required_packages)
    write_package_info(final_result$installed, "installed_packages_info_updated.txt")
    
    if (length(final_result$missing) == 0) {
      cat("🎉 All required packages are now installed!\n")
    } else {
      cat("⚠️  Some packages still missing:", paste(final_result$missing, collapse = ", "), "\n")
    }
  } else {
    cat("Installation skipped. You can run this script again later to install missing packages.\n")
  }
} else {
  cat("🎉 All required packages are already installed!\n")
}

cat("\nPackage check complete.\n")