# Debug control functions for the simmap overlap pipeline
# Kate Snyder / Claude
# 2025-07-01

# Global debug flag - can be set by user
.SIMMAP_DEBUG <- FALSE

# Function to set debug mode
setDebugMode <- function(debug = TRUE) {
  .SIMMAP_DEBUG <<- debug
  if (debug) {
    cat("Debug mode enabled for simmap overlap pipeline\n")
  } else {
    cat("Debug mode disabled for simmap overlap pipeline\n")
  }
}

# Function to get current debug mode
getDebugMode <- function() {
  return(.SIMMAP_DEBUG)
}

# Debug print function - only prints if debug mode is on
debugPrint <- function(...) {
  if (.SIMMAP_DEBUG) {
    cat("DEBUG:", ..., "\n")
  }
}

# Validation debug function - always prints validation errors
validationPrint <- function(...) {
  cat("VALIDATION:", ..., "\n")
}