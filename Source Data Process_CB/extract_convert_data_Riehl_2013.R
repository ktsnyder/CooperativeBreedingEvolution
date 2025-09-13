#Coded by Kate Snyder
#Started: 7/7/2020
#Last Edited: 7/15/2021 - changed pdf input name, removed code at the bottom that was not related to Riehl

#setwd("Unaltered from publication")  #you want to set your working directory to whichever folder has the pdf in it

if (!requireNamespace("pdftools", quietly = TRUE)) {
  stop("Package 'pdftools' is required but not installed. Please run: 
  install.packages('pdftools')")
}
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  stop("Package 'tidyverse' is required but not installed. Please run: 
  install.packages('tidyverse')")
}

library(pdftools)
library(tidyverse)

# Check for required PDF file
pdf_file <- "Riehl supp data table.pdf"
if (!file.exists(pdf_file)) {
  stop("Required file '", pdf_file, "' not found")
}

cbtext0 <- pdf_text(pdf_file) %>% readr::read_lines() # Read PDF and convert to character vector


# String manipulation functions for text processing
removelines <- str_detect(cbtext0, ":") # make a logical vector of all lines that contain a space (i.e. all those with species names)
cbtext1 <- cbtext0[!removelines] # subset to only those for which the above is FALSE
cbtext1 <- cbtext1[19:length(cbtext1)] # take out the first 10 rows (not data)
cbdf0 <- as.data.frame(cbtext1)
colnames(cbdf0) <- "x"
cbdf1 <- separate(cbdf0, col = "x",into = c("SpeciesName", "Dispersal", "Frequency", "Data", "O.F", "Kin", "Social_breeding_system_when_cooperative", "Category", "Source","Extra1","Extra2","Extra3","Extra4", "Extra5","Extra6","Extra7","Extra8", "Extra9", "Extra10","Extra11","Extra12","Extra13","Extra14", "Extra15","Extra16","Extra17","Extra18","Extra19","Extra20","Extra21","Extra22","Extra23","Extra24", "Extra25","Extra26","Extra27","Extra28","Extra29","Extra30","Extra31","Extra32","Extra33","Extra34", "Extra35","Extra36","Extra37","Extra38", "Extra39", "Extra40","Extra41","Extra42","Extra43","Extra44", "Extra45","Extra46","Extra47","Extra48", "Extra49", "Extra50"), sep = "\\s\\s", extra = "merge", fill = "right")
write.csv(cbdf1, file = "Riehl 2013 supp data columns lines.csv")


#### For Riehl 2013, did the rest of the steps in excel
####
#keeplines <- str_detect(cbtext1, "[aeiou]") # make a logical vector of all lines that contain a vowel
#cbtext2 <- cbtext1[keeplines] #subsets to just those lines