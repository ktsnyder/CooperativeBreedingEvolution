#Coded by Kate Snyder
#Started: 7/7/2020
#Last Edited: 8/12/2021 - combine Clade names that contain " & " (e.g. Titryidae & Tyrannidae)

#setwd("Unaltered from publication")  #you want to set your working directory to whichever folder has the pdf in it

#you may need to install these, but they might not work on the older version of R you have
library(pdftools)
library(tidyverse)

# Check for required PDF file
pdf_file <- file.path("Unaltered from publication", "Jetz Supplemental Data.pdf")
if (!file.exists(pdf_file)) {
  stop("Required file '", pdf_file, "' not found")
}

cbtext0 <- pdf_text(pdf_file) %>% readr::read_lines() # Read PDF and convert to character vector


# String manipulation functions for text processing
removelines <- str_detect(cbtext0, "Page ") # make a logical vector of all lines that contain the word "Page " (i.e. the footer of each page that says the page number, etc)
cbtext1 <- cbtext0[!removelines] # subset to only those for which the above is FALSE

keeplines <- str_detect(cbtext1, "[aeiou]") # make a logical vector of all lines that contain a vowel
cbtext2 <- cbtext1[keeplines] #subsets to just those lines

removelines2 <- str_detect(cbtext2, "Source") # logical vector of each page header
cbtext3 <- cbtext2[!removelines2]

###getting rid of spaces where possible - need to do because I eventually split based on spaces
cbtext3 <- str_replace_all(cbtext3, " & ", "")
# cbtext2 <- str_replace_all(cbtext2, "Male only",  "Male_only")

cbtext <- cbtext3[12:length(cbtext3)] # take out the first 11 rows (not data)

cbtext <- str_squish(cbtext) # removes all the many extra spaces

cbsplit <- str_split(cbtext, " ") # this turns the data into a list of vectors, one vector per list element. Vector elements separated by the spaces (and spaces are removed)

cbmat <- set.seed(10) # new matrix
for (i in 1:length(cbsplit)) {  # cycle through each list element
  temprow <- cbsplit[[i]][1:7] # pull out just the list element that we're on in the for loop. temprow is a vector
  species <- paste(temprow[3], temprow[4], sep = "_")
  addrow <- c(species,temprow)
  cbmat <- rbind(cbmat,addrow) # add the new row to the bottom of the new matrix
}
colnames(cbmat) <- c("Species", "Group", "Clade", "NameGenus", "NameSpecies", "Realm", "System", "Source") # set column names of the new matrix

rownames(cbmat) <- NULL


completeddf <- as.data.frame(cbmat)
#completeddf$Source[which(completeddf$Source == "--")] <- NA #replace missing sources with NA

write.csv(completeddf, file = "Jetz2011_data_clean.csv") # ta-da!
