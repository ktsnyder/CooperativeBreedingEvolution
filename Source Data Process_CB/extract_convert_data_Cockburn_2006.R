#Coded by Kate Snyder
#Started: 7/7/2020
#Last Edited: 7/15/2021 - just to change name of pdf file

#setwd("~/Documents/Creanza Lab/Comparative Evolution/Cooperative Breeding")  #you want to set your working directory to whichever folder has the pdf in it

#you may need to install these, but they might not work on the older version of R you have
library(pdftools)
library(tidyverse)

cbtext0 <- pdf_text("Cockburn supp table.pdf") %>% readr::read_lines() #read the pdf. %>% is part of the tidyverse I think, and basically "pipes" something into something else, i.e. a function. This line is thus basically the same as "cbtext0 <- readr::read_lines(pdf_text("Prevalence of...")).


#all of the functions that begin with "str_" are part of the "stringr" package (part of the tidyverse), which is for working with character objects using RegEx (Regular Expression) format
removelines <- str_detect(cbtext0, ":") # make a logical vector of all lines that contain a colon (i.e. the family headers AND the footer of each page that says the page number, etc)
cbtext1 <- cbtext0[!removelines] # subset to only those for which the above is FALSE

keeplines <- str_detect(cbtext1, "[aeiou]") # make a logical vector of all lines that contain a vowel
cbtext2 <- cbtext1[keeplines] #subsets to just those lines

###getting rid of spaces where possible - need to do because I eventually split based on spaces
cbtext2 <- str_replace_all(cbtext2, "Female only", "Female_only")
cbtext2 <- str_replace_all(cbtext2, "Male only",  "Male_only")

cbtext <- cbtext2[11:length(cbtext2)] # take out the first 10 rows (not data)
cbtext <- str_replace(cbtext, "1. Species", "Genus Species")
cbtext <- str_replace(cbtext, "2. Known", "KnownParentalCare")
cbtext <- str_replace(cbtext, "3. Inferred", "InferredParentalCare")
cbtext <- str_replace(cbtext, "4. Region", "Region")
cbtext <- str_replace(cbtext, "5. Source", "Source")
cbtext_fixcols <- str_remove_all(cbtext, "[:digit:]\\. ")  # remove the "1." "2." etc from headers
cbtext_fixsources <- str_replace(cbtext_fixcols, ", ", ",") # getting rid of spaces after commas, which appear in the "Source" column

cbtext <- str_squish(cbtext_fixsources) # removes all the many extra spaces

cbsplit <- str_split(cbtext, " ") # this turns the data into a list of vectors, one vector per list element. Vector elements separated by the spaces (and spaces are removed)

cbmat <- set.seed(10) # new matrix
for (i in 2:length(cbsplit)) {  # cycle through each list element
  temprow <- cbsplit[[i]][1:6] # pull out just the list element that we're on in the for loop. temprow is a vector (only keeping first 6)
  if (str_detect(temprow[1], "\\.")) { #if the first vector element (Genus) contains a ".", replace with the last full Genus encountered in the for loop (set below, after the "else")
    temprow[1] <- lastfullgenus
  } else {
    lastfullgenus <- temprow[1]
  }
  cbmat <- rbind(cbmat,temprow) # add the new row to the bottom of the new matrix
}
colnames(cbmat) <- cbsplit[[1]][1:6] # set the first vector of the list (just the first 6 values) as the column names of the new matrix

# combine Genus and Species
SpeciesScientific <- paste(cbmat[,1], cbmat[,2], sep="_")
completedmat <- cbind(SpeciesScientific, cbmat)
rownames(completedmat) <- NULL


completeddf <- as.data.frame(completedmat)
completeddf$Source[which(completeddf$Source == "--")] <- NA #replace missing sources with NA

write.csv(completeddf, file = "Cockburn2006_data.csv") # ta-da!
