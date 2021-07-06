#Coded by Kate Snyder
#Started: 7/7/2020
#Last Edited: 7/25/2020

#setwd("~/Documents/Creanza Lab/Comparative Evolution/Cooperative Breeding")  #you want to set your working directory to whichever folder has the pdf in it

#you may need to install these, but they might not work on the older version of R you have
library(pdftools)
library(tidyverse)

cbtext0 <- pdf_text("Evolutionary routes to non-kin cooperative breeding in birds supp data 2.pdf") %>% readr::read_lines() #read the pdf. %>% is part of the tidyverse I think, and basically "pipes" something into something else, i.e. a function. This line is thus basically the same as "cbtext0 <- readr::read_lines(pdf_text("Prevalence of...")).


#all of the functions that begin with "str_" are part of the "stringr" package (part of the tidyverse), which is for working with character objects using RegEx (Regular Expression) format
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
cbtext2 <- read.csv("Riehl 2013 supp data text lines_removedblanks.csv")

###getting rid of spaces where possible - need to do because I eventually split based on spaces
cbtext2[,2] <- str_replace_all(cbtext2[,2], "Female only", "Female_only")
cbtext2[,2] <- str_replace_all(cbtext2[,2], "Male only",  "Male_only")


#cbtext <- cbtext2[11:length(cbtext2)] # take out the first 10 rows (not data)
# cbtext <- str_replace(cbtext, "1..Species", "Species")
# cbtext <- str_replace(cbtext, "2..Known", "KnownParentalCare")
# cbtext <- str_replace(cbtext, "3..Inferred", "InferredParentalCare")
# cbtext <- str_replace(cbtext, "4..Region", "Region")
# cbtext <- str_replace(cbtext, "5..Source", "Source")
# cbtext2[1,] <- str_replace(cbtext2[1,], "[:digit:]\\.\\.", ",")  # remove the "1." "2." etc from headers
cbtext_fixsources <- str_replace(cbtext2[,2], ", ", ",") # getting rid of spaces after commas, which appear in the "Source" column

cbtext <- str_squish(cbtext_fixsources) # removes all the many extra spaces

cbsplit <- str_split(cbtext, " ") # this turns the data into a list of vectors, one vector per list element. Vector elements separated by the spaces (and spaces are removed)

cbmat <- set.seed(10) # new matrix
for (i in 1:length(cbsplit)) {  # cycle through each list element
  temprow <- cbsplit[[i]][1:6] # pull out just the list element that we're on in the for loop. temprow is a vector
  if (str_detect(temprow[1], "\\.")) { #if the first vector element (Genus) contains a ".", replace with the last full Genus encountered in the for loop (set below, after the "else")
    temprow[1] <- lastfullgenus
  } else {
    lastfullgenus <- temprow[1]
  }
  cbmat <- rbind(cbmat,temprow) # add the new row to the bottom of the new matrix
}
colnames(cbmat) <- c("Genus","Species","KnownParentalCare","InferredParentalCare", "Region", "Source") # set the column names of the new matrix

# combine Genus and Species
SpeciesScientific <- paste(cbmat[,1], cbmat[,2], sep="_")
completedmat <- cbind(SpeciesScientific, cbmat)
rownames(completedmat) <- NULL


completeddf <- as.data.frame(completedmat)
completeddf$Source[which(completeddf$Source == "--")] <- NA #replace missing sources with NA

write.csv(completeddf, file = "Riehl2013_data.csv") # ta-da!
