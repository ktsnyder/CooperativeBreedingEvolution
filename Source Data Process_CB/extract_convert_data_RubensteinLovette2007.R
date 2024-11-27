#Coded by Kate Snyder
#Started: 7/7/2020
#Last Edited: 8/12/2021

#setwd("Unaltered from publication")  #you want to set your working directory to whichever folder has the pdf in it

#you may need to install these, but they might not work on the older version of R you have
library(pdftools)
library(tidyverse)

cbtext0 <- pdf_text("Rubenstein supp data.pdf") %>% readr::read_lines() #read the pdf. %>% is part of the tidyverse I think, and basically "pipes" something into something else, i.e. a function. This line is thus basically the same as "cbtext0 <- readr::read_lines(pdf_text("Prevalence of...")).

cbtext1 <- cbtext0[8:52] # take out the first 7 rows (not data) and lines 53:65 (not data)

cbtext1[1] <- str_replace(cbtext1[1], "Social System", "Social_System")

cbtext2 <- str_squish(cbtext1) # removes all the many extra spaces

cbsplit <- str_split(cbtext2, " ") # this turns the data into a list of vectors, one vector per list element. Vector elements separated by the spaces (and spaces are removed)

cbmat <- set.seed(10) # new matrix
for (i in 1:length(cbsplit)) {  # cycle through each list element
  temprow <- cbsplit[[i]] # pull out just the list element that we're on in the for loop. temprow is a vector
  temprow <- str_replace(temprow, "Cooperativea", "Cooperative")
  species <- paste(temprow[1], temprow[2], sep = "_")
  addrow <- c(species,temprow)
  cbmat <- rbind(cbmat,addrow) # add the new row to the bottom of the new matrix
}
cbmat1 <- cbmat[,c(1,4,5,6,7,8)]
colnames(cbmat1) <- c("Species", "Social_System", "Habitat", "Predictability (P)", "Constancy (C)", "Contingency (M)") # set column names of the new matrix

rownames(cbmat1) <- NULL


completeddf <- as.data.frame(cbmat1)
#completeddf$Source[which(completeddf$Source == "--")] <- NA #replace missing sources with NA

write.csv(completeddf, file = "RubensteinLovette2007_data_clean.csv") # ta-da!
