########
#Coded by Kate T. Snyder
#Last Modified 6-18-2024
#Built using R Version 4.0.2
#
#ape_5.7-1  phytools_1.9-16
########
########

################## Subset bird data and generate associated trees##################
 
#### Inputs: ----
## columns: vector of character values matching names of columns in dataframe. Will subset dataframe to only the subset of species which do not have NA values in any of the specified columns
## cladesubsetcolumn: a character object (aka a character vector of length 1) that is the column name of whichever column lists classification (e.g. Family) for a species
## cladesubsetvalue: a character object (aka a character vector of length 1) that is the clade/family (e.g. "Corvidae") that you want to analyze fdata for
## newdata: can be either a dataframe object (= the name of the dataframe in R environment, not in quotations) or a character object (= the name of the data file ending in .csv, in quotations)
## newtree: can be either an object of class "phylo" or a character object (= the name of the file ending in ".nex", in quotations)
## islog: defaults to FALSE, meaning that none of the data will be log-transformed. If we do want to log-transform some data, this should be a character object/vector (= the name of the column(s) of the data that you want to take the log() of, in quotations)
## suppresswarning: TRUE (default) or FALSE. Determines whether certain warnings that occur while running the code will be printed in the console or not
## skinnydata: TRUE or FALSE (default). Return dataframe subsetted just to columns mentioned plus species (including columns and cladesubsetcolumn)
 
#### Outputs (list of 2 items): ----
## subsettree - subsetted tree (.nex format)
## subsetdf - subsetted dataframe
## NonMatchedSpecies - vector of any species in data that are not in the tree

## 
#### Examples: ----
#outputall <- subsetbirddata(columns = c("PolygynyOrMonog", "EPP10threshold"), newdata = "20200123_SongDatabaseUpdate2019_R.csv", newtree = "2019-10-22matezilla2treeHack_nondicho.nex")
#outputnonpass <- subsetbirddata(columns = c("PolygynyOrMonog", "EPP10threshold"), cladesubsetcolumn = "oscine", cladesubsetvalue = "nonpasserine", newdata = "20200123_SongDatabaseUpdate2019_R.csv", newtree = "2019-10-22matezilla2treeHack_nondicho.nex")
#subsettreedata(columns = c("Final.polygyny", "Syllable.rep.final"), cladesubsetcolumn = "Family", cladesubsetvalue = c("Acrocephalidae", "Icteridae"), islog = "Syllable.rep.final", suppresswarning = FALSE, skinnydata = TRUE)


#### Function ----
subsettreedata <- function(columns = NULL, cladesubsetcolumn = NULL, cladesubsetvalue = NULL, newdata = FALSE, newtree = FALSE, islog = FALSE, suppresswarning = TRUE, skinnydata = FALSE) {
  require(phytools)
  
  output <- list()
  
  suppressWarnings(  # suppressing the warning that appears after last "if" if newdata is length > 1
    if (is.data.frame(newdata)) {
      alldatadfos <- newdata
    } else if (is.character(newdata)) {
      alldatadfos <- as.data.frame(read.csv(newdata))
    } else if (newdata == FALSE) {
      alldatadfos <- as.data.frame(read.csv("Data_R.csv"))
    }
  ) #end suppressWarnings
  
  
  
  if (is.character(newtree)) {
    fulltree <- read.nexus(newtree)
  } else if (is.list(newtree)) {
    fulltree <- newtree
  } else if (newtree == FALSE) {
    fulltree <- read.nexus("ConsensusPasserineTreeHackett4_1000_OscineSubset.nex")
  } else {
    print("The tree is not a tree.")
  }
  
  alldatadfos[,1] <- as.character(alldatadfos[,1])
  colnames(alldatadfos)[1] <- "species"
  
  if (class(fulltree) == "multiPhylo") {
    subsettedMultiTree <- fullMultiTree <- fulltree
    fulltree <- fulltree[[1]]
  } else {
    subsettedMultiTree = NULL
  }
  
  
  ## check that all species in dataframe match a species in phylogeny
  matchtree <- alldatadfos[,1] %in% fulltree$tip.label
  alldatadfos <- alldatadfos[which(matchtree),]
  nonmatchedspecies <- alldatadfos[,1][!matchtree]
  numnonmatching <- length(nonmatchedspecies)
  
  if (suppresswarning == FALSE) {
    print(paste("There were", numnonmatching, "species in the data table that did not match a species in the phylogeny. These species were removed from the data:"))
    print(paste(nonmatchedspecies))
  }
  
  if (is.null(cladesubsetcolumn)) {
    alldatadf <- alldatadfos
  } else {
    includespecies <- alldatadfos[,1][which(alldatadfos[,cladesubsetcolumn] %in% cladesubsetvalue)]
    includetips <- match(includespecies, fulltree$tip.label)
    includetree <- keep.tip(fulltree, tip=includetips)
    fulltree <- includetree
    includerows <- alldatadfos[,1] %in% includespecies
    alldatadf <- alldatadfos[which(includerows),]
  }
  
  ####
  subsetdf <- alldatadf
  subsettree <- fulltree  # if multitree input, this is one tree pared down by clade subset, won't actually use this
  
  if (!is.null(columns)) {
    for (i in 1:length(columns)) {
      whichcol <- columns[i]
      subsetdf <- subsetdf[!is.na(subsetdf[,whichcol]),]
      havedatavec <- subsettree$tip.label %in% as.character(subsetdf[,1])
      matingtips <- which(havedatavec == TRUE)
      dropfortree <- which(havedatavec == FALSE)
      subsettree <- drop.tip(subsettree, tip = dropfortree)
      if (!is.null(subsettedMultiTree)) {
        keepspecies <- subsetdf[,1]
        OneTreeFromMultiTreeInput <- subsettedMultiTree[[1]]
        AllSpecies <- OneTreeFromMultiTreeInput$tip.label
        DropSpeciesMultiTree <- AllSpecies[which(!AllSpecies %in% keepspecies)]
        subsettedMultiTree <- lapply(subsettedMultiTree, drop.tip, tip = DropSpeciesMultiTree)
        class(subsettedMultiTree) <- "multiPhylo"
      }
      
      #   print(subsettree)
    }  # end for (i in length(columns))
  }  # end if !is.null(columns)
  
  if (islog != FALSE) {
    for (j in 1:length(islog)) {  #islog should be a character vector with column names of columns to be log-transformed
      whichcol <- islog[j]
      subsetdf[,whichcol] <- log(subsetdf[,whichcol])
    }
  } # end if (islog != FALSE)
  
  if (skinnydata == TRUE) {
    subsetdf <- subsetdf[, c("species", cladesubsetcolumn, columns)]
  }
  
  
  ditree <- multi2di(subsettree)
  ditree$edge.length[ditree$edge.length == 0] <- 0.0000000000000000001
  
  output$subsettree <- ditree
  output$subsetdf <- subsetdf
  output$NonMatchedSpecies <- nonmatchedspecies
  
  if (class(subsettedMultiTree) == "multiPhylo") {
    output$subsettree <- subsettedMultiTree
  }
  
  return(output)
}


