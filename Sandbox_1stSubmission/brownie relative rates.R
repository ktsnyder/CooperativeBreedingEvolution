## Process brownie outputs into table
# 1/23/2024
# Kate Snyder
# 
# edited 2/1/2024 to include coop breed tie2noncoop, other song features
# 6/14/2024 - MESSED AROUND ADDING COLUMNS. MIGHT BREAK (I think just the lower section ## summary table by filelist is iffy now)
# 6/18/2024 - made into function



BrownieRelativeRates <- function(BrownieOutputFolder = "OutputFolder", otherlabel = "") {
  
  require(stringr)
  require(dplyr)
  
  source("test_trait_overlap_simmaps.R") 
  
  filelist = list.files(BrownieOutputFolder)
  filelist = filelist[which(str_detect(filelist, "rownie"))]
  filelist = filelist[which(str_detect(filelist, ".csv"))]
  
  browniesummary = set.seed(10)
  for (i in 1:length(filelist)) {
    tempfile = filelist[i]
    
    df = read.csv(file.path(BrownieOutputFolder, tempfile))
    Ncolspresent = sum(c("DiscreteTrait", "ContinuousTrait", "convergence", "ERloglik", "ARDloglik") %in% colnames(df), na.rm = T)
    
    if (Ncolspresent == 5) {
      
      temptrait = df$DiscreteTrait[1]
      tempsongtrait = df$ContinuousTrait[1]
      
      statelabels = getLabels(temptrait)
      
      #subset results for plotting
      browniedf <- df[df$convergence == "Optimization has converged.",]
      browniedf <- browniedf[!is.na(browniedf$convergence),]
      
      #calculate overall mean pval
      ERloglikmean <- mean(browniedf$ERloglik)
      ARDloglikmean <- mean(browniedf$ARDloglik)
      ERARDPval = round(pchisq(2*(ARDloglikmean-ERloglikmean),1,lower.tail=FALSE), digits = 4) #testing whether the two rates of continuous trait evolution are significantly different
      
      nsims = length(df[,1])
      n1greaterthan0 = sum(df$ARDRate0 < df$ARDRate1)
      frac1greaterthan0 = n1greaterthan0/nsims
      if (frac1greaterthan0 > 0.5) {
        fasterRate = statelabels[2]
        slowerRate = statelabels[1]
        fractionGreater = frac1greaterthan0
      } else {
        fasterRate = statelabels[1]
        slowerRate = statelabels[2]
        fractionGreater = 1-frac1greaterthan0
      }
      nPvalUnder0.05 = sum(df$Pval < 0.05)
      medianPval = median(df$Pval)
      
      if ("numSpecies" %in% colnames(df)) {
        numSpeciesInSubset = df$numSpecies[1]
      } else {
        numSpeciesInSubset = NA
      }
      
      if (str_detect(tempfile, "jacked")) {
        splitfilename = str_split(tempfile, "jacked")
        jackedfam = str_remove(splitfilename[[1]][2],".csv")
      } else {
        jackedfam = NA
      }
      
      temprow = c(tempfile, temptrait, tempsongtrait, nsims, fasterRate, slowerRate, fractionGreater, nPvalUnder0.05, medianPval, Nspecies, jackedfam, ERloglikmean, ARDloglikmean, ERARDPval)
      browniesummary = rbind(browniesummary, temprow)
      browniesummary = as.data.frame(browniesummary)
      colnames(browniesummary) = c("File","DiscreteTrait", "ContinuousTrait", "Nsims", "HigherRate", "LowerRate", "FractionOfSims", "num_Significant_Pval", "median_Pval", "numSpeciesInSubset", "removedFamily", "ERloglikmean", "ARDloglikmean", "BrowniePvalFromMeanLogLiks")
    } else {
      print("The file does not appear to contain all of the necessary columns, skipping.")
    }
  } 
  write.csv(browniesummary, paste0(Sys.Date()," Brownie summary table_", otherlabel, ".csv"), row.names = FALSE)
}
