## scratch process multitree brownie

require(dplyr)
library(tidyr)
require(ggplot2)

setwd("/Users/kate/Desktop/CooperativeBreedingEvolution")

#files = c("2024-04-28 brownie multitree 400trees 20simsPerTree _Hackett4Oscine_fulltreeQ_ HighConfidence_Coop Syll.song.final.csv", "2024-04-28 brownie multitree 400trees 20simsPerTree _Hackett4Oscine_fulltreeQ_ HighConfidence_Coop Syllable.rep.final.csv", "2024-04-28 brownie multitree 400trees 20simsPerTree _Hackett4Oscine_fulltreeQ_ HighConfidence_Coop Song.rep.final.csv")
files = c("2024-05-16 brownie multitree 300trees 20simsPerTree _Hackett4Oscine_fulltreeQ_UpdatedSongData_ HighConfidence_Coop Song.rep.final.csv", "2024-05-16 brownie multitree 300trees 20simsPerTree _Hackett4Oscine_fulltreeQ_UpdatedSongData_ HighConfidence_Coop Syllable.rep.final.csv")

for (tempfeat in c("Syllable.rep.final", "Song.rep.final")) {
  
#  tempfile1 = list.files(pattern = paste("2024-04-28 brownie multitree 400trees 20simsPerTree _Hackett4Oscine_fulltreeQ_ HighConfidence_Coop", tempfeat))
#  tempfile2 = list.files(pattern = paste("2024-04-27 brownie multitree 100trees 20simsPerTree _Hackett4Oscine_fulltreeQ_ HighConfidence_Coop", tempfeat))
  
#  print(tempfile1)
#  print(tempfile2)
  
#  browniedf1 = read.csv(tempfile1)
#  browniedf2 = read.csv(tempfile2)
#  browniedf= rbind(browniedf1, browniedf2)

  tempfile = files[which(str_detect(files, tempfeat))]
  browniedf = read.csv(tempfile)
  
 # browniedf = browniedf[which(browniedf$TreeNum != 1),] # from when there were 301 trees
  
  browniedf$ARDRate1.GreaterThan.ARDRate0 = 0
  browniedf$ARDRate1.GreaterThan.ARDRate0[which(browniedf$ARDRate1 > browniedf$ARDRate0)] = 1
  browniedf$SignificantPvalBrownie = browniedf$Pval < 0.05
  numSimsPerTreedf = browniedf %>% group_by(TreeNum) %>% count
  numSimsPerTree = numSimsPerTreedf$n[1]
  DiscreteTrait = browniedf$DiscreteTrait[1]
  ContinuousTrait = browniedf$ContinuousTrait[1]
  state0 = "Non-cooperative"
  state1 = "Cooperative"
  nTrees = length(browniedf$DiscreteTrait)/numSimsPerTree
  nTrees = length(unique(browniedf$TreeNum))
  
  ## test for outlier trees
  kruskal.test(Pval ~ TreeNum, data = browniedf)
  
  # Identify numerical columns; let's assume these columns need outlier analysis
  numerical_cols <- c("Pval", "ERRate", "ERloglik", "ARDRate0", "ARDRate1", "ARDloglik", "ARDsimmapQ0to1", "ARDsimmapQ1to0", "ERsimmapQ", "ARDvERsimmapQ.LRtestPval", "ARDRate1.GreaterThan.ARDRate0")
  
  # Adding outlier flags for each numerical column
  browniedf_outliers <- browniedf %>%
    mutate(across(all_of(numerical_cols), ~if_else(. < quantile(., 0.25) - 1.5 * IQR(.) | 
                                                     . > quantile(., 0.75) + 1.5 * IQR(.), TRUE, FALSE), .names = "outlier_{col}"))
  
  # Review a summary of outliers by TreeNum
  outlier_summary <- browniedf_outliers %>%
    group_by(TreeNum) %>%
    summarize(across(starts_with("outlier"), sum), .groups = "drop")
  
  # Printing the outlier summary
  outlier_summary$OutlierSums = rowSums(outlier_summary[,2:12])
  
  outlierTrees = outlier_summary$TreeNum[which(outlier_summary$OutlierSums > 80)]
  
  browniedf$TreeCheck = NA
  browniedf$TreeCheck[which(browniedf$TreeNum %in% outlierTrees)] = browniedf$TreeNum[which(browniedf$TreeNum %in% outlierTrees)]
  
  
  ## median values 
  # Selecting all columns except the ones to be excluded
  columns_to_exclude <- c("DiscreteTrait", "ContinuousTrait", "k2", "convergence", "simmapnumber", "TreeNum" , "ARDRate1.GreaterThan.ARDRate0", "SignificantPvalBrownie")
  columns_to_summarize <- setdiff(colnames(browniedf), columns_to_exclude)
  
  # Grouping by 'TreeNum' and summarizing the other columns
  result <- browniedf %>%
    group_by(TreeNum) %>%
    summarize(across(.cols = all_of(columns_to_summarize),
                     .fns = list(
                       #mean = ~mean(., na.rm = TRUE),
                       median = ~median(., na.rm = TRUE)#,
                       #max = ~max(., na.rm = TRUE),
                       #min = ~min(., na.rm = TRUE)
                     )), numSimsARDRate1.GreaterThan.ARDRate0 = sum(ARDRate1.GreaterThan.ARDRate0, na.rm=T), numSimsSignificantPvalBrownie = sum(SignificantPvalBrownie, na.rm=T),
              .groups = 'drop')
  
  result$FractionSims.ARDRate0.GreaterThan.ARDRate1 = (1-result$numSimsARDRate1.GreaterThan.ARDRate0/numSimsPerTree)
  result$FractionSims.SignificantPvalBrownie = (result$numSimsSignificantPvalBrownie/numSimsPerTree)
  
  print(head(result))
  
#  write.csv(result, file = paste0("brownie multitree summary HighConfidence_Coop ", tempfeat, " ", nTrees, "trees ", numSimsPerTree, "sims fullTreeQ.csv"), row.names = F)
  
  
#  pdf(file = paste0(Sys.Date(), "_brownie all multitree HighConfidence_Coop ", tempfeat, " ", nTrees, "trees ", numSimsPerTree, "sims fullTreeQ hists.pdf"), height = 8, width = 8)
#  par(mar = c(4,4,2,1))
#  par(mfrow = c(2,2)) 
  
  MedianFractionSignificant = median(result$FractionSims.SignificantPvalBrownie)
  MedianFractionSims.ARDRate0.GreaterThan.ARDRate1 = median(result$FractionSims.ARDRate0.GreaterThan.ARDRate1)
  
  # plot1 - rates (brownie)
  for (i in 1:2) {
    if (i == 1) {
      pdf(file = paste0(Sys.Date(), "_brownie all multitree HighConfidence_Coop ", tempfeat, " ", nTrees, "trees ", numSimsPerTree, "sims fullTreeQ RateHistPlotOnly.pdf"), height = 5, width = 5)
      par(mar = c(4,4,2,1))
      par(mfrow = c(1,1)) 
    } else if (i == 2) {
      pdf(file = paste0(Sys.Date(), "_brownie all multitree HighConfidence_Coop ", tempfeat, " ", nTrees, "trees ", numSimsPerTree, "sims fullTreeQ hists.pdf"), height = 8, width = 8)
      par(mar = c(4,4,2,1))
      par(mfrow = c(2,2)) 
    }
    
  titlelabel <- paste(DiscreteTrait, ContinuousTrait, nTrees, "trees,", numSimsPerTree, "sims per tree\nPer tree: Median Fraction Sig =", MedianFractionSignificant, "| Median Fraction Rate0 > Rate1 =", MedianFractionSims.ARDRate0.GreaterThan.ARDRate1)
  
  D0 <- density(browniedf$ARDRate0)
  D1 <- density(browniedf$ARDRate1)
  DER <- density(browniedf$ERRate)
  ERresultD <- density(result$ERRate_median)
  
  length(unique(browniedf$ARDRate0))
  length(unique(browniedf$ERRate))
  
  par(cex.axis=0.8,
      cex.lab=0.8)
  
  plot(D0,col="#BFD3E6",
       xlim=c(min(c(D0$x,D1$x,DER$x)),
              max(c(D0$x,D1$x, DER$x))),
       ylim=c(min(c(D0$y,D1$y, DER$y)),
              max(c(D0$y,D1$y, DER$y))),
       main=titlelabel, 
       xlab = "", 
       ylab = "",
       cex.main = 0.6, 
       yaxt = 'n',
       lwd = 2)
  lines(D1, col="#1E00E6", lwd = 1)
  lines(DER, lty = 2)
 # lines(ERresultD, col = "green")
  title(xlab=paste("Rate of", ContinuousTrait,"evolution"),
        ylab= paste("Number of Observations across all simulations, all trees"), line = 2)
  legend("topright",legend = c(paste(state0),paste(state1), "Equal Rates"), lwd=c(2,1,1),col=c("#BFD3E6","#1E00E6", "black"), lty = c(1,1,2))
  axis(side=2, las=2)

  
  if (i == 1) {
    dev.off()
  }
  
} # end for i in 1:2 (make 2 PDFs)
  
  # plot 2 - pval density plot
  FracSignificantBrownie = sum(browniedf$Pval < 0.05)/length(browniedf$Pval)
  Dpval <- density(browniedf$Pval)
  plot(Dpval,col="black",
       xlim=c(min(Dpval$x),
              max(Dpval$x)),
       ylim=c(min(Dpval$y),
              max(Dpval$y)),
       main= paste("p-values from brownie, fraction significant:", FracSignificantBrownie),
       xlab = "", 
       ylab = "",
       cex.main = 0.6) 
  abline(v=0.05, col = "gray")
  title(xlab=paste("p-value"),
        ylab= paste("Number of Observations across all simulations, all trees"), line = 2)
  
  # plot 3 - density plot of Q rates ARD and ER
  title3label = paste("density plots of ace Q rates (ARD and ER) from each of", nTrees, "trees")
  Dard01 = density(result$ARDsimmapQ0to1_median)
  Dard10 = density(result$ARDsimmapQ1to0_median)
  DerQ = density(result$ERsimmapQ_median)
  plot(Dard01,col="#BFD3E6",
       xlim=c(min(c(Dard01$x,Dard10$x)),
              max(c(Dard01$x,Dard10$x))),
       ylim=c(min(c(Dard01$y,Dard10$y)),
              max(c(Dard01$y,Dard10$y))),
       main=title3label, 
       xlab = "", 
       ylab = "",
       cex.main = 0.6)
  lines(Dard10, col="#1E00E6")
  lines(DerQ)
  legend("topright",legend = c("ARD rate - 0 to 1", "ARD rate - 1 to 0", "ER Rate"), lwd=1,col=c("#BFD3E6","#1E00E6", "black"), lty = c(1,1,2))
  title(xlab=paste("Rate of transition between discrete states"),
        ylab= paste("Number of Observations across all trees"), line = 2)
  
  
  # plot 4 - log lik for ARD and ER rates
  FractionLRTestSignificant = sum(result$ARDvERsimmapQ.LRtestPval_median < 0.05)/nTrees
  title4label = paste("density plots of ace Q rate (ARD and ER) Log Likelihoods\n from each of", nTrees, "trees. Fraction significant:", FractionLRTestSignificant)
  Dard = density(result$ARDsimmapQ.LogLik_median)
  DerQLL = density(result$ERsimmapQ.LogLik_median)
  plot(Dard,col="blue",
       xlim=c(min(c(Dard$x,DerQLL$x)),
              max(c(Dard$x,DerQLL$x))),
       ylim=c(min(c(Dard$y,DerQLL$y)),
              max(c(Dard$y,DerQLL$y))),
       main=title4label, 
       xlab = "", 
       ylab = "",
       cex.main = 0.6)
  lines(DerQLL, col = "black")
  legend("topleft",legend = c("ARD Log.Lik", "ER Log.Lik"), lwd=1,col=c("blue", "black"), lty = c(1,1))
  title(xlab=paste("Log likelihood (ACE)"),
        ylab= paste("Number of Observations across all trees"), line = 2)
  
  dev.off()
  
} # end for tempfeat - this could be moved to after some of the next steps, currently just interested in the summary table


source("plotbrownie.R")
plotbrownie(data = browniedf, columns = c(DiscreteTrait, ContinuousTrait), discreteCategoryLabels = c("Non-cooperative", "Cooperative"), newpdf = T, otherlabel = "multitree", nsim = length(browniedf$DiscreteTrait), islog = ContinuousTrait)

# brownieall = browniedf
# browniedf = brownieall[which(brownieall$TreeNum == 233),]




# We first need to pivot the data to a long format where each row represents a median value for a variable
median_data <- result %>%
  pivot_longer(
    cols = ends_with("median"),  # Selects columns that end with 'median'
    names_to = "variable",
    values_to = "median_value"
  )

# Now, plot the data using ggplot2
ggplot(median_data, aes(x = median_value)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black") +
  facet_wrap(~ variable, scales = "free") +
  theme_minimal() +
  labs(title = paste("Histograms of Median Brownie Values across", nTrees, "trees,", numSimsPerTree, "sims per tree", DiscreteTrait, ContinuousTrait),
       x = "Median Value",
       y = "Frequency")
