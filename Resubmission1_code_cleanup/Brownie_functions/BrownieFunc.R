########
#Coded by Kate T. Snyder
#Last Modified 5-6-2020
#Forked from "MateSongFunc2.0.R" 2/11/2020
#Built using RStudio Version 1.1.453
#R Version 3.4.2 ?
#
#mnormt_1.5-5    plyr_1.8.4   geiger_2.0.6   btw_0.1
#phytools_0.6-44   R.utils_2.6.0   nortest_1.0-4
#maps_3.3.0        ape_5.1      nlme_3.1-137   nortest_1.0-4
#BayesTraitsV2
########

#Supplemental Code: BrownieFunc
#Performs phylanova, brownie, plots brownie rate distributions
#
#setwd("~/Creanza Lab/Comparative Evolution")
#
#
#5/6/2020 - renamed from "MateSysFunc" to BrownieFunc, removed heattree from this function, matingtree --> tree - uploaded to github
#NEXT: find brownie code to use (formerly in browniefunction?) 
#5/12/2020- removed "model = ARD" from args, 

# columns <- c("oscine", "oscine", "Order")
# values <- c("Oscine", "Nonpasserine", "Passeriformes")
# phylanovadf <- set.seed(10)
# for (i in 1:3) {
#   tempcolumn <- columns[i]
#   tempvalue <- values[i]
#   out1 <- matesysfunc(columns=c("PolygynyOrMonog","EPPcontinuous"), cladesubsetcolumn = tempcolumn, cladesubsetvalue = tempvalue, heattree = TRUE)
#   out2 <- matesysfunc(columns=c("MonogamyOrNot","EPPcontinuous"), cladesubsetcolumn = tempcolumn, cladesubsetvalue = tempvalue, heattree = TRUE)
#   out3 <- matesysfunc(columns=c("EPP10threshold","EPPcontinuous"), cladesubsetcolumn = tempcolumn, cladesubsetvalue = tempvalue, heattree = TRUE)
#   temprow <- c(tempvalue, out1$phylanovatitle, out1$phylanova$Pf, out2$phylanovatitle, out2$phylanova$Pf, out3$phylanovatitle, out3$phylanova$Pf)
#   phylanovadf <- rbind(phylanovadf, temprow)
#   }


browniefunc <- function(columns, cladesubsetcolumn = FALSE, cladesubsetvalue = NULL, brownie=FALSE, matensim = 1000, matingtree = NA) {
  require(R.utils)  #used in withTimeout etc
  require(phytools)
  require(ape)
  require(base)
  require(mnormt)
#  require(plyr)
#  require(geiger)
  
  source(file = "subsettreedata.R") 
 # source(file = "findQratesNewTree2.0.R") #seems to work
 # source(file = "browniefunctionNewTree2.0.R") #not done
 # source(file = "Supplement_signalfunc.R")
 # source(file = "Supplement_scatterboxes.R")
 # source(file = "Supplement_pglsResiduals.R")
  
  output <- list()
  
  subsetout <- subsetbirddata(columns = columns, cladesubsetcolumn = cladesubsetcolumn, cladesubsetvalue = cladesubsetvalue, newdata = FALSE, newtree = FALSE, islog = FALSE) 
  tree <- subsetout$subsettree
  df <- subsetout$subsetdf
  
discretetraitvec <- df[,columns[1]]
continuoustraitvec <- df[,columns[2]]
names(discretetraitvec) <- df[,1]  #species column
names(continuoustraitvec) <- df[,1]
  
    set.seed(10)
    phylanova <- phylANOVA(tree,discretetraitvec,continuoustraitvec, nsim=matensim)
    output$phylanovatitle <- paste(columns[1],columns[2],"N =", length(continuoustraitvec))
    output$phylanova <- phylanova
    
    # if (brownie == TRUE) {
    #   browniefunction(MateParam,SongParam, matensim = matensim, minmax = minmax, phylanova = phylanova$Pf)
    # } # end if brownie=T
    
  return(output)
}