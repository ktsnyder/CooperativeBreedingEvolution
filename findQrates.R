########
#Coded by Kate T. Snyder
#Last Modified 5-4-2025
#Built using RStudio Version 1.1.453
#R Version 3.5.1
#
#phytools_0.6-44
#maps_3.3.0        ape_5.1 
########
#

#findQrates - to be used to set the rates of transition between states for the building of simmaps for brownie
#Can also be used to generate simmap plots output to a PDF with the computed rates with plot = TRUE

#Example:
#qoutputpass <- findQrates(columns = c("MonogamyOrNot"), plot = TRUE, cladesubsetcolumn = "oscine", cladesubsetvalue = "nonpasserine", newdata = FALSE, newtree = "2019-10-22matezilla2treeHack_nondicho.nex")
#qoutputNonpass <- findQrates(columns = c("MonogamyOrNot"), plot = TRUE, cladesubsetcolumn = "oscine", cladesubsetvalue = "Nonpasserine", newtree = "2019-10-22matezilla2treeHack_nondicho.nex")

findQrates <- function(columns, plot=FALSE, newtree = FALSE, newdata = FALSE, cladesubsetcolumn = NULL, cladesubsetvalue = NULL, otherlabel = NULL, GlobalQrates = NULL) {
require(phytools)
source(file = "subsettreedata.R")
  
    output <- list()
  

subsetoutput <- subsettreedata(columns[1], cladesubsetcolumn = cladesubsetcolumn, cladesubsetvalue = cladesubsetvalue, newdata = newdata, newtree=newtree)  

  tree <- subsetoutput$subsettree  
  df <- subsetoutput$subsetdf
  print(length(tree$tip.label))
  discretetraitvec <- df[,columns[1]]
  names(discretetraitvec) <- df[,1] 
  
  ##########
    ERmodel <- ace(discretetraitvec,tree, type="discrete",model = "ER")
    ARDmodel <- ace(discretetraitvec,tree, type="discrete",model = "ARD")
    anovaERARD <- anova(ERmodel,ARDmodel)
    testrates=c(ERmodel$rates,ERmodel$lik.anc[1,], ARDmodel$rates,ARDmodel$lik.anc[1,])
    output$anovaERARD <- anovaERARD
    output$ERrates <- paste("ERrates",ERmodel$rates)
    output$ERlikanc <- paste("ERlik.anc",ERmodel$lik.anc[1,])
    output$ARDrates <- paste("ARDrates",ARDmodel$rates)
    output$ARDlikanc <- paste("ARDlik.anc",ARDmodel$lik.anc[1,])
  
    
    matsize <- length(ARDmodel$rates)
    if(matsize == 1){matsize<- matsize+1}
    #create matrix of correct size; make diag negative; multiply by rates
    blankqrates <- matrix(rep(1,matsize^2),matsize,matsize)
    diag(blankqrates) <- -1
    if (anovaERARD$`Pr(>|Chi|)`[2] < 0.05) {
      qrates <- c(ARDmodel$rates[2],ARDmodel$rates[1])*blankqrates  #kts reversed 2/21/2020
      qratesModel = "ARD"
    } else {
      qrates <- c(ERmodel$rates)*blankqrates
      qratesModel = "ER"
    }
    #pulls the state names and sets as col and row names
    rownames(qrates) <- dimnames(ARDmodel$lik.anc)[[2]]
    colnames(qrates) <- dimnames(ARDmodel$lik.anc)[[2]]
    output$qratesModel = qratesModel
    
    if (!is.null(GlobalQrates)) {
      qrates <- GlobalQrates
    }
    
    if (plot == TRUE) {
      subsetoutput <- subsettreedata(columns, cladesubsetcolumn = cladesubsetcolumn, cladesubsetvalue = cladesubsetvalue, newdata = newdata, newtree=newtree) 
      df <- subsetoutput$subsetdf
      tree <- subsetoutput$subsettree
      # make OutputFiles folder if not present
      subDir <- "Outputs"
      if (!dir.exists(subDir)) {
        dir.create(subDir)
      }
      
      discretetraitsimmap <- make.simmap(tree,discretetraitvec,model = qratesModel, nsim = 3) #makes 3 simmaps for viewing purposes
      discretetraitsimmapQset <- make.simmap(tree,discretetraitvec,model = qratesModel, nsim = 3, Q = qrates)
      pdf(file = file.path(subDir, paste0(Sys.Date(),columns[1], columns[2],cladesubsetvalue, otherlabel, qratesModel, " egSimmaps.pdf")),height=12,width=6)
      layout(matrix(1:6,nrow = 2,ncol=3))
      for (i in 1:3) {
        simmap <- discretetraitsimmap[[i]]
        py = c("black","red")
        pynamed <- py
        names(pynamed) <- c(0,1)
        
        # Plot simmap 
        plotSimmap(simmap,fsize=0.2, lwd = 0.8, colors = pynamed)
        #add tips
        treetiplabels <- simmap$tip.label %in% names(discretetraitvec[discretetraitvec==1]) 
        tiplabels(pch=21,bg=py[as.numeric(treetiplabels)+1], col = py[as.numeric(treetiplabels)+1], cex=0.1)
        numrates1 <- lapply(simmap[[7]],round,digits=9)
        title(paste(" ", "\nmake.simmap","Qrates:", numrates1[2],numrates1[3]),cex.main = 0.5)
        
        simmapQset <- discretetraitsimmapQset[[i]]
        plotSimmap(simmapQset,fsize=0.2,lwd=0.8, colors = pynamed)
        #add tips
        treetiplabels <- simmapQset$tip.label %in% names(discretetraitvec[discretetraitvec==1]) 
        tiplabels(pch=21,bg=py[as.numeric(treetiplabels)+1], col = py[as.numeric(treetiplabels)+1], cex=0.1)
        numrates <- lapply(qrates,round,digits=9)
        title(main=paste(" ","\n",qratesModel, "model","Qrates (output for Brownie):",numrates[2],numrates[3]),cex.main = 0.5)
      }
      dev.off()
    }
    
  output$qrates <- qrates
  return(output)

}
