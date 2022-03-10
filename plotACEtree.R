########
#Coded by Kate T. Snyder
#Last Modified 10-11-2020
#Built using RStudio Version 3.6.453
#R Version 4.0
#
#
#(not sure about package versions)
#mnormt_1.5-5    plyr_1.8.4   geiger_2.0.6   btw_0.1
#phytools_0.6-44   R.utils_2.6.0   nortest_1.0-4
#maps_3.3.0        ape_5.1      nlme_3.1-137   nortest_1.0-4
#BayesTraitsV2
########
##
##V2.1 (Feb 2020): made compatible with new subsetting method 
##5/6/2020: changed to "plotACEtree.R", added islog argument to be passed to subsettreedata, matingtree --> tree, $BirdtreeFormat --> [,1], matingtiplabels --> treetiplabels::: good, except for color bar placement and needs legend, including a way to set the legend labels
##10/8/2020 - removed setwd() to my comp evol folder
##10/10/2020 - added discretelabels to inputs for legend
##10/11/2020 - spit out ACE values; added ability to set model to something other than "ARD
##8/25/2021 - add "otherlabel" arg
##3/9/2022 - set default cladesubsetcolumn = NULL

#plot ancestral character estimation phylogenies ("heattrees")
#
# plotheattree(columns = c("PolygynyOrMonog"), passeriformesonly = TRUE, newdata = "20200108_Song Database Update 2019 - PostMerge.csv", newtree = "2019-10-22matezilla2treeHack_nondicho.nex")
# plotheattree(columns = c("PolygynyOrMonog"), passeriformesonly = FALSE, newdata = "20200108_Song Database Update 2019 - PostMerge.csv", newtree = "2019-10-22matezilla2treeHack_nondicho.nex")


plotACEtree <- function(columns, cladesubsetcolumn = NULL, cladesubsetvalue = NULL, newdata = FALSE, newtree = FALSE, islog = FALSE, discretelabels = NULL, discretemodel = "ARD", otherlabel = NULL) {
  require(phytools)

  # tipsize = 0.65
  # spfontsize = 0.25
  # fileheight = 20

  tipsize = 1
  spfontsize = 0.5
  fileheight = 20
  
source(file = "subsettreedata.R")
subsetoutput <- subsettreedata(columns = columns, cladesubsetcolumn = cladesubsetcolumn, cladesubsetvalue = cladesubsetvalue, newdata = newdata, newtree = newtree, islog = islog) 
  
tree <- subsetoutput$subsettree  
  df <- subsetoutput$subsetdf
  discretetraitvec <- df[,columns[1]]
  names(discretetraitvec) <- df[,1]
  
py <- c("black", "white")

if (islog != FALSE) {
  loglabel = "log"
} else { loglabel = NULL }

  pdf(file=paste(Sys.Date(),columns[1], loglabel, columns[2],cladesubsetvalue,otherlabel, "_ACEtree.pdf", sep=""), height = fileheight, width = 15)

  if (length(columns) > 1) {  
    continuoustraitvec <- df[,columns[2]]
    names(continuoustraitvec) <- df[,1]  #species names 
    
    contmap <- contMap(tree,x=continuoustraitvec,fsize=c(0.4,1),lwd=1,plot=FALSE,cols = heat.colors(100),label.offset=50) 
    plot(setMap(contmap,colors = heat.colors(100)),legend = FALSE,lwd=3.6, direction = "rightwards", srt = 40, label.offset = 5000,fsize=spfontsize, outline = FALSE)
  
  add.color.bar(23,heat.colors(100),title=paste(loglabel,columns[2]),lwd=4,fsize=0.90,subtitle="",prompt = FALSE, lims = c(min(continuoustraitvec),max(continuoustraitvec)), x = 0.5, y = 10)
   
  } else if (length(columns == 1)) {
  plotTree(tree,legend = FALSE,lwd=0.6, type = "fan", srt = 40, label.offset = 5000,fsize=spfontsize, outline = FALSE)
  }
  
  treetiplabels <- tree$tip.label %in% names(discretetraitvec[discretetraitvec==1])
  
  py <- c("black", "white")
  circles=ace(x=discretetraitvec,phy=tree,type="discrete",model=discretemodel)
  nodelabels(thermo=circles$lik.anc,piecol=py, height = 1.2, width = 1.2, horiz = TRUE, frame = "circle")
  
  print(circles)

  legend("bottomleft", legend = discretelabels, cex = 0.9, fill=py, bty="n")
  
  tiplabels(pch=21,bg=py[as.numeric(treetiplabels)+1], cex=tipsize)
  dev.off()
}
