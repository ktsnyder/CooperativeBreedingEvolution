### Adapt btw::Discrete 
### Code imported from btw v1 
### Created 3/17/2022
### Kate Snyder
### Last edited 5/31/2023: add plotdiscrete from btw v1; add LRtestV1 from v1
### Last edited 6/26/2023: added ability to handle bayestraits V3 (may need to be slightly adjusted), Added Seed-setting

DiscreteKTS <- function (tree, data, mode = "ML", dependent = FALSE, res = NULL, 
          resall = NULL, mrca = NULL, fo = NULL, mlt = 10, it = 1e+05, 
          bi = 5000, sa = 100, pr = NULL, pa = NULL, hp = NULL, hpall = NULL, 
          rj = NULL, rjhp = NULL, silent = TRUE, KeepBTInputFiles = FALSE, Seed = NULL) 
{
  if (class(tree) == "phylo") {
    tree$node.label = NULL
    treelabs = tree$tip.label
  }
  else if (class(tree) == "multiPhylo") {
    treelabs = attributes(tree)$TipLabel
  }
  else {
    stop("Tree must be of class phylo or multiPhylo")
  }
  if (!(class(data[, 1]) %in% c("character", "factor"))) {
    stop("First column of data should contain species names.")
  }
  if (length(setdiff(treelabs, data[, 1])) > 0) {
    stop(paste("No match found in the data:", paste(setdiff(tree$tip.label, 
                                                            data[, 1]), collapse = ", ")))
  }
  if (length(setdiff(data[, 1], treelabs)) > 0) {
    stop(paste("No match found in the phylogeny:", paste(setdiff(data[, 
                                                                      1], tree$tip.label), collapse = ", ")))
  }
  if (length(setdiff(treelabs, data[, 1])) > 0 | length(setdiff(data[, 
                                                                     1], treelabs)) > 0) {
    stop("Species in your phylogeny and data must match up exactly.")
  }
  if (ncol(data) > 3) {
    stop("Too many columns in data: BayesTraits can only analyze one or two discrete traits.")
  }
  if (!exists(".BayesTraitsPath") | !file.exists(.BayesTraitsPath)) {
    stop("Must define '.BayesTraitsPath' to be the path to BayesTraitsV2 on your computer. For example: .BayesTraitsPath <- User/Desktop/BayesTraitsV2")
  }
  if (mode == "Bayesian") {
    mode = 2
  }
  else {
    mode = 1
  }
  if (ncol(data) == 2) {
    model = 1
  }
  else (model = 3)
  input = c(model, mode)
  if (!is.null(res)) {
    for (i in 1:length(res)) {
      input = c(input, paste("Restrict", res[i]))
    }
  }
  if (dependent == FALSE) {
    input = c(input, "res q12 q34")
    input = c(input, "res q21 q43")
    input = c(input, "res q13 q24")
    input = c(input, "res q31 q42")
  }
  if (!is.null(resall)) {
    input = c(input, paste("resall", resall))
  }
  
  if (!is.null(Seed)) { # kts added 6/27/23
    input = c(input, paste("Se", Seed))
  }
  
  if (!is.null(mrca)) {
    for (i in 1:length(mrca)) {
      input = c(input, paste("mrca", paste("mrcaNode", 
                                           i, sep = ""), mrca[i]))
    }
  }
  if (!is.null(fo)) {
    for (i in 1:length(fo)) {
      input = c(input, paste("Fossil", paste("fossilNode", 
                                             i, sep = ""), fo[i]))
    }
  }
  if (mode == 1) {
    input = c(input, paste("mlt", as.numeric(mlt)))
  }
  if (mode == 2) {
    input = c(input, paste("it", format(it, scientific = F)))
    input = c(input, paste("bi", format(bi, scientific = F)))
    input = c(input, paste("sa", format(sa, scientific = F)))
    if (!is.null(pr)) {
      for (i in 1:length(pr)) {
        input = c(input, paste("prior", pr[i]))
      }
    }
    if (!is.null(pa)) {
      input = c(input, paste("pa", pa))
    }
    if (!is.null(rj)) {
      input = c(input, paste("rj", rj))
    }
    if (!is.null(hp)) {
      for (i in 1:length(hp)) {
        input = c(input, paste("hp", hp[i]))
      }
    }
    if (!is.null(hpall)) {
      input = c(input, paste("Hpall", hpall))
    }
    if (!is.null(rjhp)) {
      input = c(input, paste("rjhp", rjhp))
    }
  }
  input = c(input, paste("lf ./BTout.txt"))  # kts deleted ".log.txt". # lf = log file
  input = c(input, "run")
  write(input, file = "./inputfile.txt")
  ape::write.nexus(tree, file = "./tree.nex", translate = T)
  write.table(data, file = "./data.txt", quote = F, col.names = F, 
              row.names = F)
  system(paste(.BayesTraitsPath, "./tree.nex", "./data.txt", 
               "< ./inputfile.txt"), ignore.stdout = silent)
  
  require(stringr)
  
  if (str_detect(.BayesTraitsPath, "V4")) {
    Skip = grep("Tree No", scan(file = "./BTout.txt.Log.txt", what = "c",  # kts capitalized Log
                                quiet = T, sep = "\n", blank.lines.skip = FALSE)) - 1
    Results = read.table("./BTout.txt.Log.txt", skip = Skip, sep = "\t", 
                         quote = "\"", header = TRUE)
    Results = Results[, -ncol(Results)]
    system(paste("rm ./BTout.txt.Log.txt"))

  } else if (str_detect(.BayesTraitsPath, "V2")) {
    Skip = grep("Tree No", scan(file = "./BTout.txt", what = "c",  
                                quiet = T, sep = "\n", blank.lines.skip = FALSE)) - 1
    Results = read.table("./BTout.txt", skip = Skip, sep = "\t", 
                         quote = "\"", header = TRUE)
    Results = Results[, -ncol(Results)]
    #system(paste("rm ./BTout.txt"))
  } else if (str_detect(.BayesTraitsPath, "V3")) {
    Skip = grep("Tree No", scan(file = "./BTout.txt.Log.txt", what = "c",  # kts capitalized Log
                                quiet = T, sep = "\n", blank.lines.skip = FALSE)) - 1
    Results = read.table("./BTout.txt.Log.txt", skip = Skip, sep = "\t", 
                         quote = "\"", header = TRUE)
    Results = Results[, -ncol(Results)]
    system(paste("rm ./BTout.txt.Log.txt"))
  }
  if (KeepBTInputFiles == FALSE) {
  system(paste("rm ./inputfile.txt"))
  system(paste("rm", "./tree.nex"))
  system(paste("rm", "./data.txt"))
  }
  return(Results)
}



plotdiscrete = function (model, estimates=TRUE, main=NULL) {
  
  # CREATE TRANSITION RATE MATRIX
  mat = matrix(0, 4, 4, dimnames=list(c("00", "01", "10", "11"), c("00", "01", "10", "11")))
  for (i in 1:4) {
    for (j in 1:4) {
      if (i != j  && sum(i,j) != 5) {
        mat[i,j] = mean(model[,grep(tail(paste("q", i, j, sep=""), 1), names(model))])
      } else mat[i,j] = NaN
    }
  }
  
  # PLOT
  rates <- labs <- round(c(mat[1,2],mat[2,1],mat[2,4],mat[4,2],mat[4,3],mat[3,4],mat[3,1],mat[1,3]),2) 
  if (estimates == F) {labs=rep("",8)}
  if (is.null(main)) {main = ""}
  plot(c(0,100), c(0,100), type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", main=main)
  text(x=c(20, 80, 20, 80), y=c(80, 80, 20, 20), labels=c("00", "01", "10", "11"), cex=4)
  text(x=c(50, 50, 93, 67, 50, 50, 7, 33), y=c(93, 67, 50, 50, 7,33, 50, 50), labels=labs, cex=0.75)
  arrows(x0=c(35, 65, 85, 75, 65, 35, 15, 25), y0=c(85, 75, 65, 35, 15, 25, 35, 65), x1=c(65, 35, 85, 75, 35, 65, 15, 25), y1=c(85, 75, 35, 65, 15, 25, 65, 35), lwd=rates/max(rates, na.rm=T)*15)
  
}


lrtestV1 = function (model1, model2) {
  if (nrow(model1) != nrow(model2)) {stop("Objects must contain the same number of models.")}
  Lhs = c(mean(model1$Lh), mean(model2$Lh))
  ind = sort(Lhs, index.return=TRUE)$ix
  max = which.max(c(model1$Lh, model2$Lh))
  Lh1 = Lhs[ind[1]]
  Lh2 = Lhs[ind[2]]
  LRstat = c()
  pval = c()
  for (n in 1:length(Lh1)) {
    lrs = 2*(Lh1[n] - Lh2[n])
    if (Lh1[n] < Lh2[n]) {lrs = -lrs}
    pv = pchisq(lrs, df=1, lower.tail=F) 
    LRstat = c(LRstat, lrs)
    pval = c(pval, pv)
  }
  return(data.frame(model1$Lh, model2$Lh, LRstat, pval))
}
