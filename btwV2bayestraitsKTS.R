## Adapting btw v2 function "bayestraits" to use any version of BayesTraits
## Kate Snyder
## 6/27/2023
## Edited 6/29/2023 - variable logfile name; BTdirpath, BTdir added
## Edited 7/5/2023 - MCMC output files to dir instead of BTdir; added parse_scheduleKTS derived from btw::parse_schedule to be functional with BTv4
## Edited 7/6/2023 - OutputFolderPath and my_outdir, my_suboutdir, my_outdir_fullpath 
## Edited 7/17/2023 - new function: parse_TestPrior_log

bayestraitsKTS <- function (data = NULL, tree = NULL, commands = NULL, silent = TRUE, 
          remove_files = TRUE, BTversionNum = "V4", BTdirpath = NULL, OutputFolderPath = NULL) {
  require(btw)
  source("btwV2-parse_scheduleKTS.R")
  if (BTversionNum == "V3") {
    BTversion = "BayesTraitsV3"
  } else if (BTversionNum == "V4") {
    BTversion = "BayesTraitsV4"
  } else if (BTversionNum == "V2") {
    BTversion = "BayesTraitsV2"
  }
  if (is.null(BTdirpath)) {
    BTdir = getwd()
  } else {
    BTdir = BTdirpath
  }
  
  if (!inherits(data, "data.frame")) {
    stop("Data frame containing species data must be supplied")}
  if (inherits(tree, "phylo")) {
    treelabs <- tree$tip.label
  } else if (inherits(tree, "multiPhylo")) {
    treelabs = attributes(tree)$TipLabel
  } else {
    stop("Tree must be of class phylo or multiPhylo")
  }
  if (class(commands) != "character") {
    stop("Character vector containing BayesTraits commands must be supplied.")}
  if (!(class(data[[1]]) %in% c("character", "factor")))  {
    stop("First column of data must contain species names.")}
  if (length(setdiff(treelabs, data[[1]])) > 0) {
    stop(paste("No match found in the data:", paste(setdiff(tree$tip.label, 
                                                            data[[1]]), collapse = ", ")))}
  if (length(setdiff(data[[1]], treelabs)) > 0) {
    stop(paste("No match found in the phylogeny:", paste(setdiff(data[[1]], 
                                                                 tree$tip.label), collapse = ", ")))}
  if (length(setdiff(treelabs, data[[1]])) > 0 | length(setdiff(data[[1]], 
                                                                treelabs)) > 0) {
    stop("Species in your phylogeny and data must match up exactly.")}
  if (.Platform$OS.type == "windows") {
    windows = TRUE
  } else if (Sys.info()["sysname"] == "Darwin") {
    windows = FALSE
  } else {
    stop("Operating system not supported")
  }
  if (windows) {
    if (!(BTversion %in% list.files())) {
      stop(paste(BTversion, "is not in your current working directory."))}
  } else if (!(BTversion %in% list.files(BTdir)))  {
    stop(paste(BTversion, "is not in your designated working directory:", BTdir))}
  
  my_wd <- getwd()
  if (is.null(OutputFolderPath)) {
    my_outdir = my_wd
  } else if (!dir.exists(OutputFolderPath)) {
    dir.create(OutputFolderPath)
    my_outdir = paste0(OutputFolderPath)
  } else {
    my_outdir = paste0(OutputFolderPath)
  }
  if (commands[1] == "3") {
    ModelName = "Discrete-Dependent"
  } else if (commands[1] == "2") {
    ModelName = "Discrete-Independent"
  } 
  if (commands[2] == "1") {
    ModelMethod = "ML"
  } else if (commands[2] == "2") {
    ModelMethod = "MCMC"
  }
  my_suboutdir = paste(ModelName, ModelMethod, sep="_")
  my_suboutdir_fullpath = paste0(my_outdir,"/",my_suboutdir)
  print(my_suboutdir_fullpath)
  if (!dir.exists(my_suboutdir_fullpath)) {
    dir.create(my_suboutdir_fullpath)
  }
  
  #write(c(commands, "run"), file = "./inputfile.txt")
  #ape::write.nexus(tree, file = "./tree.nex", translate = T)
  #write.table(data, file = "./data.txt", quote = F, col.names = F,      row.names = F)
  write(c(commands, "run"), file = paste0(my_suboutdir_fullpath, "/inputfile.txt"))
  ape::write.nexus(tree, file = paste0(my_suboutdir_fullpath, "/tree.nex"), translate = T)
  write.table(data, file = paste0(my_suboutdir_fullpath, "/data.txt"), quote = F, col.names = F,  row.names = F)
  my_wd = my_suboutdir_fullpath
  
  
  if (windows) {
    if (silent) {
      invisible(shell(paste0("BayesTraits",version,".exe tree.nex data.txt < inputfile.txt"), intern = TRUE))
    }
    else {
      shell(paste0("BayesTraits",version,".exe tree.nex data.txt < inputfile.txt"))
    }
  }
  else {
    system(paste(paste0(BTdir, "/BayesTraits", BTversionNum), paste0(my_wd, "/tree.nex"), paste0(my_wd, "/data.txt"), paste0("< ", my_wd, "/inputfile.txt")), ignore.stdout = silent)
  }
  
  print(my_suboutdir_fullpath)
  print(my_wd)
  
  log <- "data.txt.Log.txt" %in% list.files(my_suboutdir_fullpath)
  logfile = "data.txt.Log.txt" # KTS added 
  if (log == FALSE) { # KTS added
    print(log)
   log <-  "data.txt.log.txt" %in% list.files(my_suboutdir_fullpath)
   logfile = "data.txt.log.txt"
   print("no log - is log TRUE now?")
   print(log)
  }
  
  schedule <- "data.txt.Schedule.txt" %in% list.files(my_wd)
  stones <- "data.txt.Stones.txt" %in% list.files(my_wd)
  ancstates <- "data.txt.AncStates.txt" %in% list.files(my_wd)
  output.trees <- "data.txt.Output.trees" %in% list.files(my_wd)
  varrates <- "data.txt.VarRates" %in% list.files(my_wd)
  TestPrior_log <- "LogAllOutputs.txt" %in% list.files(my_wd)
  if (!log) {
    print(getwd())
    print(my_suboutdir_fullpath)
    stop("Something went wrong: btw can't find a log file")
  }
  if (varrates) 
    warning("btw does not handle output from a variable rates model.")
  if (windows) {
    Log <- parse_log("data.txt.Log.txt")
    if (schedule) 
      Schedule <- parse_schedule("data.txt.Schedule.txt")
    else Schedule <- NULL
    if (stones) 
      Stones <- parse_stones("data.txt.Stones.txt")
    else Stones <- NULL
    if (ancstates) 
      AncStates <- parse_ancstates("data.txt.AncStates.txt")
    else AncStates <- NULL
    if (output.trees) 
      OutputTrees <- ape::read.nexus("data.txt.Output.trees")
    else OutputTrees <- NULL
  }
  else {
    #if (version %in% c("V3","V4")) {
    #  Log <- parse_log(paste0(dir, "/data.txt.Log.txt"))
    #} else if (version == "V2") {
    #  Log <- parse_log(paste0(dir, "/data.txt.log.txt"))
    #}
    Log <- parse_log(paste0(my_wd, "/", logfile)) # kts added/changed
    if (schedule) {
      Schedule <- parse_scheduleKTS(paste0(my_wd, "/data.txt.Schedule.txt"))
    } else {Schedule <- NULL}
    if (stones) {
      Stones <- parse_stones(paste0(my_wd, "/data.txt.Stones.txt"))
    } else {Stones <- NULL}
    if (ancstates) {
      AncStates <- parse_ancstates(paste0(my_wd, "/data.txt.AncStates.txt"))
    } else {AncStates <- NULL}
    if (output.trees) {
      OutputTrees <- ape::read.nexus(paste0(my_wd, "/data.txt.Output.trees"))
    } else {OutputTrees <- NULL}
    if (TestPrior_log) {
      TestPriorDF <- parse_TestPrior_log(file = paste0(my_wd, "/LogAllOutputs.txt"))
    } else {TestPriorDF <- NULL}
  }
  results <- list(Log = Log, Schedule = Schedule, Stones = Stones, 
                  AncStates = AncStates, OutputTrees = OutputTrees)
  if (remove_files) {
    if (windows) {
      shell(paste("DEL", "data.txt.Log.txt"))
      shell(paste("DEL", "data.txt"))
      shell(paste("DEL", "tree.nex"))
      shell(paste("DEL", "inputfile.txt"))
      if (schedule) 
        shell(paste("DEL", "data.txt.Schedule.txt"))
      if (stones) 
        shell(paste("DEL", "data.txt.Stones.txt"))
      if (ancstates) 
        shell(paste("DEL", "data.txt.AncStates.txt"))
      if (output.trees) 
        shell(paste("DEL", "data.txt.Output.trees"))
    }
    else {
      system(paste0("rm ", my_wd, "/", logfile)) # kts changed to logfile
      system(paste0("rm ", my_wd, "/data.txt"))
      system(paste0("rm ", my_wd, "/tree.nex"))
      system(paste0("rm ", my_wd, "/inputfile.txt"))
      if (schedule) 
        system(paste0("rm ", my_wd, "/data.txt.Schedule.txt"))
      if (stones) 
        system(paste0("rm ", my_wd, "/data.txt.Stones.txt"))
      if (ancstates) 
        system(paste0("rm ", my_wd, "/data.txt.AncStates.txt"))
      if (output.trees) 
        system(paste0("rm ", my_wd, "/data.txt.Output.trees"))
    }
  }
  return(results)
}



parse_TestPrior_log <-  function (file) {
  out <- scan(file = file, what = "c", quiet = T, sep = "\n")
  header_start <- grep("Sample form prior", out)
  
  LogLines = readLines(file)
  header_start <- grep("Sample form prior", LogLines)
  TestPriorLines <- LogLines[(header_start+1):(header_start+1000)]
  TestPriorMat = str_split(TestPriorLines, "\t", simplify = T)
  TestPriorSample = as.data.frame(TestPriorMat)
  return(TestPriorSample)
}
