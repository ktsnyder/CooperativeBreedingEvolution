## Adapting btw v2 function "bayestraits" to use any version of BayesTraits
## Kate Snyder
## 6/27/2023
## Edited 6/29/2023 - variable logfile name; BTdirpath, BTdir added


bayestraitsKTS <- function (data = NULL, tree = NULL, commands = NULL, silent = TRUE, 
          remove_files = TRUE, version = "V3", BTdirpath = NULL) {
  require(btw)
  if (version == "V3") {
    BTversion = "BayesTraitsV3"
  } else if (version == "V4") {
    BTversion = "BayesTraitsV4"
  } else if (version == "V2") {
    BTversion = "BayesTraitsV2"
  }
  if (is.null(BTdirpath)) {
    BTdir = getwd()
  } else {
    BTdir = BTdirpath
  }
  
  if (!inherits(data, "data.frame")) 
    stop("Data frame containing species data must be supplied")
  if (inherits(tree, "phylo")) {
    treelabs <- tree$tip.label
  }
  else if (inherits(tree, "multiPhylo")) {
    treelabs = attributes(tree)$TipLabel
  }
  else {
    stop("Tree must be of class phylo or multiPhylo")
  }
  if (class(commands) != "character") 
    stop("Character vector containing BayesTraits commands must be supplied.")
  if (!(class(data[[1]]) %in% c("character", "factor"))) 
    stop("First column of data must contain species names.")
  if (length(setdiff(treelabs, data[[1]])) > 0) 
    stop(paste("No match found in the data:", paste(setdiff(tree$tip.label, 
                                                            data[[1]]), collapse = ", ")))
  if (length(setdiff(data[[1]], treelabs)) > 0) 
    stop(paste("No match found in the phylogeny:", paste(setdiff(data[[1]], 
                                                                 tree$tip.label), collapse = ", ")))
  if (length(setdiff(treelabs, data[[1]])) > 0 | length(setdiff(data[[1]], 
                                                                treelabs)) > 0) 
    stop("Species in your phylogeny and data must match up exactly.")
  if (.Platform$OS.type == "windows") {
    windows = TRUE
  }
  else if (Sys.info()["sysname"] == "Darwin") {
    windows = FALSE
  }
  else {
    stop("Operating system not supported")
  }
  if (windows) {
    if (!(BTversion %in% list.files())) 
      stop(paste(BTversion, "is not in your current working directory."))
  }
  else if (!(BTversion %in% list.files(BTdir))) 
    stop(paste(BTversion, "is not in your designated working directory:", BTdir))
  
  dir <- getwd()
  #print(dir)
  write(c(commands, "run"), file = "./inputfile.txt")
  ape::write.nexus(tree, file = "./tree.nex", translate = T)
  write.table(data, file = "./data.txt", quote = F, col.names = F, 
              row.names = F)
  if (windows) {
    if (silent) {
      invisible(shell(paste0("BayesTraits",version,".exe tree.nex data.txt < inputfile.txt"), intern = TRUE))
    }
    else {
      shell(paste0("BayesTraits",version,".exe tree.nex data.txt < inputfile.txt"))
    }
  }
  else {
    system(paste(paste0(BTdir, "/BayesTraits", version), paste0(dir, "/tree.nex"), paste0(dir, "/data.txt"), paste0("< ", dir, "/inputfile.txt")), ignore.stdout = silent)
  }
  
  log <- "data.txt.Log.txt" %in% list.files()
  logfile = "data.txt.Log.txt" # KTS added 
  if (log == FALSE) { # KTS added
    print(log)
   log <-  "data.txt.log.txt" %in% list.files()
   logfile = "data.txt.log.txt"
  }
  
  schedule <- "data.txt.Schedule.txt" %in% list.files()
  stones <- "data.txt.Stones.txt" %in% list.files()
  ancstates <- "data.txt.AncStates.txt" %in% list.files()
  output.trees <- "data.txt.Output.trees" %in% list.files()
  varrates <- "data.txt.VarRates" %in% list.files()
  if (!log) {
    print(getwd())
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
    Log <- parse_log(paste0(dir, "/", logfile)) # kts added/changed
    if (schedule) 
      Schedule <- parse_schedule(paste0(BTdir, "/data.txt.Schedule.txt"))
    else Schedule <- NULL
    if (stones) 
      Stones <- parse_stones(paste0(BTdir, "/data.txt.Stones.txt"))
    else Stones <- NULL
    if (ancstates) 
      AncStates <- parse_ancstates(paste0(BTdir, "/data.txt.AncStates.txt"))
    else AncStates <- NULL
    if (output.trees) 
      OutputTrees <- ape::read.nexus(paste0(BTdir, "/data.txt.Output.trees"))
    else OutputTrees <- NULL
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
      system(paste0("rm ", dir, "/", logfile)) # kts changed to logfile
      system(paste0("rm ", dir, "/data.txt"))
      system(paste0("rm ", dir, "/tree.nex"))
      system(paste0("rm ", dir, "/inputfile.txt"))
      if (schedule) 
        system(paste0("rm ", BTdir, "/data.txt.Schedule.txt"))
      if (stones) 
        system(paste0("rm ", BTdir, "/data.txt.Stones.txt"))
      if (ancstates) 
        system(paste0("rm ", BTdir, "/data.txt.AncStates.txt"))
      if (output.trees) 
        system(paste0("rm ", BTdir, "/data.txt.Output.trees"))
    }
  }
  return(results)
}
