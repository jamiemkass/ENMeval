#################################################
#########	CREATE MAXENT ARGUMENTS	#############
#################################################

make.args <- function(mod.settings, algorithm, labels=FALSE) {

  if(algorithm == "maxent.jar") {
    other.args <- c("noaddsamplestobackground", "noremoveDuplicates", "noautofeature")
    args.list <- list()
    rms <- mod.settings[["rms"]]
    fcs <- mod.settings[["fcs"]]
    
    for(i in 1:length(fcs)) {
      args.list[[i]] <- other.args
      if(!grepl("L", fcs[i])) args.list[[i]] <- c(args.list[[i]], "nolinear")
      if(!grepl("Q", fcs[i])) args.list[[i]] <- c(args.list[[i]], "noquadratic")
      if(!grepl("H", fcs[i])) args.list[[i]] <- c(args.list[[i]], "nohinge")
      if(!grepl("P", fcs[i])) args.list[[i]] <- c(args.list[[i]], "noproduct")
      if(!grepl("T", fcs[i])) args.list[[i]] <- c(args.list[[i]], "nothreshold")
    }
    
    rms.lab <- rep(rms, each=length(fcs))
    rms.arg <- paste("betamultiplier=", rms.lab, sep="")
    fcs.lab <- rep(fcs, times=length(rms))
    fcs.arg <- rep(args.list, times=length(rms))

    args <- list()
    for(i in 1:length(fcs.arg)) args[[i]] <- c(rms.arg[i], fcs.arg[[i]])
    args.lab <- list(fcs.lab, rms.lab)
  }
  
  if(algorithm == "maxnet") {
    rms <- mod.settings[["rms"]]
    fcs <- mod.settings[["fcs"]]
    
    fcs.arg <- as.list(tolower(rep(fcs, times=length(rms))))
    rms.arg <- as.list(sort(rep(rms, times=length(fcs))))
    args <- mapply(c, fcs.arg, rms.arg, SIMPLIFY=FALSE)
  }
  
  if(labels == FALSE) return(args) else return(args.lab)
	
}
