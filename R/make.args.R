#################################################
#########	CREATE MAXENT ARGUMENTS	#############
#################################################

make.args <- function(mod.settings, algorithm, labels=FALSE) {

  if(algorithm %in% c("maxent.jar", "maxnet")) {
    rm <- mod.settings[["rm"]]
    fc <- mod.settings[["fc"]]
    rm.lab <- rep(rm, each=length(fc))
    fc.lab <- rep(fc, times=length(rm))
    args.lab <- list(fc.lab, rm.lab)
  }
  
  if(algorithm == "maxent.jar") {
    other.args <- c("noaddsamplestobackground", "noremoveDuplicates", "noautofeature")
    args.list <- list()
    
    for(i in 1:length(fc)) {
      args.list[[i]] <- other.args
      if(!grepl("L", fc[i])) args.list[[i]] <- c(args.list[[i]], "nolinear")
      if(!grepl("Q", fc[i])) args.list[[i]] <- c(args.list[[i]], "noquadratic")
      if(!grepl("H", fc[i])) args.list[[i]] <- c(args.list[[i]], "nohinge")
      if(!grepl("P", fc[i])) args.list[[i]] <- c(args.list[[i]], "noproduct")
      if(!grepl("T", fc[i])) args.list[[i]] <- c(args.list[[i]], "nothreshold")
    }
    
    
    rm.arg <- paste("betamultiplier=", rm.lab, sep="")
    fc.arg <- rep(args.list, times=length(rm))
    args <- list()
    for(i in 1:length(fc.arg)) args[[i]] <- c(rm.arg[i], fc.arg[[i]])
  }

  if(algorithm == "maxnet") {
    fc.arg <- as.list(tolower(rep(fc, times=length(rm))))
    rm.arg <- as.list(sort(rep(rm, times=length(fc))))
    args <- mapply(c, fc.arg, rm.arg, SIMPLIFY=FALSE)
  }
  
  if(labels == FALSE) return(args) else return(args.lab)
	
}
