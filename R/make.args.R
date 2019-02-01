#################################################
#########	CREATE MAXENT ARGUMENTS	#############
#################################################

make.args <- function(tune.args.i, mod.fun.name, occs.vals, bg.vals, other.args = NULL) {
  
  out <- list()
  
  if(mod.fun.name == "maxent") {
    out$args <- c("noaddsamplestobackground", "noremoveDuplicates", "noautofeature")
    
    if(!grepl("L", tune.args.i$fc)) out$args <- c(out$args, "nolinear")
    if(!grepl("Q", tune.args.i$fc)) out$args <- c(out$args, "noquadratic")
    if(!grepl("H", tune.args.i$fc)) out$args <- c(out$args, "nohinge")
    if(!grepl("P", tune.args.i$fc)) out$args <- c(out$args, "noproduct")
    if(!grepl("T", tune.args.i$fc)) out$args <- c(out$args, "nothreshold")
    out$args <- c(out$args, paste0("betamultiplier=", tune.args.i$rm, sep=""))
    out$x <- as.data.frame(rbind(occs.vals, bg.vals))
    out$p <- c(rep(1, nrow(occs.vals)), rep(0, nrow(bg.vals)))
  }

  if(mod.fun.name == "maxnet") {
    out$data <- rbind(occs.vals, bg.vals)
    out$p <- c(rep(1, nrow(occs.vals)), rep(0, nrow(bg.vals)))
    out$f <- maxnet::maxnet.formula(out$p, out$data, classes = tolower(tune.args.i$fc))
    out$regmult <- tune.args.i$rm
  }
  
  # need logic for BRTs
  
  # add other args
  out <- c(out, other.args)
  
  return(out)
	
}
