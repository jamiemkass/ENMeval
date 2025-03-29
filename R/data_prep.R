# remove occurrence records that are duplicates in the same grid cell
removeOccDups <- function(occs, envs) {
  # extract cell number values from occs to determine duplicates
  occs.cellNo <- terra::extract(envs, occs, cells = TRUE, ID = FALSE)
  occs.dups <- duplicated(occs.cellNo[,"cell"])
  if(sum(occs.dups) > 0) {
    msg(paste0("* Removed ", sum(occs.dups), " occurrence localities that shared the same grid cell."))
    occs <- occs[!occs.dups,]
  }
  return(occs)
}

# remove records with NA environmental values
removeNARecs <- function(x, type = c("occurrence", "background")) {
  NArecs <- which(rowSums(is.na(x)) > 0)
  if(length(NArecs) > 0) {
    msg(paste0("* Removed ", length(NArecs), " ", type, " records with NA predictor variable values."))
    x <- x[-NArecs,]
  }
  return(x)
}

# remove occurrence and background record user.grps that were filtered
# from previous cleaning steps
cleanUsrGrp <- function(user.grp, occs, bg) {
  occs.len <- length(as.numeric(rownames(occs)))
  user.grp.occs.len <- length(user.grp$occs.grp)
  occs.diff <- user.grp.occs.len - occs.len 
  if(occs.diff > 0) {
    msg(paste0("Removing ", occs.diff, " occurrence records from user.grp."))
    user.grp$occs.grp <- user.grp$occs.grp[as.numeric(rownames(occs))]
  }
  bg.len <- length(as.numeric(rownames(bg)))
  user.grp.bg.len <- length(user.grp$bg.grp)
  bg.diff <- user.grp.bg.len - bg.len 
  if(bg.diff > 0) {
    msg(paste0("Removing ", bg.diff, " background records from user.grp."))
    user.grp$bg.grp <- user.grp$bg.grp[as.numeric(rownames(bg))]
  }
  return(user.grp)
}

# make a list of categorical variable levels
catLevs <- function(d, envs, categoricals) {
  # make categorical levels list
  cat.levs <- list()
  for(i in 1:length(categoricals)) {
    if(!is.null(envs)) {
      cat.levs[[i]] <- terra::levels(envs[[categoricals[i]]])[[1]][,2]
    }else{
      cat.levs[[i]] <- levels(d[, categoricals[i]])
    }  
  }
  return(cat.levs)
}

# converts categorical variable fields to factors and coerces them to numbers for maxent.jar
numFactors <- function(x, d, i) {
  x[, categoricals[i]] <- factor(as.numeric(x[, categoricals[i]]), levels = levels(d[, categoricals[i]]))
  return(x)
}

# converts categorical variable fields to factors
regFactors <- function(x, d, i) {
  x[, categoricals[i]] <- factor(x[, categoricals[i]], levels = levels(d[, categoricals[i]]))
  return(x)
}

  

  
  