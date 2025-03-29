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

removeNARecs <- function(x, type = c("occurrence", "background")) {
  NArecs <- which(rowSums(is.na(x)) > 0)
  if(length(NArecs) > 0) {
    msg(paste0("* Removed ", length(NArecs), " ", type, " records with NA predictor variable values."))
    x <- x[-NArecs,]
  }
  return(x)
}

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
