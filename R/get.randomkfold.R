#########################################################
#########	MAKE K RANDOM EVALUATION GROUPS	#############
#########################################################

get.randomkfold <- function(occ, bg.coords, kfolds){
  occ <- as.data.frame(occ)
  bg.coords <- as.data.frame(bg.coords)
  occ.grp <- kfold(occ, kfolds)
	bg.grp <- rep(0, nrow(bg.coords))
	out <- list(occ.grp=occ.grp, bg.grp=bg.grp)
	return(out)	
}
