#########################################################
#########	MAKE K RANDOM EVALUATION GROUPS	#############
#########################################################

get.randomkfold <- function(occ, bg.coords, kfolds){
	occ.grp <- kfold(occ, kfolds)
	bg.grp <- rep(0, nrow(bg.coords))
	out <- list(occ.grp=occ.grp, bg.grp=bg.grp)
	return(out)	
}
