#################################################################
#########	MAKE K-1 JACKKNIFE EVALUATION GROUPS	#############
#################################################################

get.jackknife <- function(occ, bg.coords) {
	occ.grp <- 1:nrow(occ)
	bg.grp <- rep(0, nrow(bg.coords))
	out <- list(occ.grp=occ.grp, bg.grp=bg.grp)
	return(out)
}
