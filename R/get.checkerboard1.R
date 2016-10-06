#############################################################################
#########	MAKE CHECKERBOARD1 EVALUATION GROUPS	#############
#############################################################################

get.checkerboard1 <- function(occ, env, bg.coords, aggregation.factor){

  occ <- as.data.frame(occ)
  rownames(occ) <- 1:nrow(occ)
  bg.coords <- as.data.frame(bg.coords)
  rownames(bg.coords) <- 1:nrow(bg.coords)
  
	grid <- aggregate(env[[1]], fact=aggregation.factor[1])
	w <- gridSample(occ, grid, n=1e6, chess='white')
	b <- gridSample(occ, grid, n=1e6, chess='black')
	bgw <- gridSample(bg.coords, grid, n=1e6, chess='white')
	bgb <- gridSample(bg.coords, grid, n=1e6, chess='black')

	w$grp <- 1
	b$grp <- 2
	r <- rbind(w, b)
	occ.grp <- r[order(as.numeric(rownames(r))),]$grp

	bgw$grp <- 1
	bgb$grp <- 2
	bgr <- rbind(bgw, bgb)
	bg.grp <- bgr[order(as.numeric(rownames(bgr))),]$grp

	out <- list(occ.grp=occ.grp, bg.grp=bg.grp)
	return(out)
}
