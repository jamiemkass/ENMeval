#############################################################################
#########	MAKE CHECKERBOARD1 EVALUATION GROUPS	#############
#############################################################################

get.checkerboard1 <- function(occ, env, bg.coords, aggregation.factor){

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

# FIX IF OCC POINTS FALL INTO A SINGLE BIN
	noccgrp <- length(unique(occ.grp))
	nbggrp <- length(unique(bg.grp))
	if(noccgrp < 2 ){
		message(paste("Warning: occurrence points fall in only", noccgrp, "bin"))
		bg.grp[ ! bg.grp %in% occ.grp] <- NA
		occ.grp <- as.numeric(as.factor(occ.grp))
		bg.grp <- as.numeric(as.factor(bg.grp))
		}

	out <- list(occ.grp=occ.grp, bg.grp=bg.grp)
	return(out)
}
