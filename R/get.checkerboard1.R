#############################################################################
#########	MAKE CHECKERBOARD1 EVALUATION GROUPS	#############
#############################################################################

get.checkerboard1 <- function(occ, env, bg.coords, aggregation.factor){

	grid <- aggregate(env[[1]], fact=aggregation.factor[1])
	w <- gridSample(occ, grid, n=1e6, chess='white')
	b <- gridSample(occ, grid, n=1e6, chess='black')
	bgw <- gridSample(bg.coords, grid, n=1e6, chess='white')
	bgb <- gridSample(bg.coords, grid, n=1e6, chess='black')

	if(nrow(w) > 0) { w$grp <- 1 }
	if(nrow(b) > 0) { b$grp <- 2 }
	r <- rbind(w, b)
	occ.grp <- r[order(as.numeric(rownames(r))),]$grp

	if(nrow(bgw) > 0) { bgw$grp <- 1 }
	if(nrow(bgb) > 0) { bgb$grp <- 2 }
	bgr <- rbind(bgw, bgb)
	bg.grp <- bgr[order(as.numeric(rownames(bgr))),]$grp

# PATCH IF OCC OR BG POINTS FALL INTO A SINGLE BIN
	noccgrp <- length(unique(occ.grp))
	nbggrp <- length(unique(bg.grp))
	if(noccgrp < 2 ){
		message(paste("Warning: occurrence points fall in only", noccgrp, "bin"))
		bg.grp[ ! bg.grp %in% occ.grp] <- NA
		occ.grp <- as.numeric(as.factor(occ.grp))
		bg.grp <- as.numeric(as.factor(bg.grp))
		}

	if(length(unique(bg.grp[!is.na(bg.grp)])) != noccgrp) {
		message(paste("Error: occurrence records but no background points fall in 1 or more evaluation bin(s)"))
		stop()
		}

	out <- list(occ.grp=occ.grp, bg.grp=bg.grp)
	return(out)
}
