#############################################################################
#########	MAKE CHECKERBOARD 2 EVALUATION GROUPS	#############
#############################################################################

get.checkerboard2 <- function(occ, env, bg.coords, aggregation.factor){
	if (length(aggregation.factor) == 1) aggregation.factor <- rep(aggregation.factor, 2)
	grid <- aggregate(env[[1]], fact=aggregation.factor[1])
	grid2 <- aggregate(grid, aggregation.factor[2])
	w <- gridSample(occ, grid, n=1e4, chess='white')
	b <- gridSample(occ, grid, n=1e4, chess='black')
	ww <- gridSample(w, grid2, n=1e4, chess='white')
	wb <- gridSample(w, grid2, n=1e4, chess='black')
	bw <- gridSample(b, grid2, n=1e4, chess='white')
	bb <- gridSample(b, grid2, n=1e4, chess='black')
	bgw <- gridSample(bg.coords, grid, n=1e4, chess='white')
	bgb <- gridSample(bg.coords, grid, n=1e4, chess='black')
	bgww <- gridSample(bgw, grid2, n=1e4, chess='white')
	bgwb <- gridSample(bgw, grid2, n=1e4, chess='black')
	bgbw <- gridSample(bgb, grid2, n=1e4, chess='white')
	bgbb <- gridSample(bgb, grid2, n=1e4, chess='black')

	r <- data.frame()
	if (nrow(ww) > 0) ww$grp <- 1; r <- rbind(r, ww)
	if (nrow(wb) > 0) wb$grp <- 2; r <- rbind(r, wb)
	if (nrow(bw) > 0) bw$grp <- 3; r <- rbind(r, bw)
	if (nrow(bb) > 0) bb$grp <- 4; r <- rbind(r, bb)
	occ.grp <- r[order(as.numeric(rownames(r))),]$grp

	bgr <- data.frame()
	if (nrow(bgww) > 0) bgww$grp <- 1; bgr <- rbind(bgr, bgww)
	if (nrow(bgwb) > 0) bgwb$grp <- 2; bgr <- rbind(bgr, bgwb)
	if (nrow(bgbw) > 0) bgbw$grp <- 3; bgr <- rbind(bgr, bgbw)
	if (nrow(bgbb) > 0) bgbb$grp <- 4; bgr <- rbind(bgr, bgbb)
	bg.grp <- bgr[order(as.numeric(rownames(bgr))),]$grp

# PATCH IF OCC OR BG POINTS FALL INTO FEWER THAN FOUR BINS
	noccgrp <- length(unique(occ.grp))
	nbggrp <- length(unique(bg.grp))
	if(noccgrp < 4 ){
		message(paste("Warning: occurrence points fall in only", noccgrp, "bins"))
		bg.grp[ ! bg.grp %in% occ.grp] <- NA
		occ.grp <- as.numeric(as.factor(occ.grp))
		bg.grp <- as.numeric(as.factor(bg.grp))
		}

	if(length(unique(bg.grp[!is.na(bg.grp)])) != noccgrp) {
		message("Error: occurrence records but no background points fall in 1 or more evaluation bin(s)")
		stop()
		}

	out <- list(occ.grp=occ.grp, bg.grp=bg.grp)
	return(out)
}
