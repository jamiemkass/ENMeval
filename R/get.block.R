#########################################################
#########	MAKE BLOCK EVALUATION GROUPS	#############
#########################################################

get.block <- function(occ, bg.coords){
	
  occ <- as.data.frame(occ)
  rownames(occ) <- 1:nrow(occ)
  bg.coords <- as.data.frame(bg.coords)
  rownames(bg.coords) <- 1:nrow(bg.coords)

  # SPLIT OCC POINTS INTO FOUR SPATIAL GROUPS
	noccs <- nrow(occ)
	n1 <- ceiling(nrow(occ)/2)
	n2 <- floor(nrow(occ)/2)
	n3 <- ceiling(n1/2)
	n4 <- ceiling(n2/2)
	grpA <- occ[order(occ[, 2]),][1:n1,]
	grpB <- occ[rev(order(occ[, 2])),][1:n2,]
	grp1 <- grpA[order(grpA[, 1]),][1:(n3),]
	grp2 <- grpA[!rownames(grpA)%in%rownames(grp1),]
	grp3 <- grpB[order(grpB[, 1]),][1:(n4),]
	grp4 <- grpB[!rownames(grpB)%in%rownames(grp3),]

	# SPLIT BACKGROUND POINTS BASED ON SPATIAL GROUPS
	bvert <- mean(c(max(grp1[, 1]), min(grp2[, 1])))
	tvert <- mean(c(max(grp3[, 1]), min(grp4[, 1])))
	horz <- mean(c(max(grpA[, 2]), min(grpB[, 2])))
	bggrp1 <- bg.coords[bg.coords[, 2] <= horz & bg.coords[, 1] < bvert, ]
	bggrp2 <- bg.coords[bg.coords[, 2] <= horz & bg.coords[, 1] >= bvert, ]
	bggrp3 <- bg.coords[bg.coords[, 2] > horz & bg.coords[, 1] <= tvert, ]
	bggrp4 <- bg.coords[bg.coords[, 2] > horz & bg.coords[, 1] > tvert, ]

	r <- data.frame()
	if (nrow(grp1) > 0) grp1$grp <- 1; r <- rbind(r, grp1)
	if (nrow(grp2) > 0) grp2$grp <- 2; r <- rbind(r, grp2)
	if (nrow(grp3) > 0) grp3$grp <- 3; r <- rbind(r, grp3)
	if (nrow(grp4) > 0) grp4$grp <- 4; r <- rbind(r, grp4)
	occ.grp <- r[order(as.numeric(rownames(r))),]$grp

	bgr <- data.frame()
	if (nrow(bggrp1) > 0) bggrp1$grp <- 1; bgr <- rbind(bgr, bggrp1)
	if (nrow(bggrp2) > 0) bggrp2$grp <- 2; bgr <- rbind(bgr, bggrp2)
	if (nrow(bggrp3) > 0) bggrp3$grp <- 3; bgr <- rbind(bgr, bggrp3)
	if (nrow(bggrp4) > 0) bggrp4$grp <- 4; bgr <- rbind(bgr, bggrp4)
	bg.grp <- bgr[order(as.numeric(rownames(bgr))),]$grp

	out <- list(occ.grp=occ.grp, bg.grp=bg.grp)
	return(out)
}
