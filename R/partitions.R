#' @export

get.jackknife <- function(occs, bg) {
  occs <- as.data.frame(occs)
  rownames(occs) <- 1:nrow(occs)
  bg <- as.data.frame(bg)
  rownames(bg) <- 1:nrow(bg)
  occs.folds <- 1:nrow(occs)
  bg.folds <- rep(0, nrow(bg))
  out <- list(occs.folds=occs.folds, bg.folds=bg.folds)
  return(out)
}

#' @export

get.randomkfold <- function(occs, bg, kfolds){
  occs <- as.data.frame(occs)
  rownames(occs) <- 1:nrow(occs)
  bg <- as.data.frame(bg)
  rownames(bg) <- 1:nrow(bg)
  occs.folds <- kfold(occs, kfolds)
  bg.folds <- rep(0, nrow(bg))
  out <- list(occs.folds=occs.folds, bg.folds=bg.folds)
  return(out)	
}

#' @export

get.block <- function(occs, bg){
  occs <- as.data.frame(occs)
  rownames(occs) <- 1:nrow(occs)
  bg <- as.data.frame(bg)
  rownames(bg) <- 1:nrow(bg)
  
  # SPLIT occs POINTS INTO FOUR SPATIAL GROUPS
  noccs <- nrow(occs)
  n1 <- ceiling(nrow(occs)/2)
  n2 <- floor(nrow(occs)/2)
  n3 <- ceiling(n1/2)
  n4 <- ceiling(n2/2)
  grpA <- occs[order(occs[, 2]),][1:n1,]
  grpB <- occs[rev(order(occs[, 2])),][1:n2,]
  grp1 <- grpA[order(grpA[, 1]),][1:(n3),]
  grp2 <- grpA[!rownames(grpA)%in%rownames(grp1),]
  grp3 <- grpB[order(grpB[, 1]),][1:(n4),]
  grp4 <- grpB[!rownames(grpB)%in%rownames(grp3),]
  
  # SPLIT BACKGROUND POINTS BASED ON SPATIAL GROUPS
  bvert <- mean(c(max(grp1[, 1]), min(grp2[, 1])))
  tvert <- mean(c(max(grp3[, 1]), min(grp4[, 1])))
  horz <- mean(c(max(grpA[, 2]), min(grpB[, 2])))
  bggrp1 <- bg[bg[, 2] <= horz & bg[, 1]<bvert,]
  bggrp2 <- bg[bg[, 2] < horz & bg[, 1]>=bvert,]
  bggrp3 <- bg[bg[, 2] > horz & bg[, 1]<=tvert,]
  bggrp4 <- bg[bg[, 2] >= horz & bg[, 1]>tvert,]
  
  r <- data.frame()
  if (nrow(grp1) > 0) grp1$grp <- 1; r <- rbind(r, grp1)
  if (nrow(grp2) > 0) grp2$grp <- 2; r <- rbind(r, grp2)
  if (nrow(grp3) > 0) grp3$grp <- 3; r <- rbind(r, grp3)
  if (nrow(grp4) > 0) grp4$grp <- 4; r <- rbind(r, grp4)
  occs.folds <- r[order(as.numeric(rownames(r))),]$grp
  
  bgr <- data.frame()
  if (nrow(bggrp1) > 0) bggrp1$grp <- 1; bgr <- rbind(bgr, bggrp1)
  if (nrow(bggrp2) > 0) bggrp2$grp <- 2; bgr <- rbind(bgr, bggrp2)
  if (nrow(bggrp3) > 0) bggrp3$grp <- 3; bgr <- rbind(bgr, bggrp3)
  if (nrow(bggrp4) > 0) bggrp4$grp <- 4; bgr <- rbind(bgr, bggrp4)
  bg.folds <- bgr[order(as.numeric(rownames(bgr))),]$grp
  
  out <- list(occs.folds=occs.folds, bg.folds=bg.folds)
  return(out)
}

#' @export

get.checkerboard1 <- function(occs, env, bg, aggregation.factor){
  
  occs <- as.data.frame(occs)
  rownames(occs) <- 1:nrow(occs)
  bg <- as.data.frame(bg)
  rownames(bg) <- 1:nrow(bg)
  
  grid <- aggregate(env[[1]], fact=aggregation.factor[1])
  w <- gridSample(occs, grid, n=1e6, chess='white')
  b <- gridSample(occs, grid, n=1e6, chess='black')
  bgw <- gridSample(bg, grid, n=1e6, chess='white')
  bgb <- gridSample(bg, grid, n=1e6, chess='black')
  
  if(nrow(w) > 0) { w$grp <- 1 }
  if(nrow(b) > 0) { b$grp <- 2 }
  r <- rbind(w, b)
  occs.folds <- r[order(as.numeric(rownames(r))),]$grp
  
  if(nrow(bgw) > 0) { bgw$grp <- 1 }
  if(nrow(bgb) > 0) { bgb$grp <- 2 }
  bgr <- rbind(bgw, bgb)
  bg.folds <- bgr[order(as.numeric(rownames(bgr))),]$grp
  
  # PATCH IF occs OR BG POINTS FALL INTO A SINGLE BIN
  noccgrp <- length(unique(occs.folds))
  nbggrp <- length(unique(bg.folds))
  if(noccgrp < 2 ){
    message(paste("Warning: occurrence points fall in only", noccgrp, "bin"))
    bg.folds[ ! bg.folds %in% occs.folds] <- NA
    occs.folds <- as.numeric(as.factor(occs.folds))
    bg.folds <- as.numeric(as.factor(bg.folds))
  }
  
  if(length(unique(bg.folds[!is.na(bg.folds)])) != noccgrp) {
    message("Error: occurrence records but no background points fall in 1 or more evaluation bin(s)")
    stop()
  }
  
  out <- list(occs.folds=occs.folds, bg.folds=bg.folds)
  return(out)
}

#' @export

get.checkerboard2 <- function(occs, env, bg, aggregation.factor){
  
  occs <- as.data.frame(occs)
  rownames(occs) <- 1:nrow(occs)
  bg <- as.data.frame(bg)
  rownames(bg) <- 1:nrow(bg)
  
  if (length(aggregation.factor) == 1) aggregation.factor <- rep(aggregation.factor, 2)
  grid <- aggregate(env[[1]], fact=aggregation.factor[1])
  grid2 <- aggregate(grid, aggregation.factor[2])
  w <- gridSample(occs, grid, n=1e4, chess='white')
  b <- gridSample(occs, grid, n=1e4, chess='black')
  ww <- gridSample(w, grid2, n=1e4, chess='white')
  wb <- gridSample(w, grid2, n=1e4, chess='black')
  bw <- gridSample(b, grid2, n=1e4, chess='white')
  bb <- gridSample(b, grid2, n=1e4, chess='black')
  bgw <- gridSample(bg, grid, n=1e4, chess='white')
  bgb <- gridSample(bg, grid, n=1e4, chess='black')
  bgww <- gridSample(bgw, grid2, n=1e4, chess='white')
  bgwb <- gridSample(bgw, grid2, n=1e4, chess='black')
  bgbw <- gridSample(bgb, grid2, n=1e4, chess='white')
  bgbb <- gridSample(bgb, grid2, n=1e4, chess='black')
  
  r <- data.frame()
  if (nrow(ww) > 0) ww$grp <- 1; r <- rbind(r, ww)
  if (nrow(wb) > 0) wb$grp <- 2; r <- rbind(r, wb)
  if (nrow(bw) > 0) bw$grp <- 3; r <- rbind(r, bw)
  if (nrow(bb) > 0) bb$grp <- 4; r <- rbind(r, bb)
  occs.folds <- r[order(as.numeric(rownames(r))),]$grp
  
  bgr <- data.frame()
  if (nrow(bgww) > 0) bgww$grp <- 1; bgr <- rbind(bgr, bgww)
  if (nrow(bgwb) > 0) bgwb$grp <- 2; bgr <- rbind(bgr, bgwb)
  if (nrow(bgbw) > 0) bgbw$grp <- 3; bgr <- rbind(bgr, bgbw)
  if (nrow(bgbb) > 0) bgbb$grp <- 4; bgr <- rbind(bgr, bgbb)
  bg.folds <- bgr[order(as.numeric(rownames(bgr))),]$grp
  
  # PATCH IF occs OR BG POINTS FALL INTO FEWER THAN FOUR BINS
  noccgrp <- length(unique(occs.folds))
  nbggrp <- length(unique(bg.folds))
  if(noccgrp < 4 ){
    message(paste("Warning: occurrence points fall in only", noccgrp, "bins"))
    bg.folds[ ! bg.folds %in% occs.folds] <- NA
    occs.folds <- as.numeric(as.factor(occs.folds))
    bg.folds <- as.numeric(as.factor(bg.folds))
  }
  
  if(length(unique(bg.folds[!is.na(bg.folds)])) != noccgrp) {
    message("Error: occurrence records but no background points fall in 1 or more evaluation bin(s)")
    stop()
  }
  
  out <- list(occs.folds=occs.folds, bg.folds=bg.folds)
  return(out)
}