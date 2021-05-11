#' @title Methods to partition data for evaluation
#' 
#' @description \code{ENMeval} provides several ways to partition occurrence and background localities into bins for training and validation (or, evaluation and calibration). 
#' Users should carefully consider the objectives of their study and the influence of spatial bias when deciding on a method of data partitioning.
#' 
#' These functions are used internally to partition data during a call of \code{\link{ENMevaluate}} but can also be used independently to generate data partitions. 
#' For user-specified partitions, users can simply define groups of occurrence records and background points directly with \code{ENMevaluate}.
#' 
#' The \code{get.block} method partitions occurrence localities by finding the latitude and/or longitude lines that divide the occurrence localities into four groups of (insofar as possible) equal numbers.
#' The order and nature of the divisions can be controlled with the "orientation" parameter.
#' The default is "lat_lon", which divides first by a latitudinal line, then second by longitudinal lines.
#' This method is based on the spatial partitioning technique described in Radosavljevic & Anderson (2014), where the "lon_lon" option was used.
#' Background localities are assigned to each of the four groups based on their position with respect to these lines. 
#' While the \code{get.block} method results in (approximately) equal division of occurrence localities among four groups, the number of background localities (and, consequently, environmental and geographic space) in each group depends on the distribution of occurrence localities across the study area.
#' 
#' The \code{get.checkerboard1} and \code{get.checkerboard2} methods are variants of a checkerboard approach to partition occurrence localities. 
#' These methods use the \code{dismo::gridSample} function of the \pkg{dismo} package (Hijmans \emph{et al.} 2011) to partition records according to checkerboard grids across the study extent. 
#' The spatial grain of these grids is determined by resampling (or aggregating) the original environmental input grids based on the user-defined \code{aggregation factor} (e.g., an aggregation factor of 2 results in a checkerboard with grid cells four times as large in area as the original input grids). 
#' The \code{get.checkerboard1} method partitions data into two groups according to a single checkerboard pattern, and the \code{get.checkerboard2} method partitions data into four groups according to two nested checkerboard grids. 
#' In contrast to the \code{get.block} method, both the \code{get.checkerboard1} and \code{get.checkerboard2} methods subdivide geographic space equally but do not ensure a balanced number of occurrence localities in each group. 
#' The two \code{get.checkerboard} methods give warnings (and potentially errors) if zero points (occurrence or background) fall in any of the expected bins.
#' 
#' The \code{get.jackknife} method is a special case of \emph{k}-fold cross validation where the number of bins (\emph{k}) is equal to the number of occurrence localities (\emph{n}) in the dataset. 
#' It is suggested for occurrence datasets of relatively small sample size (generally < 25 localities) (Pearson \emph{et al.} 2007; Shcheglovitova and Anderson 2013).
#' 
#' The \code{get.randomkfold} method partitions occurrence localities randomly into a user-specified number of (\emph{k}) bins. 
#' This is equivalent to the method of \emph{k}-fold cross validation currently provided by Maxent. 
#' 
#' Users can also define custom partitions for occurrence and background data in the call to `ENMevaluate` with the "user.grp" parameter. 
#' 
#' @param occs matrix / data frame: longitude and latitude (in that order) of occurrence localities
#' @param bg matrix / data frame: longitude and latitude (in that order) of background localities
#' @param envs RasterStack: environmental predictor variables
#' @param orientation character vector: the order of spatial partitioning for the \code{get.block} method;
#' the first direction bisects the points into two groups, and the second direction bisects each of these further into two groups each, resulting in four groups; 
#' options are "lat_lon" (default), "lon_lat", "lon_lon", and "lat_lat"
#' @param aggregation.factor numeric vector: the aggregation scale for the \code{get.checkerboard1} and \code{get.checkerboard2} methods;
#' if a single number is given and \code{get.checkerboard2} partitioning method is used, the single value is used for both scales of aggregation
#' @param gridSampleN numeric: the number of points sampled from the input raster using gridSample() by the checkerboard partitioning functions
#' @param kfolds numeric: number of random \emph{k}-folds for \code{get.randomkfold} method
#' 
#' @return
#' A named list of two items:
#' \item{$occs.grp}{ A vector of bin designation for occurrence localities in the same order they were provided.}
#' \item{$bg.grp}{ A vector of bin designation for background localities in the same order they were provided.}
#' 
#' @note 
#' The \code{checkerboard1} and \code{checkerboard2} methods are designed to partition occurrence localities into two and four evaluation bins, respectively. 
#' They may give fewer bins, however, depending on where the occurrence localities fall with respect to the grid cells (e.g., all records happen to fall in the "black" squares). 
#' A warning is given if the number of bins is < 4 for the \code{checkerboard2} method, and an error is given if all localities fall into a single evaluation bin.
#' 
#' @references
#' Hijmans, R. J., Phillips, S., Leathwick, J. and Elith, J. (2011). dismo package for R. Available online at: \url{https://cran.r-project.org/package=dismo}.
#' 
#' Pearson, R. G., Raxworthy, C. J., Nakamura, M. and Peterson, A. T. (2007). Predicting species distributions from small numbers of occurrence records: a test case using cryptic geckos in Madagascar. \emph{Journal of Biogeography}, \bold{34}: 102-117. \url{https://doi.org/10.1111/j.1365-2699.2006.01594.x}
#' 
#' Radosavljevic, A., & Anderson, R. P. (2014). Making better Maxent models of species distributions: complexity, overfitting and evaluation. \emph{Journal of Biogeography}, \bold{41}: 629-643. \url{https://doi.org/10.1111/jbi.12227}
#' 
#' Shcheglovitova, M. and Anderson, R. P. (2013). Estimating optimal complexity for ecological niche models: a jackknife approach for species with small sample sizes. \emph{Ecological Modelling}, \bold{269}: 9-17. \url{https://doi.org/10.1016/j.ecolmodel.2013.08.011}
#' 
#' @author 
#' Robert Muscarella <bob.muscarella@gmail.com> and Jamie M. Kass <jkass@gc.cuny.edu>
#' 
#' @examples 
#' require(raster)
#' 
#' set.seed(1)
#' 
#' ### Create environmental extent (raster)
#' envs <- raster(matrix(nrow=25, ncol=25))
#' 
#' ### Create occurrence localities
#' set.seed(1)
#' nocc <- 25
#' xocc <- rnorm(nocc, sd=0.25) + 0.5
#' yocc <- runif(nocc, 0, 1)
#' occs <- as.data.frame(cbind(xocc, yocc))
#' 
#' ### Create background points
#' nbg <- 500
#' xbg <- runif(nbg, 0, 1)
#' ybg <- runif(nbg, 0, 1)
#' bg <- as.data.frame(cbind(xbg, ybg))
#' 
#' ### Plot points on environmental raster
#' plot(envs)
#' points(bg)
#' points(occs, pch=21, bg=2)
#' 
#' ### Block partitioning method (default orientation is "lat_lon"))
#' blk.latLon <- get.block(occs, bg)
#' plot(envs)
#' points(occs, pch=23, bg=blk.latLon$occs.grp)
#' plot(envs)
#' points(bg, pch=21, bg=blk.latLon$bg.grp)
#' # Can partition with other orientations
#' blk.latLat <- get.block(occs, bg, orientation = "lat_lat")
#' plot(envs)
#' points(occs, pch=23, bg=blk.latLat$occs.grp)
#' plot(envs)
#' points(bg, pch=21, bg=blk.latLat$bg.grp)
#' 
#' ### Checkerboard1 partitioning method with aggregation factor of 4
#' chk1.ag4 <- get.checkerboard1(occs, envs, bg, aggregation.factor = 4)
#' plot(envs)
#' points(occs, pch=23, bg=chk1.ag4$occs.grp)
#' plot(envs)
#' points(bg, pch=21, bg=chk1.ag4$bg.grp)
#' # Higher aggregation factors result in bigger checkerboard blocks
#' chk1.ag8 <- get.checkerboard1(occs, envs, bg, aggregation.factor = 8)
#' plot(envs)
#' points(occs, pch=23, bg=chk1.ag8$occs.grp)
#' plot(envs)
#' points(bg, pch=21, bg=chk1.ag8$bg.grp)
#' 
#' ### Checkerboard2 partitioning method with aggregation factors of 2, 2
#' chk2.ag2_2 <- get.checkerboard2(occs, envs, bg, c(2,2))
#' plot(envs)
#' points(occs, pch=23, bg=chk2.ag2_2$occs.grp)
#' plot(envs)
#' points(bg, pch=21, bg=chk2.ag2_2$bg.grp)
#' # Higher aggregation factors result in bigger checkerboard blocks,
#' # and can vary between hierarchical levels
#' chk2.ag4_6 <- get.checkerboard2(occs, envs, bg, c(4,6))
#' plot(envs)
#' points(occs, pch=23, bg=chk2.ag4_6$occs.grp)
#' plot(envs)
#' points(bg, pch=21, bg=chk2.ag4_6$bg.grp)
#' 
#' ### Random partitions with 4 folds
#' # Note that get.randomkkfold does not partition the background
#' krandom <- get.randomkfold(occs, bg, 4)
#' plot(envs)
#' points(occs, pch=23, bg=krandom$occs.grp)
#' plot(envs)
#' points(bg, pch=21, bg=krandom$bg.grp)
#' 
#' ### k-1 jackknife partitions
#' # Note that get.jackknife does not partition the background
#' jack <- get.jackknife(occs, bg)
#' plot(envs)
#' points(occs, pch=23, bg=rainbow(length(jack$occs.grp)))
#' plot(envs)
#' points(bg, pch=21, bg=jack$bg.grp)
#' 

#' @name partitions
NULL

#' @rdname partitions
#' 
#' @export

get.block <- function(occs, bg, orientation = "lat_lon"){
  if(!(orientation %in% c("lat_lon", "lon_lat", "lon_lon", "lat_lat"))) stop('Please enter orientation that is one of "lat_lon", "lon_lat", "lon_lon", or "lat_lat".')
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
  
  grpAB.d <- switch(orientation, lat_lon = 2, lon_lat = 1, lon_lon = 1, lat_lat = 2)
  grpA <- occs[order(occs[, grpAB.d]),][1:n1,]
  grpB <- occs[rev(order(occs[, grpAB.d])),][1:n2,]  
  
  grpNum.d <- switch(orientation, lat_lon = 1, lon_lat = 2, lon_lon = 1, lat_lat = 2)
  grp1 <- grpA[order(grpA[, grpNum.d]),][1:(n3),]
  grp2 <- grpA[!rownames(grpA ) %in% rownames(grp1),]
  grp3 <- grpB[order(grpB[, grpNum.d]),][1:(n4),]
  grp4 <- grpB[!rownames(grpB) %in% rownames(grp3),]
  
  # SPLIT BACKGROUND POINTS BASED ON SPATIAL GROUPS
  
  if(orientation == "lat_lon") {
    bvert <- mean(c(max(grp1[, 1]), min(grp2[, 1])))
    tvert <- mean(c(max(grp3[, 1]), min(grp4[, 1])))
    horz <- mean(c(max(grpA[, 2]), min(grpB[, 2])))  
    bggrp1 <- bg[bg[, 2] <= horz & bg[, 1] < bvert,]
    bggrp2 <- bg[bg[, 2] <= horz & bg[, 1] >= bvert,]
    bggrp3 <- bg[bg[, 2] > horz & bg[, 1] <= tvert,]
    bggrp4 <- bg[bg[, 2] > horz & bg[, 1] > tvert,]
  }else if(orientation == "lon_lat") {
    bhorz <- mean(c(max(grp1[, 2]), min(grp2[, 2])))
    thorz <- mean(c(max(grp3[, 2]), min(grp4[, 2])))
    vert <- mean(c(max(grpA[, 1]), min(grpB[, 1])))  
    bggrp1 <- bg[bg[, 1] <= vert & bg[, 2] < bhorz,]
    bggrp2 <- bg[bg[, 1] <= vert & bg[, 2] >= bhorz,]
    bggrp3 <- bg[bg[, 1] > vert & bg[, 2] <= thorz,]
    bggrp4 <- bg[bg[, 1] > vert & bg[, 2] > thorz,]
  }else if(orientation == "lon_lon") {
    bvert <- mean(c(max(grp1[, 1]), min(grp2[, 1])))
    tvert <- mean(c(max(grp3[, 1]), min(grp4[, 1])))
    vert <- mean(c(max(grpA[, 1]), min(grpB[, 1])))
    bggrp1 <- bg[bg[, 1] <= bvert,]
    bggrp2 <- bg[bg[, 1] > bvert & bg[, 1] <= vert,]
    bggrp3 <- bg[bg[, 1] >= vert & bg[, 1] < tvert,]
    bggrp4 <- bg[bg[, 1] >= tvert,]
  }else if(orientation == "lat_lat"){
    bhorz <- mean(c(max(grp1[, 2]), min(grp2[, 2])))
    thorz <- mean(c(max(grp3[, 2]), min(grp4[, 2])))
    horz <- mean(c(max(grpA[, 2]), min(grpB[, 2])))
    bggrp1 <- bg[bg[, 2] <= bhorz,]
    bggrp2 <- bg[bg[, 2] > bhorz & bg[, 2] <= horz,]
    bggrp3 <- bg[bg[, 2] >= horz & bg[, 2] < thorz,]
    bggrp4 <- bg[bg[, 2] >= thorz,]
  }
  
  r <- data.frame()
  if (nrow(grp1) > 0) grp1$grp <- 1; r <- rbind(r, grp1)
  if (nrow(grp2) > 0) grp2$grp <- 2; r <- rbind(r, grp2)
  if (nrow(grp3) > 0) grp3$grp <- 3; r <- rbind(r, grp3)
  if (nrow(grp4) > 0) grp4$grp <- 4; r <- rbind(r, grp4)
  occs.grp <- r[order(as.numeric(rownames(r))),]$grp
  
  bgr <- data.frame()
  if (nrow(bggrp1) > 0) bggrp1$grp <- 1; bgr <- rbind(bgr, bggrp1)
  if (nrow(bggrp2) > 0) bggrp2$grp <- 2; bgr <- rbind(bgr, bggrp2)
  if (nrow(bggrp3) > 0) bggrp3$grp <- 3; bgr <- rbind(bgr, bggrp3)
  if (nrow(bggrp4) > 0) bggrp4$grp <- 4; bgr <- rbind(bgr, bggrp4)
  bg.grp <- bgr[order(as.numeric(rownames(bgr))),]$grp
  
  out <- list(occs.grp=occs.grp, bg.grp=bg.grp)
  return(out)
}

#' @rdname partitions
#' 
#' @export

get.checkerboard1 <- function(occs, envs, bg, aggregation.factor, gridSampleN = 10000){
  if(is.null(envs)) stop("Cannot use checkerboard partitioning if envs is NULL.")
  occs <- as.data.frame(occs)
  rownames(occs) <- 1:nrow(occs)
  bg <- as.data.frame(bg)
  rownames(bg) <- 1:nrow(bg)
  
  grid <- raster::aggregate(envs[[1]], fact=aggregation.factor[1])
  w <- dismo::gridSample(occs, grid, n=gridSampleN, chess='white')
  b <- dismo::gridSample(occs, grid, n=gridSampleN, chess='black')
  bgw <- dismo::gridSample(bg, grid, n=gridSampleN, chess='white')
  bgb <- dismo::gridSample(bg, grid, n=gridSampleN, chess='black')
  
  if(nrow(w) > 0) { w$grp <- 1 }
  if(nrow(b) > 0) { b$grp <- 2 }
  r <- rbind(w, b)
  occs.grp <- r[order(as.numeric(rownames(r))),]$grp
  
  if(nrow(bgw) > 0) { bgw$grp <- 1 }
  if(nrow(bgb) > 0) { bgb$grp <- 2 }
  bgr <- rbind(bgw, bgb)
  bg.grp <- bgr[order(as.numeric(rownames(bgr))),]$grp
  
  # PATCH IF occs OR BG POINTS FALL INTO A SINGLE BIN
  noccgrp <- length(unique(occs.grp))
  nbggrp <- length(unique(bg.grp))
  if(noccgrp < 2 ){
    message(paste("Warning: occurrence points fall in only", noccgrp, "bin"))
    bg.grp[ ! bg.grp %in% occs.grp] <- NA
    occs.grp <- as.numeric(as.factor(occs.grp))
    bg.grp <- as.numeric(as.factor(bg.grp))
  }
  
  if(length(unique(bg.grp[!is.na(bg.grp)])) != noccgrp) {
    stop("Error: occurrence records but no background points fall in 1 or more evaluation bin(s)")
  }
  
  out <- list(occs.grp=occs.grp, bg.grp=bg.grp)
  return(out)
}

#' @rdname partitions
#' 
#' @export

get.checkerboard2 <- function(occs, envs, bg, aggregation.factor, gridSampleN = 10000){
  if(is.null(envs)) stop("Cannot use checkerboard partitioning if envs is NULL.")
  occs <- as.data.frame(occs)
  rownames(occs) <- 1:nrow(occs)
  bg <- as.data.frame(bg)
  rownames(bg) <- 1:nrow(bg)
  
  if (length(aggregation.factor) == 1) aggregation.factor <- rep(aggregation.factor, 2)
  grid <- raster::aggregate(envs[[1]], fact=aggregation.factor[1])
  grid2 <- raster::aggregate(grid, aggregation.factor[2])
  w <- dismo::gridSample(occs, grid, n=gridSampleN, chess='white')
  b <- dismo::gridSample(occs, grid, n=gridSampleN, chess='black')
  ww <- dismo::gridSample(w, grid2, n=gridSampleN, chess='white')
  wb <- dismo::gridSample(w, grid2, n=gridSampleN, chess='black')
  bw <- dismo::gridSample(b, grid2, n=gridSampleN, chess='white')
  bb <- dismo::gridSample(b, grid2, n=gridSampleN, chess='black')
  bgw <- dismo::gridSample(bg, grid, n=gridSampleN, chess='white')
  bgb <- dismo::gridSample(bg, grid, n=gridSampleN, chess='black')
  bgww <- dismo::gridSample(bgw, grid2, n=gridSampleN, chess='white')
  bgwb <- dismo::gridSample(bgw, grid2, n=gridSampleN, chess='black')
  bgbw <- dismo::gridSample(bgb, grid2, n=gridSampleN, chess='white')
  bgbb <- dismo::gridSample(bgb, grid2, n=gridSampleN, chess='black')
  
  r <- data.frame()
  if (nrow(ww) > 0) ww$grp <- 1; r <- rbind(r, ww)
  if (nrow(wb) > 0) wb$grp <- 2; r <- rbind(r, wb)
  if (nrow(bw) > 0) bw$grp <- 3; r <- rbind(r, bw)
  if (nrow(bb) > 0) bb$grp <- 4; r <- rbind(r, bb)
  occs.grp <- r[order(as.numeric(rownames(r))),]$grp
  
  bgr <- data.frame()
  if (nrow(bgww) > 0) bgww$grp <- 1; bgr <- rbind(bgr, bgww)
  if (nrow(bgwb) > 0) bgwb$grp <- 2; bgr <- rbind(bgr, bgwb)
  if (nrow(bgbw) > 0) bgbw$grp <- 3; bgr <- rbind(bgr, bgbw)
  if (nrow(bgbb) > 0) bgbb$grp <- 4; bgr <- rbind(bgr, bgbb)
  bg.grp <- bgr[order(as.numeric(rownames(bgr))),]$grp
  
  # PATCH IF occs OR BG POINTS FALL INTO FEWER THAN FOUR BINS
  noccgrp <- length(unique(occs.grp))
  nbggrp <- length(unique(bg.grp))
  if(noccgrp < 4 ){
    message(paste("Warning: occurrence points fall in only", noccgrp, "bins"))
    bg.grp[ ! bg.grp %in% occs.grp] <- NA
    occs.grp <- as.numeric(as.factor(occs.grp))
    bg.grp <- as.numeric(as.factor(bg.grp))
  }
  
  if(length(unique(bg.grp[!is.na(bg.grp)])) != noccgrp) {
    stop("Error: occurrence records but no background points fall in 1 or more evaluation bin(s)")
  }
  
  out <- list(occs.grp=occs.grp, bg.grp=bg.grp)
  return(out)
}

#' @rdname partitions
#' 
#' @export

get.jackknife <- function(occs, bg) {
  occs <- as.data.frame(occs)
  rownames(occs) <- 1:nrow(occs)
  bg <- as.data.frame(bg)
  rownames(bg) <- 1:nrow(bg)
  occs.grp <- 1:nrow(occs)
  bg.grp <- rep(0, nrow(bg))
  out <- list(occs.grp=occs.grp, bg.grp=bg.grp)
  return(out)
}

#' @rdname partitions
#' 
#' @export

get.randomkfold <- function(occs, bg, kfolds){
  occs <- as.data.frame(occs)
  rownames(occs) <- 1:nrow(occs)
  bg <- as.data.frame(bg)
  rownames(bg) <- 1:nrow(bg)
  occs.grp <- dismo::kfold(occs, kfolds)
  bg.grp <- rep(0, nrow(bg))
  out <- list(occs.grp=occs.grp, bg.grp=bg.grp)
  return(out)	
}
