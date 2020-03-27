#' Methods to partition data for evaluation
#' 
#' \code{ENMeval} provides several ways to partition occurrence and background localities into bins for training and testing (or, evaluation and calibration).  Users should carefully consider the objectives of their study and the influence of spatial bias when deciding on a method of data partitioning.
#' 
#' These functions are used internally to partition data during a call of \code{\link{ENMevaluate}} but can also be used independently to generate data partitions. The old \code{get.user} function from previous releases of \code{ENMeval} has been removed because users can simply define groups of occurrence records and background points directly in the call of the \code{ENMevaluate} function.
#' 
#' The \code{get.block} method partitions occurrence localities by finding the latitude and longitude that divide the occurrence localities into four groups of (insofar as possible) equal numbers.  Background localities are assigned to each of the four groups based on their position with respect to these lines.  While the \code{get.block} method results in (approximately) equal division of occurrence localities among four groups, the number of background localities (and, consequently, environmental and geographic space) in each group depends on the distribution of occurrence localities across the study area.
#' 
#' The \code{get.checkerboard1} and \code{get.checkerboard2} methods are variants of a checkerboard approach to partition occurrence localities.  These methods use the \code{dismo::gridSample} function of the \pkg{dismo} package (Hijmans \emph{et al.} 2011) to partition records according to checkerboard grids across the study extent.  The spatial grain of these grids is determined by resampling (or aggregating) the original environmental input grids based on the user-defined \code{aggregation factor} (e.g., an aggregation factor of 2 results in a checkerboard with grid cells four times as large in area as the original input grids).  The \code{get.checkerboard1} method partitions data into two groups according to a single checkerboard pattern, and the \code{get.checkerboard2} method partitions data into four groups according to two nested checkerboard grids.  In contrast to the \code{get.block} method, both the \code{get.checkerboard1} and \code{get.checkerboard2} methods subdivide geographic space equally but do not ensure a balanced number of occurrence localities in each group.  The two \code{get.checkerboard} methods give warnings (and potentially errors) if zero points (occurrence or background) fall in any of the expected bins.
#' 
#' The \code{get.jackknife} method is a special case of \emph{k}-fold cross validation where the number of bins (\emph{k}) is equal to the number of occurrence localities (\emph{n}) in the dataset.  It is suggested for datasets of relatively small sample size (generally < 25 localities) (Pearson \emph{et al.} 2007; Shcheglovitova and Anderson 2013).
#' 
#' The \code{get.randomkfold} method partitions occurrence localities randomly into a user-specified number of (\emph{k}) bins.   This is equivalent to the method of \emph{k}-fold cross valiation currently provided by Maxent. 
#' 
#' Users can also define evaluation bins \emph{a priori} directly in the call to `ENMevaluate`.  With this method, occurrence and background localities, as well as evaluation bin designation for each locality, are supplied by the user.
#' 
#' @param occs Two-column matrix or data.frame of longitude and latitude (in that order) of occurrence localities.
#' @param bg.coords Two-column matrix or data.frame of longitude and latitude (in that order) of background localities.
#' @param env RasterStack of environmental predictor variables.
#' @param aggregation.factor A vector or list of 1 or 2 numbers giving the scale for aggregation used for the \code{get.checkerboard1} and \code{get.checkerboard2} methods.  If a single number is given and \code{get.checkerboard2} partitioning method is used, the single value is used for both scales of aggregation.
#' @param kfolds Number of random \emph{k}-folds for \code{get.randomkfold} method.
# #' @param occ.grp Vector of user-defined bins for occurrence localities for \code{get.user} method.
# #' @param bg.grp Vector of user-defined bins for background localities for \code{get.user} method.
#' 
#' @return
#' A named list of two items:
#' \item{$occ.grp}{ A vector of bin designation for occurrence localities in the same order they were provided.}
#' \item{$bg.grp}{ A vector of bin designation for background localities in the same order they were provided.}
#' 
#' @note 
#' The \code{checkerboard1} and \code{checkerboard2} methods are designed to partition occurrence localities into two and four evaluation bins, respectively.  They may give fewer bins, however, depending on where the occurrence localities fall with respect to the grid cells (e.g., all records happen to fall in the "black" squares).  A warning is given if the number of bins is < 4 for the \code{checkerboard2} method, and an error is given if all localities fall into a single evaluation bin.
#' 
#' @references
#' Hijmans, R. J., Phillips, S., Leathwick, J. and Elith, J. 2011. dismo package for R. Available online at: \url{https://cran.r-project.org/package=dismo}.
#'
#' Pearson, R. G., Raxworthy, C. J., Nakamura, M. and Peterson, A. T. 2007. Predicting species distributions from small numbers of occurrence records: a test case using cryptic geckos in Madagascar. \emph{Journal of Biogeography}, \bold{34}: 102-117.
#' 
#' Shcheglovitova, M. and Anderson, R. P. (2013) Estimating optimal complexity for ecological niche models: a jackknife approach for species with small sample sizes. \emph{Ecological Modelling}, \bold{269}: 9-17.
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
#' env <- raster(matrix(nrow=25, ncol=25))
#' 
#' ### Create presence localities
#' set.seed(1)
#' nocc <- 25
#' xocc <- rnorm(nocc, sd=0.25) + 0.5
#' yocc <- runif(nocc, 0, 1)
#' occ.pts <- as.data.frame(cbind(xocc, yocc))
#' 
#' ### Create background points
#' nbg <- 500
#' xbg <- runif(nbg, 0, 1)
#' ybg <- runif(nbg, 0, 1)
#' bg.pts <- as.data.frame(cbind(xbg, ybg))
#' 
#' ### Show points
#' plot(env)
#' points(bg.pts)
#' points(occ.pts, pch=21, bg=2)
#' 
#' ### Block partitioning method
#' blk.pts <- get.block(occ.pts, bg.pts)
#' plot(env)
#' points(occ.pts, pch=23, bg=blk.pts$occ.grp)
#' plot(env)
#' points(bg.pts, pch=21, bg=blk.pts$bg.grp)
#' 
#' ### Checkerboard1 partitioning method
#' chk1.pts <- get.checkerboard1(occ.pts, env, bg.pts, 4)
#' plot(env)
#' points(occ.pts, pch=23, bg=chk1.pts$occ.grp)
#' plot(env)
#' points(bg.pts, pch=21, bg=chk1.pts$bg.grp)
#' 
#' ### Checkerboard2 partitioning method
#' chk2.pts <- get.checkerboard2(occ.pts, env, bg.pts, c(2,2))
#' plot(env)
#' points(occ.pts, pch=23, bg=chk2.pts$occ.grp)
#' plot(env)
#' points(bg.pts, pch=21, bg=chk2.pts$bg.grp)
#' 
#' ### Random k-fold partitions
#' # Note that k random does not partition the background
#' krandom.pts <- get.randomkfold(occ.pts, bg.pts, 4)
#' plot(env)
#' points(occ.pts, pch=23, bg=krandom.pts$occ.grp)
#' plot(env)
#' points(bg.pts, pch=21, bg=krandom.pts$bg.grp)
#' 
#' ### k-1 jackknife partitions
#' # Note background is not partitioned
#' jack.pts <- get.jackknife(occ.pts, bg.pts)
#' plot(env)
#' points(occ.pts, pch=23, bg=rainbow(length(jack.pts$occ.grp)))
#' plot(env)
#' points(bg.pts, pch=21, bg=jack.pts$bg.grp)
#' 
# ## User-defined partitions
# #' # Note background is not partitioned
# #' occ.grp <- c(rep(1, 10), rep(2, 5), rep(3, 10))
# #' bg.grp <- c(rep(1, 200), rep(2, 100), rep(3, 200))
# #' user.pts <- get.user(occ.grp, bg.grp)
# #' plot(env)
# #' points(occ.pts, pch=23, bg=user.pts$occ.grp)
# #' plot(env)
# #' points(bg.pts, pch=21, bg=user.pts$bg.grp)

#' @name partitions
NULL

#' @rdname partitions
#' 
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
  bggrp1 <- bg[bg[, 2] <= horz & bg[, 1] < bvert,]
  bggrp2 <- bg[bg[, 2] <= horz & bg[, 1] >= bvert,]
  bggrp3 <- bg[bg[, 2] > horz & bg[, 1] <= tvert,]
  bggrp4 <- bg[bg[, 2] > horz & bg[, 1] > tvert,]
  
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

#' @rdname partitions
#' 
#' @export

get.checkerboard1 <- function(occs, envs, bg, aggregation.factor, quiet){
  if(is.null(envs)) stop("Cannot use checkerboard partitioning if envs is NULL.")
  occs <- as.data.frame(occs)
  rownames(occs) <- 1:nrow(occs)
  bg <- as.data.frame(bg)
  rownames(bg) <- 1:nrow(bg)
  
  grid <- raster::aggregate(envs[[1]], fact=aggregation.factor[1])
  w <- dismo::gridSample(occs, grid, n=1e6, chess='white')
  b <- dismo::gridSample(occs, grid, n=1e6, chess='black')
  bgw <- dismo::gridSample(bg, grid, n=1e6, chess='white')
  bgb <- dismo::gridSample(bg, grid, n=1e6, chess='black')
  
  if(nrow(w) > 0) { w$grp <- 1 }
  if(nrow(b) > 0) { b$grp <- 2 }
  r <- rbind(w, b)
  occ.grp <- r[order(as.numeric(rownames(r))),]$grp
  
  if(nrow(bgw) > 0) { bgw$grp <- 1 }
  if(nrow(bgb) > 0) { bgb$grp <- 2 }
  bgr <- rbind(bgw, bgb)
  bg.grp <- bgr[order(as.numeric(rownames(bgr))),]$grp
  
  # PATCH IF occs OR BG POINTS FALL INTO A SINGLE BIN
  noccgrp <- length(unique(occ.grp))
  nbggrp <- length(unique(bg.grp))
  if(noccgrp < 2 ){
    msg(paste("Warning: occurrence points fall in only", noccgrp, "bin"), quiet)
    bg.grp[ ! bg.grp %in% occ.grp] <- NA
    occ.grp <- as.numeric(as.factor(occ.grp))
    bg.grp <- as.numeric(as.factor(bg.grp))
  }
  
  if(length(unique(bg.grp[!is.na(bg.grp)])) != noccgrp) {
    stop("Error: occurrence records but no background points fall in 1 or more evaluation bin(s)")
  }
  
  out <- list(occ.grp=occ.grp, bg.grp=bg.grp)
  return(out)
}

#' @rdname partitions
#' 
#' @export

get.checkerboard2 <- function(occs, envs, bg, aggregation.factor, gridSampleN = 10000, quiet){
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
  occ.grp <- r[order(as.numeric(rownames(r))),]$grp
  
  bgr <- data.frame()
  if (nrow(bgww) > 0) bgww$grp <- 1; bgr <- rbind(bgr, bgww)
  if (nrow(bgwb) > 0) bgwb$grp <- 2; bgr <- rbind(bgr, bgwb)
  if (nrow(bgbw) > 0) bgbw$grp <- 3; bgr <- rbind(bgr, bgbw)
  if (nrow(bgbb) > 0) bgbb$grp <- 4; bgr <- rbind(bgr, bgbb)
  bg.grp <- bgr[order(as.numeric(rownames(bgr))),]$grp
  
  # PATCH IF occs OR BG POINTS FALL INTO FEWER THAN FOUR BINS
  noccgrp <- length(unique(occ.grp))
  nbggrp <- length(unique(bg.grp))
  if(noccgrp < 4 ){
    msg(paste("Warning: occurrence points fall in only", noccgrp, "bins"), quiet)
    bg.grp[ ! bg.grp %in% occ.grp] <- NA
    occ.grp <- as.numeric(as.factor(occ.grp))
    bg.grp <- as.numeric(as.factor(bg.grp))
  }
  
  if(length(unique(bg.grp[!is.na(bg.grp)])) != noccgrp) {
    stop("Error: occurrence records but no background points fall in 1 or more evaluation bin(s)")
  }
  
  out <- list(occ.grp=occ.grp, bg.grp=bg.grp)
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
  occ.grp <- 1:nrow(occs)
  bg.grp <- rep(0, nrow(bg))
  out <- list(occ.grp=occ.grp, bg.grp=bg.grp)
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
  occ.grp <- dismo::kfold(occs, kfolds)
  bg.grp <- rep(0, nrow(bg))
  out <- list(occ.grp=occ.grp, bg.grp=bg.grp)
  return(out)	
}