#' @title Methods to partition data for evaluation
#' 
#' @description \pkg{ENMeval} provides several ways to partition occurrence and background localities into bins for training and validation (or, evaluation and calibration). 
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
#' The \code{get.checkerboard} methods are variants of a checkerboard approach to partition occurrence localities. 
#' These methods use the \code{spatSample} function of the \pkg{terra} package (Hijmans 2023) to partition records according to checkerboard squares generated based on the input rasters. 
#' The spatial grain of these squares is determined by resampling (or aggregating) the original environmental input grids based on the user-defined \code{aggregation factor} (e.g., an aggregation factor with value 2 results in a checkerboard with grid cells four times the area of the original input rasters). 
#' With one input aggregation factor, \code{get.checkerboard} partitions data into two groups according to a 'basic' checkerboard pattern. With two aggregation factors, \code{get.checkerboard} partitions data into four groups according to 'hierarchical', or nested, checkerboard squares (see Muscarella et al. 2014). 
#' In contrast to the \code{get.block} method, the \code{get.checkerboard} methods subdivide geographic space equally but do not ensure a balanced number of occurrence localities in each group. 
#' The \code{get.checkerboard} methods give warnings (and potentially errors) if zero points (occurrence or background) fall in any of the expected bins.
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
#' @param envs SpatRaster: environmental predictor variables
#' @param orientation character vector: the order of spatial partitioning for the \code{get.block} method;
#' the first direction bisects the points into two groups, and the second direction bisects each of these further into two groups each, resulting in four groups; 
#' options are "lat_lon" (default), "lon_lat", "lon_lon", and "lat_lat"
#' @param aggregation.factor numeric or numeric vector: the scale of aggregation for \code{get.checkerboard}; can have one value (for 'basic') or two values (for 'hierarchical') -- see Details.
#' @param gridSampleN numeric: the number of points sampled from the input raster using gridSample() by the checkerboard partitioning functions
#' @param kfolds numeric: number of random \emph{k}-folds for \code{get.randomkfold} method
#' 
#' @return
#' A named list of two items:
#' \item{$occs.grp}{ A vector of bin designation for occurrence localities in the same order they were provided.}
#' \item{$bg.grp}{ A vector of bin designation for background localities in the same order they were provided.}
#' 
#' @note 
#' The \code{checkerboard} methods are designed to partition occurrence localities into spatial evaluation bins: two ('basic', for one aggregation factor) or four ('hierarchical', for two aggregation factors). 
#' They may give fewer bins, however, depending on where the occurrence localities fall with respect to the grid cells (e.g., all records happen to fall in one group of checkerboard squares). 
#' A warning is given if the number of bins is < 4 for the hierarchical method, and an error is given if all localities fall into a single evaluation bin.
#' 
#' @references
#' Hijmans, R. J. (2023). terra: Spatial Data Analysis. Available online at: \url{https://cran.r-project.org/package=terra}.
#' 
#' Muscarella, R., Galante, P. J., Soley‐Guardia, M., Boria, R. A., Kass, J. M., Uriarte, M., & Anderson, R. P. (2014). ENMeval: An R package for conducting spatially independent evaluations and estimating optimal model complexity for Maxent ecological niche models. \emph{Methods in Ecology and Evolution}, 5(11), 1198-1205. \doi{https://doi.org/10.1111/2041-210X.12945} 
#' 
#' Pearson, R. G., Raxworthy, C. J., Nakamura, M. and Peterson, A. T. (2007). Predicting species distributions from small numbers of occurrence records: a test case using cryptic geckos in Madagascar. \emph{Journal of Biogeography}, \bold{34}: 102-117. \doi{10.1111/j.1365-2699.2006.01594.x}
#' 
#' Radosavljevic, A., & Anderson, R. P. (2014). Making better Maxent models of species distributions: complexity, overfitting and evaluation. \emph{Journal of Biogeography}, \bold{41}: 629-643. \doi{10.1111/jbi.12227}
#' 
#' Shcheglovitova, M. and Anderson, R. P. (2013). Estimating optimal complexity for ecological niche models: a jackknife approach for species with small sample sizes. \emph{Ecological Modelling}, \bold{269}: 9-17. \doi{10.1016/j.ecolmodel.2013.08.011}
#' 
#' @author 
#' Robert Muscarella <bob.muscarella@gmail.com> and Jamie M. Kass <jamie.m.kass@gmail.com>
#' 
#' @examples 
#' library(terra)
#' 
#' set.seed(1)
#' 
#' ### Create environmental extent (raster)
#' envs <- rast(matrix(nrow=25, ncol=25))
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
#' ### Checkerboard partitioning method with aggregation factor of 4
#' chk.ag4 <- get.checkerboard(occs, envs, bg, aggregation.factor = 4)
#' plot(envs)
#' points(occs, pch=23, bg=chk.ag4$occs.grp)
#' plot(envs)
#' points(bg, pch=21, bg=chk.ag4$bg.grp)
#' # Higher aggregation factors result in bigger checkerboard blocks
#' chk.ag8 <- get.checkerboard(occs, envs, bg, aggregation.factor = 8)
#' plot(envs)
#' points(occs, pch=23, bg=chk.ag8$occs.grp)
#' plot(envs)
#' points(bg, pch=21, bg=chk.ag8$bg.grp)
#' 
#' ### Hierarchical checkerboard partitioning method with aggregation factors 
#' ### of 2 and 2
#' chk.ag2_2 <- get.checkerboard(occs, envs, bg, c(2,2))
#' plot(envs)
#' points(occs, pch=23, bg=chk.ag2_2$occs.grp)
#' plot(envs)
#' points(bg, pch=21, bg=chk.ag2_2$bg.grp)
#' # Higher aggregation factors result in bigger checkerboard blocks,
#' # and can vary between hierarchical levels
#' chk.ag4_6 <- get.checkerboard(occs, envs, bg, c(4,6))
#' plot(envs)
#' points(occs, pch=23, bg=chk.ag4_6$occs.grp)
#' plot(envs)
#' points(bg, pch=21, bg=chk.ag4_6$bg.grp)
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
    bggrp3 <- bg[bg[, 1] > vert & bg[, 1] < tvert,]
    bggrp4 <- bg[bg[, 1] >= tvert,]
  }else if(orientation == "lat_lat"){
    bhorz <- mean(c(max(grp1[, 2]), min(grp2[, 2])))
    thorz <- mean(c(max(grp3[, 2]), min(grp4[, 2])))
    horz <- mean(c(max(grpA[, 2]), min(grpB[, 2])))
    bggrp1 <- bg[bg[, 2] <= bhorz,]
    bggrp2 <- bg[bg[, 2] > bhorz & bg[, 2] <= horz,]
    bggrp3 <- bg[bg[, 2] > horz & bg[, 2] < thorz,]
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

get.checkerboard <- function(occs, envs, bg, aggregation.factor,
                              gridSampleN = 10000){
  if(is.null(envs)) stop("Cannot use checkerboard partitioning if envs is NULL.")
  if(length(aggregation.factor) == 1) {
    message("Generating basic checkerboard partitions...")
  } else if(length(aggregation.factor) == 2){
    message("Generating hierarchical checkerboard partitions...")
  } else {
    stop("You can only input one (for basic) or two (for hierarchical) aggregation factors.")
  }
  if(length(names(occs)) > 2) stop("Please input occs with only two columns: longitude and latitude.")
  if(length(names(bg)) > 2) stop("Please input bg with only two columns: longitude and latitude.")
  occs.v <- as.data.frame(occs) |> terra::vect(geom = colnames(occs))
  bg.v <- as.data.frame(bg) |> terra::vect(geom = colnames(bg))

  # make checkerboards  
  grid <- terra::aggregate(envs[[1]], fact=aggregation.factor[1])
  occs.w <- terra::spatSample(occs.v, size = gridSampleN, strata = grid, 
                              chess='white')
  occs.b <- terra::spatSample(occs.v, size = gridSampleN, strata = grid, 
                              chess='black')
  bg.w <- terra::spatSample(bg.v, size = gridSampleN, strata = grid, 
                            chess='white')
  bg.b <- terra::spatSample(bg.v, size = gridSampleN, strata = grid, 
                            chess='black')
  
  # for hierarchical checkerboards, which need two aggregation factors
  if(length(aggregation.factor) == 2) {
    grid2 <- terra::aggregate(grid, fact=aggregation.factor[2])
    occs.ww <- terra::spatSample(occs.w, size = gridSampleN, strata = grid2, 
                                 chess='white')
    occs.wb <- terra::spatSample(occs.w, size = gridSampleN, strata = grid2, 
                                 chess='black')
    occs.bb <- terra::spatSample(occs.b, size = gridSampleN, strata = grid2, 
                                 chess='black')
    occs.bw <- terra::spatSample(occs.b, size = gridSampleN, strata = grid2, 
                                 chess='white')
    bg.ww <- terra::spatSample(bg.w, size = gridSampleN, strata = grid2, 
                               chess='white')
    bg.wb <- terra::spatSample(bg.w, size = gridSampleN, strata = grid2, 
                               chess='black')
    bg.bb <- terra::spatSample(bg.b, size = gridSampleN, strata = grid2, 
                               chess='black')
    bg.bw <- terra::spatSample(bg.b, size = gridSampleN, strata = grid2, 
                               chess='white')
  }
  
  # assemble groups
  if (length(aggregation.factor) == 1) {
    occs.w <- terra::crds(occs.w, df = TRUE)
    occs.b <- terra::crds(occs.b, df = TRUE)
    if(nrow(occs.w) > 0) { occs.w$grp <- 1 }
    if(nrow(occs.b) > 0) { occs.b$grp <- 2 }
    occs.r <- rbind(occs.w, occs.b)
    occs.grp <- occs.r[order(as.numeric(rownames(occs.r))),]$grp
    
    bg.w <- terra::crds(bg.w, df = TRUE)
    bg.b <- terra::crds(bg.b, df = TRUE)
    if(nrow(bg.w) > 0) { bg.w$grp <- 1 }
    if(nrow(bg.b) > 0) { bg.b$grp <- 2 }
    bg.r <- rbind(bg.w, bg.b)
    bg.grp <- bg.r[order(as.numeric(rownames(bg.r))),]$grp
  } else if (length(aggregation.factor) == 2) {
    occs.ww <- terra::crds(occs.ww, df = TRUE)
    occs.bw <- terra::crds(occs.bw, df = TRUE)
    occs.wb <- terra::crds(occs.wb, df = TRUE)
    occs.bb <- terra::crds(occs.bb, df = TRUE)
    occs.r <- data.frame()
    if (nrow(occs.ww) > 0) occs.ww$grp <- 1; occs.r <- rbind(occs.r, occs.ww)
    if (nrow(occs.wb) > 0) occs.wb$grp <- 2; occs.r <- rbind(occs.r, occs.wb)
    if (nrow(occs.bw) > 0) occs.bw$grp <- 3; occs.r <- rbind(occs.r, occs.bw)
    if (nrow(occs.bb) > 0) occs.bb$grp <- 4; occs.r <- rbind(occs.r, occs.bb)
    occs.grp <- occs.r[order(as.numeric(rownames(occs.r))),]$grp
    
    bg.ww <- terra::crds(bg.ww, df = TRUE)
    bg.bw <- terra::crds(bg.bw, df = TRUE)
    bg.wb <- terra::crds(bg.wb, df = TRUE)
    bg.bb <- terra::crds(bg.bb, df = TRUE)
    bg.r <- data.frame()
    if (nrow(bg.ww) > 0) bg.ww$grp <- 1; bg.r <- rbind(bg.r, bg.ww)
    if (nrow(bg.wb) > 0) bg.wb$grp <- 2; bg.r <- rbind(bg.r, bg.wb)
    if (nrow(bg.bw) > 0) bg.bw$grp <- 3; bg.r <- rbind(bg.r, bg.bw)
    if (nrow(bg.bb) > 0) bg.bb$grp <- 4; bg.r <- rbind(bg.r, bg.bb)
    bg.grp <- bg.r[order(as.numeric(rownames(bg.r))),]$grp
  }
  
  # PATCH IF occs OR BG POINTS FALL INTO A SINGLE BIN
  noccs.grp <- length(unique(occs.grp))
  nbg.grp <- length(unique(bg.grp))
  if(noccs.grp < 2 ){
    message(paste("Warning: occurrence points fall in only", noccs.grp, "bin"))
    bg.grp[ ! bg.grp %in% occs.grp] <- NA
    occs.grp <- as.numeric(as.factor(occs.grp))
    bg.grp <- as.numeric(as.factor(bg.grp))
  }
  
  if(length(unique(bg.grp[!is.na(bg.grp)])) != noccs.grp) {
    stop("Error: occurrence records but no background points fall in 1 or more evaluation bin(s)")
  }
  
  out <- list(occs.grp=occs.grp, bg.grp=bg.grp)
  return(out)
}

#' @export
get.checkerboard1 <- function(occs, envs, bg, aggregation.factor, gridSampleN = 10000) {
  message("This function is deprecated and will be phased out with the next ENMeval version. Please use get.checkerbaord instead.")
  return(get.checkerboard(occs, envs, bg, aggregation.factor, gridSampleN = 10000))
}

#' @export
get.checkerboard2 <- function(occs, envs, bg, aggregation.factor, gridSampleN = 10000) {
  message("This function is deprecated and will be phased out with the next ENMeval version. Please use get.checkerbaord instead.")
  return(get.checkerboard(occs, envs, bg, aggregation.factor, gridSampleN = 10000))
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
  occs.grp <- predicts::folds(occs, kfolds)
  bg.grp <- rep(0, nrow(bg))
  out <- list(occs.grp=occs.grp, bg.grp=bg.grp)
  return(out)	
}

