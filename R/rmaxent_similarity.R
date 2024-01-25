#' @title Calculate Multivariate Environmental Similarity
#'
#' @description NOTICE: This function was borrowed from the rmaxent package 
#' written by John Baumgartner (https://github.com/johnbaums/rmaxent/).
#  It is included here with John's permission to make ENMeval CRAN-compatible 
#' (dependencies on Github-only packages are not allowed for CRAN).
#' 
#' Calculate Multivariate Environmental Similarity and most dissimilar/similar 
#' variables with respect to a reference dataset, for a set of environmental 
#' variables.
#'
#' @param x a `SpatRaster`, `list`, `matrix`, or `data.frame`
#'   where each layer/column/element represents focal values of an environmental
#'   variable.
#' @param ref a `list`, `matrix`, or `data.frame` where each
#'   column/element represents reference values for an environmental variable
#'   (corresponding to those given in `x`).
#' @param full (logical) should similarity values be returned for all variables?
#'   If `FALSE` (the default), then only the minimum similarity scores
#'   across variables will be returned.
#' @return If `x` is a `SpatRaster` object, this function returns a list 
#'   containing:
#'   - `similarity`: a `SpatRaster` giving the environmental similarities for
#'   each variable in `x` (only included when `full=TRUE`); 
#'   - `similarity_min`: a `SpatRaster` giving the minimum similarity value 
#'   across all variables for each location (i.e. the MESS);
#'   - `mod`: a factor `SpatRaster` indicating which variable was most 
#'   dissimilar to its reference range (i.e. the MoD map, 
#'   Elith \emph{et al.} 2010); and
#'   - `mos`: a factor `SpatRaster` indicating which variable was most 
#'   similar to its reference range.
#'   
#'   If `x` is a `list`, `matrix`, or `data.frame`, the function will return
#'   a list as above, but with `SpatRaster` objects replaced by matrix and 
#'   vectors. If a predictor variable is missing from the "most similar" or
#'   "most different" results, it will not be included in the output, and 
#'   thus will not plot when ENMeval plotting functions are used.
#' @details `similarity` uses the MESS algorithm described in Appendix S3
#'   of Elith \emph{et al.} 2010.
#' @references 
#' Elith, J., Kearney, M., and Phillips, S. (2010) The art of modelling range-shifting species. \emph{Methods in Ecology and Evolution}, \bold{1}: 330-342. \doi{doi:10.1111/j.2041-210X.2010.00036.x}
#' 
#' @importFrom methods is
#' @importFrom stats na.omit
#' @examples
#' library(predicts)
#' library(terra)
#' occs <- read.csv(file.path(system.file(package="predicts"), 
#' "/ex/bradypus.csv"))[,2:3]
#' predictors <- rast(list.files(path = 
#' paste(system.file(package='predicts'), '/ex', sep=''), 
#' ref <- extract(predictors, occs, ID = FALSE)
#' mess <- similarity(predictors, ref, full = TRUE)
#' 
#' \dontrun{
#' library(rasterVis)
#' library(RColorBrewer)
#' levelplot(mess$mod, col.regions=brewer.pal(8, 'Set1'))
#' levelplot(mess$mos, col.regions=brewer.pal(8, 'Set1'))
#' }
#' @export

similarity <- function(x, ref, full=FALSE) {
  if(!methods::is(ref, 'data.frame')) {
    ref <- as.data.frame(ref)
  }
  if(is(x, 'SpatRaster')) {
    r <- TRUE
    if(isTRUE(full)) {
      out <- terra::rast(replicate(
        terra::nlyr(x), terra::init(x, NA)))
    } else {
      out <- terra::init(x, fun = NA)[[1]]
    }
  } else r <- FALSE
  ref <- stats::na.omit(ref)
  if(!methods::is(x, 'data.frame')) {
    x <- terra::as.data.frame(x, na.rm = FALSE)
  }
  if(is.null(dim(ref))) {
    rng <- as.data.frame(range(ref, na.rm=TRUE))
  } else {
    rng <- as.data.frame(apply(ref, 2, range, na.rm=TRUE))
  }
  pct_less <- mapply(function(x, ref) {
    findInterval(x, sort(ref))/length(ref)
  }, x, ref, SIMPLIFY=FALSE)
  sim <- mapply(function(f, rng, p) {
    ifelse(f==0, (p-rng[1])/diff(rng)*100,
           ifelse(f > 0 & f <= 0.5, f*200,
                  ifelse(f > 0.5 & f < 1, (1-f)*200,
                         (rng[2]-p)/diff(rng)*100)))
  }, pct_less, rng, x)
  
  min_sim <- if(is.matrix(sim)) apply(sim, 1, min) else(min(sim))
  
  mins <- apply(sim, 1, which.min)
  most_dissimilar_vec <- unlist(ifelse(lengths(mins)==0, NA, mins))
  maxs <- apply(sim, 1, which.max)
  most_similar_vec <- unlist(ifelse(lengths(maxs)==0, NA, maxs))
  
  if(isTRUE(r)) {
    most_dissimilar <- terra::rast(out)
    most_dissimilar[] <- most_dissimilar_vec
    most_dissimilar <- terra::as.factor(most_dissimilar)
    levels(most_dissimilar)[[1]] <- data.frame(ID=seq_len(ncol(sim)), 
                                               var=colnames(sim))
    most_similar <- terra::rast(out)
    most_similar[] <- most_similar_vec
    most_similar <- terra::as.factor(most_similar)
    levels(most_similar)[[1]] <- data.frame(ID=seq_len(ncol(sim)), 
                                            var=colnames(sim))  
    
    out_min <- terra::rast(out)
    out_min[] <- min_sim
    names(out_min) <- "minimum_similarity"
    if(isTRUE(full)) {
      out[] <- sim
      list(similarity=out, similarity_min=out_min, mod=most_dissimilar, 
           mos=most_similar)
    } else list(similarity_min=out_min, mod=most_dissimilar, mos=most_similar)
  } else {
    if(isTRUE(full)) {
      list(similarity=sim, similarity_min=min_sim, 
           mod=most_dissimilar_vec, mos=most_similar_vec)
    } else list(similarity_min=min_sim, mod=most_dissimilar_vec, 
                mos=most_similar_vec)
  }
}
