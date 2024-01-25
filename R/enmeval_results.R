#' Example ENMevaluation object.
#'
#' An example ENMevaluation object produced after running ENMevaluate() with 
#' feature classes L, LQ, LQH, H and regularization multipliers 1 through 5.
#' The SpatRaster data in the predictions slot is wrapped with `terra::wrap` to 
#' preserve the connections to the raster data. To use it, do the following:
#' `e@predictions <- terra::unwrap(e@predictions)`.
#'
#' @format An ENMevaluation object
"enmeval_results"