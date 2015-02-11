ENMevaluate <- function (occ, env, bg.coords = NULL, occ.grp = NULL, bg.grp = NULL, 
                         RMvalues = seq(0.5, 4, 0.5), fc = c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"), 
                         categoricals = NULL, n.bg = 10000, method = NULL, overlap = FALSE, 
                         aggregation.factor = c(2, 2), kfolds = NA, bin.output = FALSE, 
                         rasterPreds = TRUE, updateProgress = NULL) 
{
  if (is.null(method)) {
    stop("Evaluation method needs to be specified.")
  }
  if (is.null(bg.coords)) {
    bg.coords <- randomPoints(env[[1]], n = n.bg)
  }
  maxent.args <- make.args(RMvalues, fc)
  args.lab <- make.args(RMvalues, fc, labels = TRUE)
  occ <- as.data.frame(occ)
  colnames(occ) <- c("LON", "LAT")
  bg.coords <- as.data.frame(bg.coords)
  colnames(bg.coords) <- c("LON", "LAT")
  message(ifelse(method == "jackknife", "Doing evaluations using k-1 jackknife...", 
                 ifelse(method == "checkerboard1", "Doing evaluations using checkerboard 1...", 
                        ifelse(method == "checkerboard2", "Doing evaluations using checkerboard 2...", 
                               ifelse(method == "block", "Doing evaluations using spatial blocks...", 
                                      ifelse(method == "randomkfold", "Doing random k-fold evaluation groups...", 
                                             ifelse(method == "user", "Doing user-defined evaluation groups...", 
                                                    "Error: You need to specify an accepted evaluation method. Check the documentation.")))))))
  results <- tuning(occ, env, bg.coords, occ.grp, bg.grp, method, 
                    maxent.args, args.lab, categoricals, aggregation.factor, 
                    kfolds, bin.output, rasterPreds, updateProgress)
  if (overlap == TRUE) {
    if (length(maxent.args) > 1) {
      message("Calculating niche overlap")
      overlap.mat <- calc.niche.overlap(results@predictions, 
                                        "D")
      results@overlap <- overlap.mat
    }
    else {
      message("Need >1 settings to do niche overlap")
    }
  }
  return(results)
}