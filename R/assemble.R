assemble.stats.parts <- function(nk, train.stats.all, val.stats.all, tune.names, tune.tbl) {
  # if jackknife cross-validation (leave-one-out), correct variance for
  # non-independent samples (Shcheglovitova & Anderson 2013)
  if(partitions == "jackknife") {
    sum.list <- list(avg = mean, sd = ~sqrt(corrected.var(., nk)))
  }else{
    # if other partition method, leave sd as is
    sum.list <- list(avg = mean, sd = sd)
  } 
  
  # if there is one partition, or if using an testing dataset, do not take summary statistics
  if(nk == 1 | partitions == "testing") sum.list <- list(function(x) {x})
  
  # if tune.tbl exists, make tune.args column a factor to keep order after using dplyr functions
  if(!is.null(tune.tbl)) val.stats.all$tune.args <- factor(val.stats.all$tune.args, levels = tune.names)
  
  # calculate summary statistics
  cv.stats.sum <- val.stats.all |> 
    dplyr::group_by(tune.args) |>
    dplyr::select(-fold) |> 
    dplyr::summarize_all(sum.list) |>
    dplyr::ungroup() 
  
  # change names (replace _ with .)
  names(cv.stats.sum) <- gsub("(.*)_(.*)", "\\1.\\2", names(cv.stats.sum))
  # order columns alphabetically
  cv.stats.sum <- cv.stats.sum[, order(colnames(cv.stats.sum))]
  
  # if tune.tbl exists
  if(!is.null(tune.tbl)) {
    # make tune.args column in training stats factor too for smooth join
    train.stats.all$tune.args <- factor(train.stats.all$tune.args, levels = tune.names)
    # eval.stats is the join of tune.tbl, training stats, and cv stats
    eval.stats <- tune.tbl |> dplyr::left_join(train.stats.all, by = "tune.args") |>
      dplyr::left_join(cv.stats.sum, by = "tune.args")
  }else{
    # if tune.tbl does not exist, eval.stats is the binding of training stats to cv stats
    train.stats.all$tune.args <- NULL
    cv.stats.sum$tune.args <- NULL
    eval.stats <- dplyr::bind_cols(train.stats.all, cv.stats.sum)
  }
  return(eval.stats)
}

assemble.stats.noParts <- function(train.stats.all, val.stats.all, tune.names, tune.tbl) {
  # make tune.args column in training stats factor too for smooth join
  train.stats.all$tune.args <- factor(train.stats.all$tune.args, levels = tune.names)
  # if no partitions assigned, eval.stats is the join of tune.tbl to training stats
  eval.stats <- dplyr::left_join(tune.tbl, train.stats.all, by = "tune.args") 
  if(nrow(val.stats.all) > 0) eval.stats <- dplyr::left_join(eval.stats, val.stats.all, by = "tune.args")
  if("fold" %in% names(eval.stats)) eval.stats <- eval.stats |> dplyr::select(-fold)
  return(eval.stats)
}

assemble.overlap <- function(mod.full.pred.all, overlapStat) {
  nr <- terra::nlyr(mod.full.pred.all)
  if(nr == 0) {
    message("Warning: calculate range overlap without model prediction rasters.")
  }else if(nr == 1) {
    message("Warning: only 1 model prediction raster found. Need at least 2 rasters to calculate range overlap. Increase number of tuning arguments and run again.") 
  }else{
    ls <- list()
    for(ovStat in overlapStat) {
      message(paste0("Calculating range overlap for statistic ", ovStat, "..."))
      # turn negative values to 0 for range overlap calculations
      predictions.noNegs <- terra::rast(lapply(mod.full.pred.all, function(x) {x[x<0] <- 0; x}))
      overlap.mat <- calc.niche.overlap(predictions.noNegs, ovStat)
      ls[[ovStat]] <- overlap.mat
    }
  }
  return(ls)
}
