assign.partitions <- function(occs, bg, envs, partitions, partition.settings) {
  switch(partitions, 
         jackknife = get.jackknife(occs, bg),
         randomkfold = get.randomkfold(occs, bg, partition.settings$kfolds),
         block = get.block(occs, bg, partition.settings$orientation),
         checkerboard = get.checkerboard(occs, envs, bg, partition.settings$aggregation.factor),
         user = user.grp,
         testing = list(occs.grp = rep(0, nrow(occs)), bg.grp = rep(0, nrow(bg))),
         none = list(occs.grp = rep(0, nrow(occs)), bg.grp = rep(0, nrow(bg))))
}

assign.partition.msg <- function(partitions, partition.settings) {
  switch(partitions,
         jackknife = "* Model evaluations with k-1 jackknife (leave-one-out) cross validation...",
         randomkfold = paste0("* Model evaluations with random ", 
                              partition.settings$kfolds, 
                              "-fold cross validation..."),
         block =  paste0("* Model evaluations with spatial block (4-fold) cross validation and ", 
                         partition.settings$orientation, " orientation..."),
         checkerboard = ifelse(length(partition.settings$aggregation.factor) == 1, 
                               "* Model evaluations with basic checkerboard (2-fold) cross validation...",
                               "* Model evaluations with hierarchical checkerboard (4-fold) cross validation..."),
         user = paste0("* Model evaluations with user-defined ", 
                       length(unique(user.grp$occs.grp)), 
                       "-fold cross validation..."),
         testing = "* Model evaluations with testing data...",
         none = "* Skipping model evaluations (only calculating full model statistics)...")
}

  

  
  