#' @title Compare model accuracy metrics of Ecological Niche Models (ENMs) built with different set of predictors.
#' #' @description \code{ENMnulls.test()} Iteratively builds null ENMs for "k" sets of user-specified model
#' settings bsaed on "k" input ENMevaluation objects, from which all other analysis 
#' settings are extracted.Summary statistics of the performance metrics for the null ENMs 
#' are taken (averages and standard deviations) and effect sizes and p-values are calculated by 
#' comparing these summary statistics to the empirical values of the performance metrics 
#' (i.e., from the model built with the empirical data). This is an extension of {ENMnulls()} for comparisons
#' of two or more predictor variable sets.
#' 
#' @param e.list: list of ENMevaluation objects to be compared
#' @param mod.settings.list named list: model settings corresponding to ENMevaluation objects in e.list 
#' that specify the settings used for building null models
#' @param no.iter numeric: number of null model iterations.
#' @param eval.stats character: model accuarcy metrics that will be used to calculate null model statistics. 
#' Can be “auc.val”, “auc.diff”, “cbi.val”, “or.mtp”, “or.10p”.
#' @param user.enm ENMdetails object: if implementing a user-specified model. 
#' @user.eval.type character: if implementing a user-specified model -- either "knonspatial", "kspatial", 
#' "testing" or "none".
#' @param userStats.signs named list: user-defined evaluation statistics attributed with either 1 or -1 
#' to designate whether the expected difference between empirical and null models is positive or negative; 
#' this is used to calculate the p-value of the z-score when comparing two predictor variable sets. Default is NULL.
#' @param removeMxTemp boolean: if TRUE, delete all temporary data generated when using maxent.jar for modeling
#' @param parallel boolean: if TRUE, use parallel processing.
#' @param numCores numeric: number of cores to use for parallel processing; if NULL, all available cores will be used.
#' @param parallelType character: either "doParallel" or "doSNOW" (default: "doSNOW").
#' @param quiet boolean: if TRUE, silence all function messages (but not errors).
#' 
#' #' @details This null ENM technique extends the implementation in Bohl \emph{et al.} (2019)and Kass \emph{et al.} (2020),
#' which follows the original methodology of Raes & ter Steege (2007). Here we evaluate if observed differences in accuracy metric values 
#' (e.g., AUC, omission rates, CBI) of empirical models built with different sets of predictor variable are greater than expected 
#' at random. This is done by building the null distributions of the difference in accuracy metrics
#  employing the same withheld validation data used to evaluate the empirical models. Please see the vignette for a brief example: <
#' 
#' This function avoids using raster data to speed up each iteration, and instead samples null occurrences from the 
#' partitioned background records. Thus, you should avoid running this when your background records are not well 
#' sampled across the study extent, as this limits the extent that null occurrences can be sampled from.
#' 
#' @references 
#' Bohl, C. L., Kass, J. M., & Anderson, R. P. (2019). A new null model approach to quantify performance and significance for ecological niche models of species distributions. \emph{Journal of Biogeography}, \bold{46}: 1101-1111. \url{https://doi.org/10.1111/jbi.13573}
#' 
#' Kass, J. M., Anderson, R. P., Espinosa-Lucas, A., Juárez-Jaimes, V., Martínez-Salas, E., Botello, F.,  Tavera, G., Flores-Martínez, J. J., & Sánchez-Cordero, V. (2020). Biotic predictors with phenological information improve range estimates for migrating monarch butterflies in Mexico. \emph{Ecography}, \bold{43}: 341-352. \url{https://doi.org/10.1111/ecog.04886}
#'
#' Raes, N., & ter Steege, H. (2007). A null-model for significance testing of presence-only species distribution models. \emph{Ecography}, \bold{30}: 727-736. \url{https://doi.org/10.1111/j.2007.0906-7590.05041.x} 
#' 
#' #' @return An \code{ENMnull}An ENMnull object with slots containing evaluation summary statistics 
#' for the null models and their cross-validation results, as well as differences in results between the 
#' empirical and null models. This comparison table includes T-statistics for pairwise comparisons (T-test) 
#' and F-statistic (ANOVA) of these differences and their associated p-values (under a normal distribution). 

ENMnulls.test <- 	function(e.list, mod.settings.list, no.iter, 
                           eval.stats = c("auc.val", "auc.diff","cbi.val","or.mtp","or.10p"),
                           user.enm = NULL, user.eval.type = NULL, user.Stats.signs = NULL,
                           removeMxTemp = TRUE, parallel = FALSE, numCores = NULL, 
                           parallelType = "doSnow", quiet = FALSE){
  
  ## checks
  #more than one input enm
  if(length(e.list) == 1){
    stop("Please input more than one ENM to run comparisons.")
  }
  
  # model settings are all single entries for each enm treatment
  for(i in 1:length(mod.settings.list)){
    if(!all(sapply(mod.settings.list[[i]], length) == 1)){
      stop("Please input a single set of model settings.")
    }
  }
  # model settings are correct for input algorithm and are entered in the right order -- 
  #if not, put them in the right order, else indexing models later will fail because the model 
  # name will be incorrect
  for(i in 1:length(e.list)){
    if(e.list[[i]]@algorithm %in% c("maxent.jar", "maxnet")) {
      if(length(mod.settings[[i]]) != 2){
      stop("Please input two complexity settings (fc [feature classes] and rm [regularization
           multipliers]) for mod.settings for maxent.jar and maxnet models.")
        }
    if(all(names(mod.settings[[i]]) %in% c("fc", "rm"))) {
      if(!all(names(mod.settings[[i]]) == c("fc", "rm"))) {
        mod.settings <- mod.settings[c("fc", "rm")]
      }
    }else{
      stop('Please input only "fc" (feature classes) and "rm" (regularization multipliers) for
           mod.settings for maxent.jar and maxnet models.')
    }
  #}else if(e@algorithm == "bioclim") {
  #  if(length(mod.settings) != 1) {
  #    stop("Please input one complexity setting (tails) for mod.settings for BIOCLIM models.")
  #  }
  #  if(!all(names(mod.settings) == "tails")) {
  #    stop('Please input only "tails" for mod.settings for BIOCLIM models.')
  #  }
    }
  }
  
    # assign evaluation type based on partition method
    if(is.null(user.eval.type)) {
      eval.type <- switch(e.list[[1]]@partition.method,
                          randomkfold = "knonspatial",
                          jackknife = "knonspatial",
                          block = "kspatial",
                          checkerboard1 = "kspatial",
                          checkerboard2 = "kspatial",
                          testing = "testing",
                          none = "none")  
    }else{
      eval.type <- user.eval.type
    }
  #Create list of background SWD dataframes for each ENM treatment
  bg.list <- lapply(e.list, function(e){
    #Get bg points env values from ENMevaluation object
    bg<- e@bg
    #Get bg partitions groups
    bg.grp <- e@bg.grp
    #Concatenate bg points, env values & partition groups
    bg <- cbind(bg, bg.grp)
    return(bg)
  })
  
  #Create list of "i" sets of null "n" occurrences for each ENMevaluation object in e.list
  
  #Iteratively apply ENMulls
  #Extract relevant model accuracy metrics
  #Run statistical tests
  #Run post-hoc test (one-tailed z-test) for pairwise differences
}
