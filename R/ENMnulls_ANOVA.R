#' @title Compare model accuracy metrics of Ecological Niche Models (ENMs) built with different set of predictors.
#' @description \code{ENMnulls_ANOVA()} Iteratively builds null ENMs for "k" sets of user-specified model
#' settings based on "k" input ENMevaluation objects, from which all other analysis 
#' settings are extracted.Summary statistics of the performance metrics for the null ENMs 
#' are taken (averages and standard deviations) and effect sizes and p-values are calculated by 
#' comparing these summary statistics to the empirical values of the performance metrics 
#' (i.e., from the model built with the empirical data). This is an extension of {ENMnulls()} for comparisons
#' of two or more predictor variable sets.
#' 
#' @param e.list list of ENMevaluation objects to be compared
#' @param mod.settings.list named list: model settings corresponding to ENMevaluation objects in e.list 
#' that specify the settings used for building null models
#' @param no.iter numeric: number of null model iterations.
#' @param eval.stats character: model accuarcy metrics that will be used to calculate null model statistics. 
#' Can be one of “auc”, “cbi”, “or.mtp”, “or.10p”.
#' @param user.eval.type character: if implementing a user-specified model -- either "knonspatial", "kspatial", 
#' "testing" or "none".
#' @param alternative a character string indicating the type of test for post-hoc analyses. Can be one of "two.sided" (default),
#' "greater", "less".  
#' @param userStats.signs named list: user-defined evaluation statistics attributed with either 1 or -1 
#' to designate whether the expected difference between empirical and null models is positive or negative; 
#' this is used to calculate the p-value of the z-score when comparing two predictor variable sets. Default is NULL.
#' 
#' @details This null ENM technique extends the implementation in Bohl \emph{et al.} (2019)and Kass \emph{et al.} (2020),
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
#' @return An \code{ENMnull}An ENMnull object with slots containing evaluation summary statistics 
#' for the null models and their cross-validation results, as well as differences in results between the 
#' empirical and null models. This comparison table includes T-statistics for pairwise comparisons (T-test) 
#' and F-statistic (ANOVA) of these differences and their associated p-values (under a normal distribution). 
#' @export

ENMnulls_ANOVA <- function(e.list, mod.settings.list, 
                           eval.stats = c("auc.val","auc.diff","cbi.val","or.mtp","or.10p"),
                           alternative = "two.sided",
                           no.iter, user.eval.type = NULL, userStats.signs = NULL) {
  
  # z equal number of treatments
  z <- length(e.list)
  
  # assign directionality of sign for evaluation stats
  signs <- c(list("auc.val" = 1, "auc.train" = 1, "cbi.val" = 1, "cbi.train" = 1,
                  "auc.diff" = -1, "or.10p" = -1, "or.mtp" = -1), userStats.signs)
  
  # check that e.list occs and bg are the same
  e.occs <- list()
  e.bg <- list()
  for(i in 1:z) {
    e.occs.i <- e.list[[i]]@occs[,c("longitude", "latitude")] %>% dplyr::arrange(longitude, latitude)
    e.occs[[i]] <- unname(apply(e.occs.i, 1, paste, collapse = "|"))
    e.bg.i <- e.list[[i]]@bg[,c("longitude", "latitude")] %>% dplyr::arrange(longitude, latitude)
    e.bg[[i]] <- unname(apply(e.bg.i, 1, paste, collapse = "|"))
  }
  if(!all(length(e.occs[[1]]) == sapply(e.occs, length))) {
    stop("At least one occurrence dataset has a different number of rows from the others.
         All ENMevaluation objects must have the same number of occurrence records.")
  }
  if(!all(length(e.bg[[1]]) == sapply(e.bg, length))) {
    stop("At least one background dataset has a different number of rows from the others.
         All ENMevaluation objects must have the same number of background records.")
  }
  e.occs.df <- dplyr::bind_cols(e.occs)
  e.bg.df <- dplyr::bind_cols(e.bg)
  e.occs.allEq <- sum(apply(e.occs.df, 1, function(x) x[1] %in% x)) == nrow(e.occs.df)
  e.bg.allEq <- sum(apply(e.bg.df, 1, function(x) x[1] %in% x)) == nrow(e.bg.df)
  if(e.occs.allEq != TRUE) {
    stop("At least one occurrence record has coordinates that differ between datasets.
         Occurrence record coordinates must be identical among datasets.")
  }
  if(e.bg.allEq != TRUE) {
    stop("At least one background record has coordinates that differ between datasets.
         Background record coordinates must be identical among datasets.")
  }
  
  # Check that partition.method is the same or just includes "user".
  e.part.methods <- unique(sapply(e.list, function(x) x@partition.method))
  if(length(e.part.methods) == 2 & "user" %in% e.part.methods) {
    message('Partition methods include "user": make sure partitions are identical across ENMevaluation objects.')
  }else if(length(e.part.methods) != 1) {
    stop("When statistically comparing ENMevaluation objects, they cannot use different partition methods.")
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
  
  
  # make a list to hold the same null occs for each treatment
  null.occs.treatments <- list()
  # The first list will be for the first ENMevaluation object. This one will 
  # have the predictor values of e1 in the table.
  null.occs.treatments[[1]] <- list()
  for(i in 1:no.iter) {
    null.occs.treatments[[1]][[i]] <- make_null_occs(e.list[[1]], eval.type)  
  }
  # To avoid duplication of work, we keep these values for e1 and associate the 
  # predictor values of e2, e3, ... based on the same coordinates.
  for(x in 2:z) {
    e.x <- e.list[[x]]
    null.occs.treatments[[x]] <- list()
    for(i in 1:no.iter) {
      null.occs.xy.e1.i <- null.occs.treatments[[1]][[i]][,c("longitude", "latitude")]
      occs.and.bg.x <- cbind(rbind(e.x@occs, e.x@bg), grp = c(e.x@occs.grp, e.x@bg.grp))
      null.occs.treatments[[x]][[i]] <- merge(occs.and.bg.x, null.occs.xy.e1.i, by = c("longitude", "latitude"))
    }
  }
  
  # Now, make ENMnull objects for each treatment, which have the same null occs.
  # user.eval.type must be set here because partition.method may vary across
  # treatments -- for example, one may run ENMevaluate on one predictor treatment
  # with some partitions settings, then want to input those same settings for
  # other treatment and thus specify partition.method as "user". However, the
  # user.eval.type is the same.
  ENMnull.list <- list()
  for(i in 1:z) {
    message(paste0("* Predictor treatment ", i, ":"))
    ENMnull.list[[i]] <- ENMnulls(e.list[[i]], mod.settings.list[[i]], no.iter, 
                                  input.random.data = null.occs.treatments[[i]],
                                  user.enm = NULL, user.eval.type = user.eval.type, 
                                  userStats.signs = NULL, removeMxTemp = TRUE, 
                                  parallel = FALSE, numCores = NULL, 
                                  parallelType = "doSNOW", quiet = FALSE)
  }
  
  # get empirical model evaluation statistics for comparison
  emp.diff <- lapply(1:z, function(x) ENMnull.list[[x]]@null.emp.results[1,] %>%
                       dplyr::select(dplyr::contains(eval.stats)) %>%
                       dplyr::rename_with(~paste0(., ".avg")))
  
  names(emp.diff) <- LETTERS[1:z]
  
  null.results.all <- lapply(1:z, function(x) ENMnull.list[[x]]@null.results %>%
                               dplyr::select(dplyr::contains(eval.stats) & dplyr::ends_with("avg")))
  names(null.results.all) <- LETTERS[1:z]
  
  # Make code names for different combinations of treatments and assign them.
  # Then calculate the differences in evaluation stats.
  comb <- combn(LETTERS[1:z], 2)
  null.results.diff.list <- list()
  emp.diff.list <- list()
  for(i in 1:ncol(comb)) {
    comb.i <- comb[,i]
    if(alternative == "two.sided"){
      emp.diff.list[[i]] <- abs(emp.diff[[comb.i[1]]] - emp.diff[[comb.i[2]]])
      null.results.diff.list[[i]] <- abs(null.results.all[[comb.i[1]]] - null.results.all[[comb.i[2]]])
      # nulls.diff <- null.results.all %>% dplyr::transmute(`null.B-A` = abs(B - A)) %>% dplyr::mutate(iter = 1:nrow(.))  
    }else{
      emp.diff.list[[i]] <- emp.diff[[comb.i[1]]] - emp.diff[[comb.i[2]]]
      null.results.diff.list[[i]] <- null.results.all[[comb.i[1]]] - null.results.all[[comb.i[2]]]
      # nulls.diff <- nulls %>% dplyr::transmute(`null.B-A` = B - A) %>% dplyr::mutate(iter = 1:nrow(.))  
    }
    null.results.diff.list[[i]] <- null.results.diff.list[[i]] %>% dplyr::mutate(iter = 1:nrow(.),
                                                                                 comb = paste(comb.i, collapse = "_"))
    emp.diff.list[[i]] <- dplyr::mutate(emp.diff.list[[i]], comb = paste(comb.i, collapse = "_"))
  }
  
  # Combine into one table
  emp.diff.comb <- dplyr::bind_rows(emp.diff.list)%>% 
    dplyr::rename(treatment.dif = comb)%>%
    dplyr::mutate(treatment.dif = as.factor(treatment.dif))
  
  null.results.diff.comb <- dplyr::bind_rows(null.results.diff.list) %>% 
    dplyr::rename(treatment.dif = comb)
  
  # Get averages and sds
  nulls.diff.avg <- null.results.diff.comb %>% dplyr::select(-iter) %>% 
    dplyr::group_by(treatment.dif) %>% dplyr::summarise_all(mean, na.rm = TRUE)
  
  nulls.diff.avg.pool <- null.results.diff.comb %>% dplyr::select(-iter, -treatment.dif) %>% 
    dplyr::summarise_all(mean, na.rm = TRUE)
  
  nulls.diff.sd <- null.results.diff.comb %>% dplyr::select(-iter) %>% 
    dplyr::group_by(treatment.dif) %>% dplyr::summarise_all(sd, na.rm = TRUE)
  
  nulls.diff.sd.pool <- null.results.diff.comb %>% dplyr::select(-iter, -treatment.dif) %>% 
    dplyr::summarise_all(sd, na.rm = TRUE)
  
  #Estimate statistical differences among null model treatments
  
  if (z == 2){
    # Use the rstatix package to implement paired t-test with bonferroni correction for each 
    # metric in eval.stats
    pairwise.nulls.ls <- list()
    
    for(i in eval.stats) {
      pairwise.mod <- as.formula(paste(eval.stats.i, "treatment.dif", sep = "~"))
      
      pairwise.nulls.ls[[i]] <- null.results.diff.comb %>%
        rstatix::pairwise_t_test(pairwise.mod, paired = TRUE, 
                                 alternative = alternative)
      
      emp.diff.comb <- emp.diff.comb %>% 
        tidyr::pivot_longer(cols = paste0(eval.stats, ".avg"), names_to = "metric", values_to = "emp.diff")
      nulls.diff.avg <- nulls.diff.avg %>% 
        tidyr::pivot_longer(cols = paste0(eval.stats, ".avg"), names_to = "metric", values_to = "null.diff.mean")
      nulls.diff.sd <- nulls.diff.sd %>% 
        tidyr::pivot_longer(cols = paste0(eval.stats, ".avg"), names_to = "metric", values_to = "null.diff.sd")
      
      empNull.stats <- dplyr::inner_join(emp.diff.comb, nulls.diff.avg, by = c("treatment.dif", "metric")) %>%
        dplyr::inner_join(nulls.diff.sd, by = c("treatment.dif", "metric")) %>%
        dplyr::group_by(treatment.dif, metric) %>%
        dplyr::mutate(zscore =((emp.diff - null.diff.mean)/ null.diff.sd))
      
      # find statistics that need a positive pnorm, and those that need a negative pnorm
      p.pos <- names(signs[sapply(signs, function(x) x == 1)]) %>% paste0(., ".avg")
      p.neg <- names(signs[sapply(signs, function(x) x == -1)]) %>% paste0(., ".avg")
      
      # obtain p-values
      empNull.stats[empNull.stats$metric %in% p.pos, "pvalue"] <- pnorm(empNull.stats[empNull.stats$metric %in% p.pos, "zscore", drop = TRUE], lower.tail = FALSE)
      empNull.stats[empNull.stats$metric %in% p.neg, "pvalue"] <- pnorm(empNull.stats[empNull.stats$metric %in% p.neg, "zscore", drop = TRUE])
    }
  } else if(z >= 3) {
    # Use the rstatix package to implement repeated measures anova for each 
    # metric in eval.stats
    anova.nulls.ls <- list()
    pairwise.nulls.ls <- list()
    for(i in eval.stats) {
      eval.stats.i <- paste0(i, ".avg")
      anova.nulls.ls[[i]] <- rstatix::anova_test(data  = null.results.diff.comb, 
                                                 dv = eval.stats.i,
                                                 wid = iter, 
                                                 within = treatment.dif)
      
      #post-hoc tests to examine pairwise differences among predictor sets
      pairwise.mod <- as.formula(paste(eval.stats.i, "treatment.dif", sep = "~"))
      
      pairwise.nulls.ls[[i]] <- null.results.diff.comb %>%
        rstatix::pairwise_t_test(pairwise.mod, paired = TRUE, 
                                 p.adjust.method = "bonferroni", 
                                 alternative = alternative)
      
      
      #Testing if empirical evaluation metrics are significantly different from null differences.
      #Calculate z-score and p-values for empirical against null treatment differences 
      if(anova.nulls.ls[[i]]$ANOVA$p < 0.05){
        #create output dataframe
        emp.diff.comb <- emp.diff.comb %>% 
          tidyr::pivot_longer(cols = paste0(eval.stats, ".avg"), names_to = "metric", values_to = "emp.diff")
        nulls.diff.avg <- nulls.diff.avg %>% 
          tidyr::pivot_longer(cols = paste0(eval.stats, ".avg"), names_to = "metric", values_to = "null.diff.mean")
        nulls.diff.sd <- nulls.diff.sd %>% 
          tidyr::pivot_longer(cols = paste0(eval.stats, ".avg"), names_to = "metric", values_to = "null.diff.sd")
        
        empNull.stats <- dplyr::inner_join(emp.diff.comb, nulls.diff.avg, by = c("treatment.dif", "metric")) %>%
          dplyr::inner_join(nulls.diff.sd, by = c("treatment.dif", "metric")) %>%
          dplyr::group_by(treatment.dif, metric) %>%
          dplyr::mutate(zscore =((emp.diff - null.diff.mean)/ null.diff.sd))
        
      } else if (anova.nulls.ls[[i]]$ANOVA$p >= 0.05){
        emp.diff.comb <- emp.diff.comb %>% 
          tidyr::pivot_longer(cols = paste0(eval.stats, ".avg"), names_to = "metric", values_to = "emp.diff")
        nulls.diff.avg.pool <- nulls.diff.avg.pool %>% 
          tidyr::pivot_longer(cols = paste0(eval.stats, ".avg"), names_to = "metric", values_to = "null.diff.mean.pool")
        nulls.diff.sd.pool <- nulls.diff.sd.pool %>% 
          tidyr::pivot_longer(cols = paste0(eval.stats, ".avg"), names_to = "metric", values_to = "null.diff.sd.pool")
        
        empNull.stats <- dplyr::inner_join(emp.diff.comb, nulls.diff.avg.pool, by = "metric")%>%
          dplyr::inner_join(nulls.diff.sd.pool, by = "metric") %>%
          dplyr::mutate(zscore =((emp.diff - null.diff.mean.pool)/ null.diff.sd.pool))
      }
    }
  }
  pairwise.nulls <- dplyr::bind_rows(pairwise.nulls.ls)
  
  
  # find statistics that need a positive pnorm, and those that need a negative pnorm
  p.pos <- names(signs[sapply(signs, function(x) x == 1)]) %>% paste0(., ".avg")
  p.neg <- names(signs[sapply(signs, function(x) x == -1)]) %>% paste0(., ".avg")
  
  # obtain p-values
  empNull.stats[empNull.stats$metric %in% p.pos, "pvalue"] <- pnorm(empNull.stats[empNull.stats$metric %in% p.pos, "zscore", drop = TRUE], lower.tail = FALSE)
  empNull.stats[empNull.stats$metric %in% p.neg, "pvalue"] <- pnorm(empNull.stats[empNull.stats$metric %in% p.neg, "zscore", drop = TRUE])
  
  # apply Bonferroni corrections to adjust p-values
  empNull.stats <- empNull.stats %>% mutate(p.adj = stats::p.adjust(pvalue, "bonferroni", n = z))
  
  if (z == 2 ){
    return(list(pairwise.nulls = pairwise.nulls, emp.nulls =empNull.stats, 
                nulls.diff = null.results.diff.comb))
    
  } else if (z >= 3){
    return(list(anova.nulls = anova.nulls.ls, emp.nulls = empNull.stats, pair.nulls = pairwise.nulls, 
                nulls.diff = null.results.diff.comb))  
  }
  
}
