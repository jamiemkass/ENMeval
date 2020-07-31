#' @title Partition group plots
#' @description Plot occurrence partition groups over an environmental predictor raster.
#' @param e ENMevaluation object
#' @param envs Raster of an environmental predictor variable used to build the models in "e"
#' @param pts Matrix or data frame of coordinates for occurrence or background data
#' @param pts.grp Numeric vector of partition groups corresponding to data in "pts"
#' @param pts.type Character specifying which to plot: occurrences ("occs") or background ("bg"), with default "occs"
#' @details This function serves as a quick way to visualize occurrence or background partitions over the extent of an environmental predictor raster.
#' It can be run with an existing ENMevaluate object, or alternatively with occurrence or background coordinates and the corresponding partitions.
#' @export

evalplot.grps <- function(e = NULL, envs, pts = NULL, pts.grp = NULL, pts.type = "occs") {
  if(!is.null(e)) {
    pts.plot <- switch(pts.type, occs = cbind(e@occs, partition = e@occs.grp),
                  bg = cbind(e@bg, partition = e@bg.grp))  
    names(pts.plot)[1:2] <- c("longitude", "latitude")
  }else{
    if(!is.null(pts) & !is.null(pts.grp)) {
      # make sure pts is a data frame with the right column names
      pts <- as.data.frame(pts)
      names(pts) <- c("longitude", "latitude")
      pts.plot <- cbind(pts, partition = factor(pts.grp))
    }else{
      stop("If inputting point data instead of an ENMevaluation object, make sure to also input the partition groups (pts.grp).")
    }
  }
  
  grp.n <- length(unique(pts.plot$partition))
  if(grp.n > 9) {
    theme.custom <- ggplot2::guides(color = FALSE)
    pt.cols <- rainbow(grp.n)
  }else{
    theme.custom <- NULL
    pt.cols <- RColorBrewer::brewer.pal(9, "Set1")
  }
  
  envs.df <- raster::as.data.frame(envs, xy = TRUE)
  names(envs.df)[3] <- "value"
  ggplot2::ggplot() + ggplot2::geom_raster(data = envs.df, ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::geom_point(data = pts.plot, ggplot2::aes(x = longitude, y = latitude, color = partition)) +
    ggplot2::scale_color_manual(values = pt.cols) +
    ggplot2::scale_fill_distiller(palette = "Greys", na.value = "white") + ggplot2::theme_classic() + 
    ggplot2::coord_equal() + theme.custom
}

#' @title MESS plots for partition groups
#' @description Plot multivariate environmental similarity surface (MESS) estimates of test occurrences on the study extent 
#' by partition group with density curves or rasters
#' @param e ENMevaluation object
#' @param envs RasterStack of environmental predictor variables used to build the models in "e"; categorical variables should be 
#' removed before input, as they cannot be used to calculate MESS
#' @param pts.type Character specifying which to calculate MESS based on: occurrences ("occs") or background ("bg"), with default "occs"
#' @param plot.type Character specifying which to plot: MESS density curves ("density") or MESS rasters ("raster")
#' @details There are two variations for this plot. If "density", density curves are plotted showing the MESS estimates for each partition group. 
#' If "raster", rasters are plotted showing the geographical MESS estimates for each partition group. 
#' As implemented here, MESS calculates the similarity between environmental values associated with the 
#' test occurrences (per partition group) and those associated with the entire study extent (specified by the extent 
#' of the input RasterStack "envs"). Higher negative values indicate greater environmental difference between the test occurrences
#' and the study extent, and higher positive values indicate greater similarity. This function uses the `dismo::mess()` function 
#' to calculate MESS. See the below reference for more details. 
#' @return A ggplot of MESS calculations for data partitions.
#' @references Elith J., M. Kearney M., and S. Phillips, 2010. The art of modelling range-shifting species. Methods in Ecology and Evolution 1:330-342.
#' @export

evalplot.grps.mess <- function(e, envs, pts.type = "occs", plot.type = "density") {
  names(e@occs)[1:2] <- c("longitude","latitude")
  pts <- switch(pts.type, occs = dplyr::bind_cols(e@occs[,c("longitude","latitude")], grp = e@occs.grp),
                bg = dplyr::bind_cols(e@bg[,c("longitude","latitude")], grp = e@bg.grp))
  pts.x <- raster::extract(envs, pts[,c("longitude","latitude")])
  vals <- data.frame(pts.x, grp = pts$grp) 
  test.mss <- list()
  ras.mss <- list()
  nk <- length(unique(e@occs.grp))
  for(k in 1:nk) {
    test.vals <- vals %>% dplyr::filter(grp == k) %>% dplyr::select(-grp)
    train.xy <- pts %>% dplyr::filter(grp != k) %>% dplyr::select(-grp)
    # test.ext <- as(raster::extent(sp::bbox(sp::SpatialPoints(test.xy))), "SpatialPolygons")
    # envs.mess.train <- raster::mask(envs.mess, test.ext, inverse = TRUE)
    mss <- enmeval.mess(envs, test.vals)
    ras.mss[[k]] <- mss
    mss.x <- raster::extract(mss, train.xy)
    test.mss[[k]] <- data.frame(mess.value = mss.x, grp = k)
  }
  if(plot.type == "raster") {
    rs.mss <- raster::stack(ras.mss)
    rs.mss.df <- raster::as.data.frame(rs.mss, xy = TRUE) %>% 
      tidyr::pivot_longer(cols = 3:ncol(.), names_to = "ras", values_to = "mess.value")
    rs.mss.df$ras <- gsub("mess", "grp", rs.mss.df$ras)
    ggplot2::ggplot() + ggplot2::geom_raster(data = rs.mss.df, ggplot2::aes(x = x, y = y, fill = mess.value)) +
      ggplot2::scale_fill_viridis_c(na.value = "white") +
      ggplot2::geom_point(data = pts, ggplot2::aes(x = longitude, y = latitude, shape = grp)) +
      ggplot2::facet_wrap(ggplot2::vars(ras), ncol = 2) +
      ggplot2::theme_classic()
  }else{
    test.mss.df <- dplyr::bind_rows(test.mss)
    test.mss.df$grp <- factor(test.mss.df$grp)
    ggplot2::ggplot(test.mss.df, ggplot2::aes(x = mess.value, fill = grp)) + 
      ggplot2::geom_density() +
      ggplot2::facet_grid(ggplot2::vars(grp)) + 
      ggplot2::theme_classic() +
      ggplot2::theme(
        strip.background = ggplot2::element_blank(),
        strip.text.y = ggplot2::element_blank()
      )  
  }
}
#' @title ENMevaluation statistics plot
#' @description Plot evaluation statistics over tuning parameter ranges to visualize differences in performance
#' @param e ENMevaluation object
#' @param stats Character vector of names of statistics from results table to be plotted; if more than
#' one statistic is specified, the plot will be faceted
#' @param x Character of variable to be plotted on x-axis
#' @param col Character of variable used to assign symbology colors 
#' @param dodge (Optional) numeric value specifying the dodge width for points and lines (this improves visibility when there is high overlap) 
#' @param error.bars (Optional) boolean specifying whether or not to plot error bars (defaults to TRUE)
#' @details In this plot, the x-axis represents a tuning parameter range, while the y-axis represents the average of a statistic over all partitions.
#' Different colors represent the categories or values of another tuning parameter. 
#' Error bars represent the standard deviation of a statistic around the mean. 
#' Currently, this function can only plot two tuning parameters at a time.
#' @return A ggplot of evaluation statistics. 
#' @export

evalplot.stats <- function(e, stats, x, col, dodge = NULL, error.bars = TRUE) {
  exp <- paste(paste0("*", stats), collapse = "|")
  res <- e@results %>% 
    tidyr::pivot_longer(cols = auc.train:nparam, names_to = "metric", values_to = "value") %>%
    dplyr::filter(grepl(exp, metric))
  avgs <- res %>% 
    dplyr::filter(grepl("avg", metric)) %>%
    dplyr::rename(avg = value) %>%
    dplyr::mutate(metric = gsub(".avg", "", metric))
  sds <- res %>% 
    dplyr::filter(grepl("sd", metric)) %>%
    dplyr::rename(sd = value) %>%
    dplyr::mutate(metric = gsub(".sd", "", metric))
  join.names <- names(avgs)
  join.names <- join.names[join.names != "avg"]
  res.avgs <- dplyr::left_join(avgs, sds, by = join.names) %>%
    dplyr::mutate(lower = avg - sd, upper = avg + sd,
                  metric = factor(metric, levels = stats))
  
  if(nrow(res.avgs) > 0) {
    if(is.null(dodge)) dodge <- 0.2
    p <- ggplot2::ggplot(res.avgs, ggplot2::aes_string(x = x, y = "avg", color = col, group = col)) + 
      ggplot2::geom_point(position=ggplot2::position_dodge(width=dodge)) + 
      ggplot2::geom_line(position=ggplot2::position_dodge(width=dodge)) + 
      ggplot2::facet_wrap(ggplot2::vars(metric), scales = "free_y", nrow = length(stats)) + 
      ggplot2::theme_bw()  
    if(error.bars == TRUE) {
      p + ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper), width = 0.5, 
                                                      position=ggplot2::position_dodge(width=dodge))
    }else{
      p
    }
    
  }else{
    if(is.null(dodge)) dodge <- 0
    ggplot2::ggplot(res, ggplot2::aes_string(x = x, y = "value", color = col, group = col)) + 
      ggplot2::geom_point(position=ggplot2::position_dodge(width=dodge)) + 
      ggplot2::geom_line(position=ggplot2::position_dodge(width=dodge)) + 
      ggplot2::facet_wrap(ggplot2::vars(metric), scales = "free_y", nrow = length(stats)) + 
      ggplot2::theme_bw()
  }
}

#' @title ENMnullSims statistics plot
#' @description Plot histogram of evaluation statistics for null ENM simulations
#' @param e.null ENMnull object
#' @param stats Character vector of statistics from results table to be plotted; if more than
#' one statistic is specified, the histogram plot will be faceted
#' @param plot.type Character specifying the plot type: either "violin" or "histogram"
#' @details There are two variations for this plot, but both show null quantiles (0.01, 0.05, 0.5, 0.95, 0.99). 
#' For violin plots, the null distribution is displayed as a vertical shape (i.e., the violin) with horizontal lines showing 
#' the quantiles and the real value is plotted as a red point along the vertical axis. 
#' For histogram plots, the null distribution is displayed as a histogram with vertical lines showing the quantiles 
#' and the real value as a vertical red line on the distribution.
#' @return A ggplot of null model statistics. 
#' @export

evalplot.nulls <- function(e.null, stats, plot.type) {
  exp <- paste(paste0("*", stats), collapse = "|")
  null.res <- e.null@null.results %>% 
    tidyr::pivot_longer(cols = auc.train:nparam, names_to = "metric", values_to = "value") %>%
    dplyr::filter(grepl(exp, metric)) %>%
    dplyr::select(metric, value)
  null.avgs <- null.res %>% 
    dplyr::filter(grepl("avg", metric) | metric %in% stats) %>%
    dplyr::rename(avg = value) %>%
    dplyr::mutate(metric = gsub(".avg", "", metric))
  # null.sds <- null.res %>% 
  #   dplyr::filter(grepl("sd", metric)) %>%
  #   dplyr::rename(sd = value) %>%
  #   dplyr::mutate(metric = gsub(".sd", "", metric))
  # null.res.avgs <- dplyr::bind_cols(null.avgs, null.sds %>% dplyr::select(sd))
  
  real.res <- e.null@real.vs.null.results %>% 
    dplyr::slice(1) %>%
    tidyr::pivot_longer(cols = stats, names_to = "metric", values_to = "value") %>%
    dplyr::select(statistic, metric, value) %>%
    tidyr::pivot_wider(names_from = statistic, values_from = value) %>%
    dplyr::rename(avg = real.mean)
  
  # stat.nullResults.name <- ifelse(stat == "auc.train", "auc.train", paste0(stat, ".avg"))
  # null.stats <- round(e.null@null.results[,stat.nullResults.name, drop = FALSE], 3)
  # real.stat <- round(e.null@real.vs.null.results %>% 
                       # dplyr::filter(statistic == "real.mean") %>% 
                       # dplyr::pull(stat), 3)
  # stat.max <- ifelse(real.stat < max(null.stats), max(null.stats), real.stat)
  # all.stats <- rbind(null.stats, real.stat)
  # vlines <- data.frame(name = c("95 quantile", "99 quantile", "real value"),
                       # value = c(quantile(null.stats[[stat.nullResults.name]], 0.95),
                                 # quantile(null.stats[[stat.nullResults.name]], 0.99),
                                 # real.stat))
  if(plot.type == "violin") {
    ggplot2::ggplot(null.avgs, ggplot2::aes(x = metric, y = avg)) + 
      ggplot2::geom_violin(draw_quantiles = c(0.01, 0.05, 0.5, 0.95, 0.99)) +
      ggplot2::geom_point(data = real.res, ggplot2::aes(y = avg), color = "red") +
      ggplot2::theme_bw()  
  }else if(plot.type == "histogram") {
    stats.all <- rbind(null.avgs, real.res)
    vlines <- null.avgs %>% dplyr::group_by(metric) %>% 
      dplyr::summarize(`0.01 quantile` = quantile(avg, 0.01),
                       `0.05 quantile` = quantile(avg, 0.05),
                       `0.50 quantile` = quantile(avg, 0.5),
                       `0.95 quantile` = quantile(avg, 0.95),
                       `0.99 quantile` = quantile(avg, 0.99)) %>%
      tidyr::pivot_longer(cols = `0.01 quantile`:`0.99 quantile`, names_to = "quantile", values_to = "value")
    vlines <- rbind(vlines, real.res %>% dplyr::mutate(quantile = "real value") %>% dplyr::rename(value = avg))
    ggplot2::ggplot(mapping = ggplot2::aes(x = avg)) + 
      ggplot2::geom_histogram(data = null.avgs, fill = "gray80") +
      ggplot2::geom_vline(data = vlines, ggplot2::aes(xintercept = value, color = quantile, linetype = quantile)) +
      ggplot2::scale_color_manual(values = c(`0.01 quantile` = "purple", 
                                             `0.05 quantile` = "blue",
                                             `0.50 quantile` = "blue",
                                             `0.95 quantile` = "blue",
                                             `0.99 quantile` = "purple",
                                             `real value` = "red")) +
      ggplot2::scale_linetype_manual(values = c(`0.01 quantile` = "dotted", 
                                                `0.05 quantile` = "dashed",
                                                `0.50 quantile` = "solid",
                                                `0.95 quantile` = "dashed", 
                                                `0.99 quantile` = "dotted",
                                                `real value` = "solid")) +
      ggplot2::facet_wrap(ggplot2::vars(metric), scales = "free_x", ncol = 1) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.title=ggplot2::element_blank(), 
                     legend.position="bottom")
  }
}
