#' @title Partition group plots
#' @description Plot occurrence partition groups over an environmental predictor raster.
#' @param e ENMevaluation object
#' @param envs Raster of an environmental predictor variable used to build the models in "e"
#' @param pts Matrix or data frame of coordinates for occurrence or background data
#' @param pts.grp Numeric vector of partition groups corresponding to data in "pts"
#' @param pts.type Character specifying which to plot: occurrences ("occs") or background ("bg"), with default "occs"
#' @param pts.size Character specifying which to plot: occurrences ("occs") or background ("bg"), with default "occs"
#' @details This function serves as a quick way to visualize occurrence or background partitions over the extent of an environmental predictor raster.
#' It can be run with an existing ENMevaluate object, or alternatively with occurrence or background coordinates and the corresponding partitions.
#' @export

evalplot.grps <- function(e = NULL, envs, pts = NULL, pts.grp = NULL, pts.type = "occs", pts.size = 1.5) {
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
  
  if(raster::nlayers(envs) > 1) {
    message("Plotting first raster in stack...")
    envs <- envs[[1]]
  }
  envs.df <- raster::as.data.frame(envs, xy = TRUE)
  names(envs.df)[3] <- "value"
  g <- ggplot2::ggplot() + ggplot2::geom_raster(data = envs.df, ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::geom_point(data = pts.plot, ggplot2::aes(x = longitude, y = latitude, color = partition), size = pts.size) +
    ggplot2::scale_color_manual(values = pt.cols) +
    ggplot2::scale_fill_distiller(palette = "Greys", na.value = "white") + ggplot2::theme_classic() + 
    ggplot2::coord_equal() + theme.custom
  return(g)
}

#' @title MESS plots for partition groups
#' @description Plot environmental similarity of predictor variable extent with respect to occurrence partitions
#' @param e ENMevaluation object
#' @param envs RasterStack: environmental predictor variables used to build the models in "e"; categorical variables should be 
#' removed before input, as they cannot be used to calculate MESS
#' @param envs.var Character: the name of a predictor variable to plot similarities for; if left NULL, calculations are done
#' with respect to all variables
#' @param sim.type Character: either "mess" for Multivariate Environmental Similarity Surface, "most_diff" for most different variable,
#' or "most_sim" for most similar variable; uses similarity function from package rmaxent
#' @param pts.type Character: the reference to calculate MESS based on: occurrences ("occs") or background ("bg"), with default "occs"
#' @param plot.type Character specifying which to plot: MESS histograms ("histogram") or MESS rasters ("raster"); default is histogram
#' @param sim.palette Character: the RColorBrewer palette name to use for plotting; if left NULL, the defaults are
#' "Set1" for discrete variables and reverse "RdYlBu" for continuous variables
#' @param palette.dir Numeric: direction for palette plotting; either 1 or -1
#' @param hist.bins Numeric: number of histogram bins for histogram plots; default is 30
#' @details There are two variations for this plot. If "histogram", histograms are plotted showing the MESS estimates for each partition group. 
#' If "raster", rasters are plotted showing the geographical MESS estimates for each partition group. 
#' With sim.type option "mess", the similarity between environmental values associated with the 
#' validation occurrences (per partition group) and those associated with the entire study extent (specified by the extent 
#' of the input RasterStack "envs") are calculated, and the minimum similarity per grid is returned. 
#' Higher negative values indicate greater environmental difference between the validation occurrences
#' and the study extent, and higher positive values indicate greater similarity. This function uses the `rmaxent::similarity()` function 
#' to calculate the similarities. See the below reference for details on MESS. 
#' @return A ggplot of MESS calculations for data partitions.
#' @references Elith J., M. Kearney M., and S. Phillips, 2010. The art of modelling range-shifting species. Methods in Ecology and Evolution 1:330-342.
#' @export

evalplot.envSim.hist <- function(e = NULL, envs = NULL, occs = NULL, bg = NULL, occs.grp = NULL, 
                                 bg.grp = NULL, envs.var = NULL, pts.type = c("occs", "bg"), 
                                 sim.type = c("mess", "most_diff", "most_sim"),
                                 hist.bins = 30, return.tbl = FALSE) {
  # remove categorical rasters
  nr <- raster::nlayers(envs)
  cats <- numeric(nr)
  for(n in 1:nr) cats[n] <- is.factor(envs[[n]])
  if(sum(cats) > 0) {
    rem.ras <- which(cats == 1)
    message(paste("Ignoring categorical raster", names(envs)[rem.ras], "..."))
    envs <- envs[[-rem.ras]]
  }
  
  if(!is.null(e)) {
    pts <- rbind(e@occs, e@bg) %>% mutate(type = c(rep(1, nrow(e@occs)), rep(0, nrow(e@bg))),
                                          partition = c(e@occs.grp, e@bg.grp))
  }else{
    pts <- switch(pts.type, occs = occs, bg = bg)
    pts.grp <- switch(pts.type, occs = occs.grp, bg = bg.grp)
    pts.z <- raster::extract(envs, pts)
    pts.plot <- bind_cols(pts, pts.z, partition = factor(pts.grp))
  }
  names(pts.plot)[1:2] <- c("longitude","latitude")
  
  test.sim <- list()
  nk <- length(unique(pts.grp))
  
  for(k in 1:nk) {
    test.z <- pts.plot %>% dplyr::filter(partition == k) %>% dplyr::select(-longitude, -latitude, -partition)
    train.z <- pts.plot %>% dplyr::filter(partition != k) %>% dplyr::select(-longitude, -latitude, -partition)
    if(!is.null(envs.var)) {
      test.z <- test.z %>% select(all_of(envs.var))
      train.z <- train.z %>% select(all_of(envs.var))
    }
    sim <- rmaxent::similarity(train.z, test.z)
    mss <- switch(sim.type, mess = sim$similarity_min, most_diff = sim$mod, most_sim = sim$mos)  
    
    test.sim[[k]] <- data.frame(mss, partition = k)
    names(test.sim[[k]])[1] <- sim.type  
  }
  
  if(nk > 9) {
    theme.custom <- ggplot2::guides(color = FALSE)
    pt.cols <- rainbow(grp.n)
  }else{
    theme.custom <- NULL
    pt.cols <- RColorBrewer::brewer.pal(nk, "Set1")
  }
  
  plot.df <- dplyr::bind_rows(test.sim)
  plot.df$partition <- factor(plot.df$partition)
  if(sim.type != "mess") {
    envs.tbl <- data.frame(sort(unique(plot.df[,1])), names(envs))
    names(envs.tbl) <- c(sim.type, "env.var")
    envs.tbl$env.var <- factor(envs.tbl$env.var)
    plot.df <- plot.df %>% left_join(envs.tbl, by = sim.type)
    title <- ifelse(sim.type == "most_diff", "Most different environmental variable\n(Differences based on training groups with respect to validation group)", 
                    "Most similar environmental variable\n(Differences based on training groups with respect to validation group)")
    g <- ggplot2::ggplot(plot.df, ggplot2::aes(x = env.var, fill = partition)) + 
      ggplot2::stat_count() +
      ggplot2::facet_grid(ggplot2::vars(partition)) + 
      ggplot2::scale_fill_manual(values = pt.cols) +
      ggplot2::theme_classic() +
      ggplot2::geom_vline(xintercept = 0) +
      ggplot2::ggtitle(title) +
      ggplot2::theme(
        strip.background = ggplot2::element_blank(),
        strip.text.y = ggplot2::element_blank()
      )
  }else{
    g <- ggplot2::ggplot(plot.df, ggplot2::aes(x = mess, fill = partition)) + 
      ggplot2::geom_histogram(bins = hist.bins) +
      ggplot2::facet_grid(ggplot2::vars(partition)) + 
      ggplot2::scale_fill_manual(values = pt.cols) +
      ggplot2::theme_classic() +
      ggplot2::geom_vline(xintercept = 0) +
      ggplot2::ggtitle("Multivariate environmental similarity\n(Differences based on training groups with respect to validation group)") +
      ggplot2::theme(
        strip.background = ggplot2::element_blank(),
        strip.text.y = ggplot2::element_blank()
      )  
  }
  g
  
  if(return.tbl == TRUE) {
    return(plot.df)
  }else{
    return(g)  
  }
}

evalplot.envSim.map <- function(e = NULL, envs = NULL, occs = NULL, bg = NULL, occs.grp = NULL, 
                                bg.grp = NULL, envs.var = NULL, pts.type = c("occs", "bg"), 
                                pts.size = 1.5, mess.cols = c("red","yellow","blue"), mess.naCol = "white",
                                sim.palette = NULL,
                                sim.type = c("mess", "most_diff", "most_sim"), bb.buf = NULL,
                                return.tbl = FALSE) {
  # remove categorical rasters
  nr <- raster::nlayers(envs)
  cats <- numeric(nr)
  for(n in 1:nr) cats[n] <- is.factor(envs[[n]])
  if(sum(cats) > 0) {
    rem.ras <- which(cats == 1)
    message(paste("Ignoring categorical raster", names(envs)[rem.ras], "..."))
    envs <- envs[[-rem.ras]]
  }
  
  if(!is.null(e)) {
    pts <- rbind(e@occs, e@bg) %>% mutate(type = c(rep(1, nrow(e@occs)), rep(0, nrow(e@bg))),
                                          partition = c(e@occs.grp, e@bg.grp))
  }else{
    pts <- switch(pts.type, occs = occs, bg = bg)
    pts.grp <- switch(pts.type, occs = occs.grp, bg = bg.grp)
    pts.z <- raster::extract(envs, pts)
    pts.plot <- bind_cols(pts, pts.z, partition = factor(pts.grp))
  }
  names(pts.plot)[1:2] <- c("longitude","latitude")
  # pts$id <- row.names(pts)
  
  ras.sim <- list()
  nk <- length(unique(pts.grp))
  
  for(k in 1:nk) {
    test.z <- pts.plot %>% dplyr::filter(partition == k) %>% dplyr::select(-longitude, -latitude, -partition)
    if(!is.null(envs.var)) {
      test.z <- test.z %>% select(all_of(envs.var))
    }
    sim <- rmaxent::similarity(envs, test.z)
    sim.sel <- switch(sim.type, mess = sim$similarity_min, most_diff = sim$mod, most_sim = sim$mos)  
    
    ras.sim[[k]] <- sim.sel
  }
  
  rs.sim <- raster::stack(ras.sim)
  names(rs.sim) <- gsub("layer|mess", "partition", names(rs.sim))
  plot.df <- raster::as.data.frame(rs.sim, xy = TRUE) %>%
    tidyr::pivot_longer(cols = 3:ncol(.), names_to = "ras", values_to = "value")
  if(!is.null(bb.buf)) {
    plot.df <- plot.df %>% filter(x > min(pts$longitude) - bb.buf, x < max(pts$longitude) + bb.buf,
                                  y > min(pts$latitude) - bb.buf, y < max(pts$latitude) + bb.buf)
  }
  if(sim.type != "mess") {
    if(is.null(sim.palette)) sim.palette <- "Set1"
    plot.df$ras <- gsub("_var", "", plot.df$ras)
    plot.df$value <- factor(plot.df$value)
    title.part <- switch(sim.type, most_diff = "Most different", most_sim = "Most similar")
    title <- paste(title.part, "environmental variable\n(Differences based on entire extent with respect to validation group)")
  }
  if(sim.type == "mess") title <- "Multivariate environmental similarity\n(Differences based on entire extent with respect to validation group)"
  
  g <- ggplot2::ggplot() +
    ggplot2::geom_raster(data = plot.df, ggplot2::aes(x = x, y = y, fill = value)) +
    # geom_sf(data=west.shp2,fill="gray40",color=NA) +
    ggplot2::geom_point(data = pts.plot, ggplot2::aes(x = longitude, y = latitude, shape = partition), color = "black", size = pts.size) +
    ggplot2::facet_wrap(ggplot2::vars(ras), ncol = 2) +
    ggplot2::ggtitle(title) +
    ggplot2::theme_classic()
  if(sim.type != "mess") {
    g <- g + ggplot2::scale_fill_brewer(palette = sim.palette, na.value = "white", breaks = levels(plot.df$value))
  }else{
    g <- g + ggplot2::scale_fill_gradientn(na.value = mess.naCol,
                                           colors = mess.cols,
                                           values = scales::rescale(c(min(plot.df$value, na.rm=TRUE), 0, 0.01,
                                                                      max(plot.df$value, na.rm=TRUE))))
  }
  g
  
  if(return.tbl == TRUE) {
    return(plot.df)
  }else{
    return(g)  
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

evalplot.stats <- function(e, stats, x, col, dodge = NULL, error.bars = TRUE, facet.labs = NULL, metric.levs = NULL) {
  exp <- paste(paste0("*", stats), collapse = "|")
  res <- e@results %>% 
    tidyr::pivot_longer(cols = auc.train:ncoef, names_to = "metric", values_to = "value") %>%
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
  if(!is.null(facet.labs)) labeller <- as_labeller(facet.labs) else labeller <- NULL
  if(!is.null(metric.levs)) res$metric <- factor(res$metric, levels = metric.levs)
  if(!is.null(metric.levs)) res.avgs$metric <- factor(res.avgs$metric, levels = metric.levs)
  
  if(nrow(res.avgs) > 0) {
    if(is.null(dodge)) dodge <- 0.2
    p <- ggplot2::ggplot(res.avgs, ggplot2::aes_string(x = x, y = "avg", color = col, group = col)) + 
      ggplot2::geom_point(position=ggplot2::position_dodge(width=dodge)) + 
      ggplot2::geom_line(position=ggplot2::position_dodge(width=dodge)) + 
      ggplot2::facet_wrap(ggplot2::vars(metric), scales = "free_y", nrow = length(stats), labeller = labeller) + 
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
      ggplot2::facet_wrap(ggplot2::vars(metric), scales = "free_y", nrow = length(stats), labeller = labeller) + 
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

evalplot.nulls <- function(e.null, stats, plot.type, facet.labs = NULL, metric.levs = NULL) {
  exp <- paste(paste0("*", stats), collapse = "|")
  null.res <- e.null@null.results %>% 
    tidyr::pivot_longer(cols = auc.train:ncoef, names_to = "metric", values_to = "value") %>%
    dplyr::filter(grepl(exp, metric)) %>%
    dplyr::select(metric, value)
  null.avgs <- null.res %>% 
    dplyr::filter(grepl("avg", metric) | metric %in% stats) %>%
    dplyr::rename(avg = value) %>%
    dplyr::mutate(metric = gsub(".avg", "", metric))
  if(!is.null(metric.levs)) null.avgs$metric <- factor(null.avgs$metric, levels = metric.levs)
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
  
  if(!is.null(facet.labs)) labeller <- as_labeller(facet.labs) else labeller <- NULL
  
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
    vlines <- rbind(vlines, real.res %>% dplyr::mutate(quantile = "empirical value") %>% dplyr::rename(value = avg))
    ggplot2::ggplot(mapping = ggplot2::aes(x = avg)) + 
      ggplot2::geom_histogram(data = null.avgs, fill = "gray80") +
      ggplot2::geom_vline(data = vlines, ggplot2::aes(xintercept = value, color = quantile, linetype = quantile)) +
      ggplot2::scale_color_manual(values = c(`0.01 quantile` = "purple", 
                                             `0.05 quantile` = "blue",
                                             `0.50 quantile` = "blue",
                                             `0.95 quantile` = "blue",
                                             `0.99 quantile` = "purple",
                                             `empirical value` = "red")) +
      ggplot2::scale_linetype_manual(values = c(`0.01 quantile` = "dotted", 
                                                `0.05 quantile` = "dashed",
                                                `0.50 quantile` = "solid",
                                                `0.95 quantile` = "dashed", 
                                                `0.99 quantile` = "dotted",
                                                `empirical value` = "solid")) +
      ggplot2::facet_wrap(ggplot2::vars(metric), scales = "free_x", ncol = 1, labeller = labeller) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.title=ggplot2::element_blank(), 
                     legend.position="bottom")
  }
}
