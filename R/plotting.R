#' @title Partition group plots
#' @description Plot partition groups over an environmental predictor raster
#' @param e ENMevaluation object
#' @param envs Raster of an environmental predictor variable used to build the models in "e"
#' @param pts.type character specifying which to plot: occurrences ("occs") or background ("bg"), with default "occs"
#' @export

plot.grps <- function(e = NULL, pts = NULL, pts.grp = NULL, envs, pts.type = "occs") {
  if(!is.null(e)) {
    pts.plot <- switch(pts.type, occs = cbind(e@occ.pts, grp = e@occ.grp),
                  bg = cbind(e@bg.pts, grp = e@bg.grp))  
  }else{
    if(!is.null(pts) & !is.null(pts.grp)) {
      # make sure pts is a data frame with the right column names
      pts <- as.data.frame(pts)
      names(pts) <- c("longitude", "latitude")
      pts.plot <- cbind(pts, grp = factor(pts.grp))
    }else{
      stop("If inputting point data and not an ENMevaluation object, make sure to also input the partition groups (pts.grp).")
    }
  }
  
  if(length(unique(pts.plot$grp)) > 10) {
    theme.custom <- ggplot2::guides(color = FALSE)
  }else{
    theme.custom <- NULL
  }
  
  envs.df <- raster::as.data.frame(envs, xy = TRUE)
  names(envs.df)[3] <- "value"
  ggplot2::ggplot() + ggplot2::geom_raster(data = envs.df, ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::geom_point(data = pts.plot, ggplot2::aes(x = longitude, y = latitude, color = grp)) +
    ggplot2::scale_fill_distiller(palette = "Greys", na.value = "white") + ggplot2::theme_classic() + theme.custom
}

#' @title MESS plots for partition groups
#' @description Plot multivariate environmental similarity surface (MESS) estimates of test occurrences on the study extent 
#' by partition group with density curves or rasters
#' @param e ENMevaluation object
#' @param envs RasterStack of environmental predictor variables used to build the models in "e"; categorical variables should be 
#' removed before input, as they cannot be used to calculate MESS
#' @param pts.type character specifying which to calculate MESS on: occurrences ("occs") or background ("bg"), with default "occs"
#' @param plot.type character specifying which to plot: MESS density curves ("density") or MESS rasters ("raster")
#' @details As implemented here, MESS calculates the similarity between environmental values associated with the 
#' test occurrences (per partition group) and those associated with the entire study extent (specified by the extent 
#' of the input RasterStack "envs"). Higher negative values indicate greater environmental difference between the test occurrences
#' and the study extent, and higher positive values indicate greater similarity. This function uses the `dismo::mess()` function 
#' to calculate MESS. See the below reference for more details.
#' @return If "density", density curves showing the MESS estimates for each partition group. If "raster", rasters 
#' showing the geographical MESS estimates for each partition group.
#' @references Elith J., M. Kearney M., and S. Phillips, 2010. The art of modelling range-shifting species. Methods in Ecology and Evolution 1:330-342.
#' @export

plot.grps.mess <- function(e, envs, pts.type = "occs", plot.type = "density") {
  pts <- switch(pts.type, occs = dplyr::bind_cols(e@occ.pts, grp = e@occ.grp),
                bg = dplyr::bind_cols(e@bg.pts, grp = e@bg.grp))
  pts.x <- raster::extract(envs, pts[,1:2])
  vals <- data.frame(pts.x, grp = pts[,3]) 
  test.mss <- list()
  ras.mss <- list()
  nk <- length(unique(e@occ.grp))
  for(k in 1:nk) {
    test.vals <- vals %>% dplyr::filter(grp == k) %>% dplyr::select(-grp)
    train.xy <- pts %>% dplyr::filter(grp != k) %>% dplyr::select(-grp)
    # test.ext <- as(raster::extent(sp::bbox(sp::SpatialPoints(test.xy))), "SpatialPolygons")
    # envs.mess.train <- raster::mask(envs.mess, test.ext, inverse = TRUE)
    mss <- dismo::mess(envs, test.vals)
    ras.mss[[k]] <- mss
    mss.x <- raster::extract(mss, train.xy)
    test.mss[[k]] <- data.frame(mess.value = mss.x, grp = k)
  }
  if(plot.type == "raster") {
    rs.mss <- raster::stack(ras.mss)
    rs.mss.df <- raster::as.data.frame(rs.mss, xy = TRUE) %>% 
      tidyr::pivot_longer(cols = 3:ncol(.), names_to = "ras", values_to = "mess.value")
    rs.mss.df$ras <- unique(gsub("mess", "grp", rs.mss.df$ras))
    ggplot2::ggplot() + ggplot2::geom_raster(data = rs.mss.df, ggplot2::aes(x = x, y = y, fill = mess.value)) +
      ggplot2::geom_point(data = pts, ggplot2::aes(x = longitude, y = latitude, color = grp)) +
      ggplot2::facet_wrap(ggplot2::vars(ras), ncol = 2) +
      ggplot2::scale_fill_viridis_c(na.value = "white") + ggplot2::theme_classic()
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
#' @param stats vector of names of statistics from results table to be plotted; if more than
#' one statistic is specified, the plot will be faceted
#' @param x character of variable to be plotted on x-axis
#' @param col character of variable used to assign symbology colors 
#' @return A plot of evaluation statistics, with x representing a tuning parameter range, y
#' representing the average of a statistic over all partitions, and colors representing another
#' tuning parameter's values. Error bars represent the standard deviation of a statistic around the 
#' mean. Currently, this function only can handle two tuning parameters at a time.
#' @export

plot.eval <- function(e, stats, x, col) {
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
  res2 <- dplyr::left_join(avgs, sds) %>%
    dplyr::mutate(lower = avg - sd, upper = avg + sd)
  
  ggplot2::ggplot(res2, ggplot2::aes_string(x = x, y = "avg", color = col, group = col)) + 
    ggplot2::geom_point() + ggplot2::geom_line() + 
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper), width = 0.1) +
    ggplot2::facet_wrap(ggplot2::vars(metric), scales = "free_y", nrow = length(stats)) + ggplot2::theme_bw()
}
