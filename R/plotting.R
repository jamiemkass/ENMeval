#' @export
plot.grps <- function(e, r, pts.type) {
  pts <- switch(pts.type, occs = dplyr::bind_cols(e@occ.pts, grp = e@occ.grp),
                bg = dplyr::bind_cols(e@bg.pts, grp = e@bg.grp))
  
  r.df <- raster::as.data.frame(r, xy = TRUE)
  names(r.df)[3] <- "value"
  ggplot2::ggplot() + ggplot2::geom_raster(data = r.df, ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::geom_point(data = pts, ggplot2::aes(x = longitude, y = latitude, color = grp)) +
    ggplot2::scale_fill_distiller(palette = "Greys", na.value = "white") + ggplot2::theme_classic()
}

#' @export
plot.grps.mess <- function(e, envs, categoricals, pts.type, plot.type = "density") {
  pts <- switch(pts.type, occs = dplyr::bind_cols(e@occ.pts, grp = e@occ.grp),
                bg = dplyr::bind_cols(e@bg.pts, grp = e@bg.grp))
  pts.x <- raster::extract(envs, pts[,1:2])
  vals <- data.frame(pts.x, grp = pts[,3]) %>% dplyr::select(-categoricals)
  i <- which(names(envs) == categoricals)
  envs.mess <- envs[[-i]]
  test.mss <- list()
  ras.mss <- list()
  nk <- length(unique(e@occ.grp))
  for(k in 1:nk) {
    test.vals <- vals %>% dplyr::filter(grp == k) %>% dplyr::select(-grp)
    train.xy <- pts %>% dplyr::filter(grp != k) %>% dplyr::select(-grp)
    # test.ext <- as(raster::extent(sp::bbox(sp::SpatialPoints(test.xy))), "SpatialPolygons")
    # envs.mess.train <- raster::mask(envs.mess, test.ext, inverse = TRUE)
    mss <- dismo::mess(envs.mess, test.vals)
    ras.mss[[k]] <- mss
    mss.x <- raster::extract(mss, train.xy)
    test.mss[[k]] <- data.frame(mess.value = mss.x, grp = k)
  }
  if(plot.type == "rasters") {
    rs.mss <- raster::stack(ras.mss)
    rs.mss.df <- raster::as.data.frame(rs.mss, xy = TRUE) %>% 
      tidyr::pivot_longer(cols = 3:ncol(rs.mss.df), names_to = "ras", values_to = "mess.value")
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

#' @export
plot.eval <- function(e, x, color, vars) {
  exp <- paste(paste0("*", vars), collapse = "|")
  res <- e@results %>% 
    tidyr::pivot_longer(cols = auc.train:nparam, names_to = "metric", values_to = "value") %>%
    dplyr::filter(grepl(exp, metric)) %>% 
    tidyr::separate(metric, sep = "_", fill = "right", c("metric", "statistic")) %>%
    tidyr::replace_na(list(statistic = "avg")) %>%
    tidyr::pivot_wider(names_from = statistic, values_from = value) %>%
    dplyr::mutate(lower = avg - sd, upper = avg + sd)
  
  ggplot2::ggplot(res, ggplot2::aes_string(x = x, y = "avg", color = color, group = color)) + 
    ggplot2::geom_point() + ggplot2::geom_line() + 
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper), width = 0.1) +
    ggplot2::facet_wrap(ggplot2::vars(metric), scales = "free_y", nrow = length(vars)) + ggplot2::theme_bw()
}
