#' @title Partition group plots
#' @description Plot occurrence partition groups over an environmental predictor raster.
#' @param e ENMevaluation object
#' @param envs SpatRaster: environmental predictor variable used to build the models in "e"
#' @param pts matrix / data frame: coordinates for occurrence or background data
#' @param pts.grp numeric vector: partition groups corresponding to data in "pts"
#' @param ref.data character: plot occurrences ("occs") or background ("bg"), with default "occs"
#' @param pts.size numeric: custom point size for ggplot
#' @param return.tbl boolean: if TRUE, return the data frames used to make the ggplot instead of the plot itself
#' @details This function serves as a quick way to visualize occurrence or background partitions over the extent of an environmental predictor raster.
#' It can be run with an existing ENMevaluation object, or alternatively with occurrence or background coordinates and the corresponding partitions.
#' @export

evalplot.grps <- function(e = NULL, envs, pts = NULL, pts.grp = NULL, ref.data = "occs", pts.size = 1.5, return.tbl = FALSE) {
  if(!is.null(e)) {
    pts.plot <- switch(ref.data, occs = cbind(e@occs, partition = e@occs.grp),
                       bg = cbind(e@bg, partition = e@bg.grp))  
    if(e@partition.method == "testing") {
      pts.plot <- pts.plot |> dplyr::mutate(partition = as.numeric(as.character(partition))) |> 
        dplyr::bind_rows(e@occs.testing |> dplyr::mutate(partition = 1)) |>
        dplyr::mutate(partition = factor(partition))
    } 
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
  
  if(terra::nlyr(envs) > 1) {
    message("Plotting first raster in stack...")
    envs <- envs[[1]]
  }
  envs.df <- terra::as.data.frame(envs, xy = TRUE, na.rm = FALSE)
  names(envs.df)[3] <- "value"
  g <- ggplot2::ggplot() + ggplot2::geom_raster(data = envs.df, ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::geom_point(data = pts.plot, ggplot2::aes(x = longitude, y = latitude, color = partition), size = pts.size) +
    ggplot2::scale_color_manual(values = pt.cols) +
    ggplot2::scale_fill_distiller(palette = "Greys", na.value = "white") + ggplot2::theme_classic() + 
    ggplot2::coord_equal() + theme.custom
  
  if(return.tbl == TRUE) {
    return(list(envs.df = envs.df, pts.plot = pts.plot))
  }else{
    return(g)  
  }
}




plot.sim.dataPrep <- function(e, envs, occs.z, bg.z, occs.grp, bg.grp, ref.data, categoricals, occs.testing.z, quiet) {
  
  if(!is.null(e) & any(!is.null(occs.z), !is.null(bg.z), !is.null(occs.grp), !is.null(bg.grp))) {
    stop("* If inputting an ENMevaluation object, leave occs.z, bg.z, occs.grp, and bg.grp NULL. These are read from the object.")
  }
  
  if(is.null(envs)) {
    if(is.null(e) & any(is.null(occs.z), is.null(bg.z), is.null(occs.grp), is.null(bg.grp))) {
      stop("* If inputting occurrence and background data instead of an ENMevaluation object, please input occs.z, bg.z, occs.grp, and bg.grp.")
      if(!quiet) message("* Similarity values calculated by contrasting occurrences with background.")
    }
  }else{
    if(is.null(e)) {
      if(ref.data == "occs") {
        if(any(is.null(occs.z), is.null(occs.grp))) {stop("* If inputting occurrence data instead of an ENMevaluation object, please input occs.z and occs.grp.")}
      }else if (ref.data == "bg") {
        if(any(is.null(bg.z), is.null(bg.grp))) {stop("* If inputting background data instead of an ENMevaluation object, please input bg.z and bg.grp.")}
      }
    }
    if(!quiet) message("* Similarity values calculated by contrasting occurrences with all cell values in raster extent.")
  }
  
  # assign variables from ENMevaluation object
  if(!is.null(e)) {
    occs.z <- e@occs
    bg.z <- e@bg
    occs.grp <- as.numeric(as.character(e@occs.grp))
    bg.grp <- as.numeric(as.character(e@bg.grp))
  }else{
    occs.grp <- as.numeric(as.character(occs.grp))
    bg.grp <- as.numeric(as.character(bg.grp))
  }
  
  if(ref.data == "bg" & length(unique(bg.grp)) == 1) stop('If background is not partitioned (non-spatial), do not assign ref.data to "bg".')
  
  if(any(is.null(occs.z), is.null(occs.grp))) {
    pts.plot <- bg.z |> dplyr::mutate(type = rep(0, nrow(bg.z)), partition = factor(bg.grp))  
  }else if(any(is.null(bg.z), is.null(bg.grp))) {
    pts.plot <- occs.z |> dplyr::mutate(type = rep(1, nrow(occs.z)), partition = factor(occs.grp))  
  }else{
    pts.plot <- rbind(occs.z, bg.z) |> as.data.frame() |>
      dplyr::mutate(type = c(rep(1, nrow(occs.z)), rep(0, nrow(bg.z))), partition = factor(c(occs.grp, bg.grp)))
  }
  names(pts.plot)[1:2] <- c("longitude","latitude")
  
  # find factor rasters or columns and identify them as categoricals
  if(!is.null(envs)) {
    categoricals <- unique(c(categoricals, names(envs)[which(terra::is.factor(envs))]))
    if(length(categoricals) == 0) categoricals <- NULL
  }else{
    categoricals <- unique(c(categoricals, names(occs.z)[which(sapply(occs.z, is.factor))]))
    if(length(categoricals) == 0) categoricals <- NULL
  }
  
  # if categoricals argument was specified, convert these columns to factor class
  if(!is.null(categoricals)) {
    for(i in 1:length(categoricals)) {
      if(!quiet) message(paste0("* Ignoring categorical variable ", categoricals[i], "..."))
      pts.plot[, categoricals[i]] <- NULL
    }
  }
  
  if(all(unique(pts.plot$partition) == 0)) {
    if(ref.data == "bg") stop('If using fully withheld testing data, input ref.data as "occs".')
    if(is.null(e) & is.null(occs.testing.z)) stop("If using fully withheld testing data, input either an ENMevaluation object or occs.testing.z.")
    if(!is.null(e)) occs.testing.z <- e@occs.testing
    names(occs.testing.z)[1:2] <- c("longitude","latitude")
    if(!is.null(categoricals)) occs.testing.z[[categoricals]] <- NULL
    occs.testing.z <- occs.testing.z |> dplyr::mutate(type = 1, partition = 2)
    pts.plot$partition <- as.numeric(as.character(pts.plot$partition))
    pts.plot[pts.plot$type == 1, "partition"] <- 1
    pts.plot <- dplyr::bind_rows(pts.plot, occs.testing.z) |> dplyr::mutate(partition = factor(partition))
  }
  
  return(pts.plot)
}

#' @title Similarity histogram plots for partition groups
#' @description Plots environmental similarity of reference partitions (occurrences or 
#' background) to remaining data (occurrence and background for all other partitions). This 
#' function does not use raster data, and thus only calculates similarity values for data used 
#' in model training. Further, this function does not calculate similarity for categorical 
#' variables.
#' @details When fully withheld testing groups are used, make sure to input either an 
#' ENMevaluation object or the argument occs.testing.z. In the resulting plot, partition 1 
#' refers to the training data, while partition 2 refers to the fully withheld testing group.
#' @param e ENMevaluation object
#' @param occs.z data frame: longitude, latitude, and environmental predictor variable values for occurrence records, in that order (optional);
#' the first two columns must be named "longitude" and "latitude"
#' @param occs.grp numeric vector: partition groups for occurrence records (optional)
#' @param bg.z data frame: longitude, latitude, and environmental predictor variable values for background records, in that order (optional);
#' the first two columns must be named "longitude" and "latitude"
#' @param bg.grp numeric vector: partition groups for background records (optional)
#' @param ref.data character: the reference to calculate MESS based on occurrences ("occs") or background ("bg"), with default "occs"
#' @param sim.type character: either "mess" for Multivariate Environmental Similarity Surface, "most_diff" for most different variable,
#' or "most_sim" for most similar variable; uses similarity function from package rmaxent
#' @param categoricals character vector: names of categorical variables in input SpatRaster or data frames to be removed from the analysis;
#' these must be specified as this function was intended for use with continuous data only; these must be specified when inputting tabular data instead of an ENMevaluation object 
#' @param envs.vars character vector: names of a predictor variable to plot similarities for; if left NULL, calculations are done
#' with respect to all variables (optional) 
#' @param occs.testing.z data frame: longitude, latitude, and environmental predictor variable values for fully withheld testing records, 
#' in that order; this is for use only with the "testing" partition option when an ENMevaluation object is not input (optional)
#' @param hist.bins numeric: number of histogram bins for histogram plots; default is 30
#' @param return.tbl boolean: if TRUE, return the data frames of similarity values used to make the ggplot instead of the plot itself
#' @param quiet boolean: if TRUE, silence all function messages (but not errors)
#' @details Histograms are plotted showing the environmental similarity estimates for each 
#' partition group. The similarity between environmental values associated with the 
#' validation occurrence or background records per partition group and those associated with 
#' the remaining data (occurrences and background) are calculated, and the minimum similarity 
#' per grid is returned. For option "mess", higher negative values indicate greater 
#' environmental difference between the validation occurrences and the study extent, and higher 
#' positive values indicate greater similarity. This function uses the `similarity()` function 
#' from the package `rmaxent` (https://github.com/johnbaums/rmaxent/), updated for `terra`, to calculate the 
#' similarities. Please see the below reference for details on MESS.
#' @return A ggplot of environmental similarities between the occurrence or background data 
#' for each partition and the rest of the data (all other occurrences and background data).
#' @references 
#' Baumgartner J, Wilson P (2021). _rmaxent: Tools for working with Maxent in R_. R package version 0.8.5.9000, <URL: https://github.com/johnbaums/rmaxent>.
#' Elith, J., Kearney, M., and Phillips, S. (2010) The art of modelling range-shifting species. \emph{Methods in Ecology and Evolution}, \bold{1}: 330-342. \doi{doi:10.1111/j.2041-210X.2010.00036.x}
#' 
#' @export

evalplot.envSim.hist <- function(e = NULL, occs.z = NULL, bg.z = NULL, occs.grp = NULL, 
                                 bg.grp = NULL, ref.data = "occs", sim.type = "mess", 
                                 categoricals = NULL, envs.vars = NULL, occs.testing.z = NULL,
                                 hist.bins = 30, return.tbl = FALSE, quiet = FALSE) {
  
  pts.plot <- plot.sim.dataPrep(e, envs = NULL, occs.z, bg.z, occs.grp, bg.grp, ref.data, categoricals, occs.testing.z, quiet)
  
  envs.names <- pts.plot |> dplyr::select(-longitude, -latitude, -partition, -type) |> names()
  
  if(!is.null(envs.vars)) {
    if(!quiet) message(paste0("* Similarity values calculated based only on ", paste(envs.vars, collapse = ", "), "."))
    envs.rem <- envs.names[-which(envs.names %in% envs.vars)]
    pts.plot <- pts.plot |> dplyr::select(-dplyr::all_of(envs.rem))
  }
  
  test.sim <- list()
  nk <- length(unique(pts.plot$partition[pts.plot$partition != 0]))
  if(nk == sum(pts.plot$type)) {
    stop("This plotting function is not available for jackknife (leave-one-out) partitions.")
  }
  
  for(k in 1:nk) {
    test.z <- pts.plot |> dplyr::filter(partition == k) |> dplyr::select(-longitude, -latitude, -partition)
    if(ref.data == "occs") {
      test.z <- test.z |> dplyr::filter(type == 1) |> dplyr::select(-type)
    }else if (ref.data == "bg") {
      test.z <- test.z |> dplyr::filter(type == 0) |> dplyr::select(-type)
    }
    train.z <- pts.plot |> dplyr::filter(partition != k) |> dplyr::select(-longitude, -latitude, -partition, -type)
    
    sim <- tryCatch({
      similarity(train.z, test.z)
    }, error = function(cond) {
      message('Error: there may be a categorical variable present in the predictor variable data. Please make sure to declare all categorical variables with the "categoricals" argument.')
      # Choose a return value in case of error
      return(NULL)
    })
    sim.sel <- switch(sim.type, mess = sim$similarity_min, most_diff = sim$mod, most_sim = sim$mos)  
    
    test.sim[[k]] <- data.frame(partition = k, sim.sel)
    names(test.sim[[k]])[2] <- sim.type  
  }
  
  if(nk > 9) {
    theme.custom <- ggplot2::guides(color = FALSE)
    pt.cols <- rainbow(nk)
  }else{
    theme.custom <- NULL
    pt.cols <- RColorBrewer::brewer.pal(nk, "Set1")
  }
  
  plot.df <- dplyr::bind_rows(test.sim)
  plot.df$partition <- factor(plot.df$partition)
  
  plot.text <- paste("\n(Values represent environmental similarity between", 
                     switch(ref.data, occs = "occurrence", bg = "background"),
                     "partitions and all other background partitions.)")
  
  if(sim.type != "mess") {
    facts <- sort(unique(plot.df[,2]))
    envs.tbl <- data.frame(facts, envs.names[facts])
    names(envs.tbl) <- c(sim.type, "env.var")
    envs.tbl$env.var <- factor(envs.tbl$env.var, levels = envs.tbl$env.var)
    plot.df <- plot.df |> dplyr::left_join(envs.tbl, by = sim.type)
    plot.df[[sim.type]] <- NULL
    names(plot.df)[2] <- sim.type
    title <- paste(switch(sim.type, most_diff = "Most different", most_sim = "Most similar"), "environmental variable")
    
    g <- ggplot2::ggplot(plot.df, ggplot2::aes_string(x = sim.type, fill = "partition")) + 
      ggplot2::stat_count() +
      ggplot2::facet_grid(ggplot2::vars(partition)) + 
      ggplot2::scale_fill_manual(values = pt.cols) +
      ggplot2::theme_classic() +
      ggplot2::geom_vline(xintercept = 0) +
      ggplot2::ggtitle(paste(title, plot.text, collapse = "\n")) +
      ggplot2::theme(strip.background = ggplot2::element_blank(), strip.text.y = ggplot2::element_blank())
  }else if(sim.type == "mess"){
    g <- ggplot2::ggplot(plot.df, ggplot2::aes(x = mess, fill = partition)) + 
      ggplot2::geom_histogram(bins = hist.bins) +
      ggplot2::facet_grid(ggplot2::vars(partition)) + 
      ggplot2::scale_fill_manual(values = pt.cols) +
      ggplot2::theme_classic() +
      ggplot2::geom_vline(xintercept = 0) +
      ggplot2::ggtitle(paste("Multivariate environmental similarity", plot.text, collapse = "\n")) +
      ggplot2::theme(strip.background = ggplot2::element_blank(), strip.text.y = ggplot2::element_blank())  
  }
  g
  
  if(return.tbl == TRUE) {
    return(plot.df)
  }else{
    return(g)  
  }
}

#' @title Similarity maps for partition groups
#' @description Maps environmental similarity of reference partitions (occurrences or 
#' background) to all cells with values in the predictor variable rasters. This function uses 
#' raster data, and thus cannot map similarity values using only tables of environmental values f
#' or occurrences or background. Further, this function does not calculate similarity for 
#' categorical variables.
#' @details When fully withheld testing groups are used, make sure to input either an ENMevaluation 
#' object or the argument occs.testing.z. In the resulting plot, partition 1 refers to the training data,
#' while partition 2 refers to the fully withheld testing group.
#' @param e ENMevaluation object (optional) 
#' @param envs SpatRaster: environmental predictor variables used to build the models in "e"; categorical variables should be 
#' removed before input or identified with the argument "categoricals", as they cannot be used to calculate MESS
#' @param occs.z data frame: longitude, latitude, and environmental predictor variable values for occurrence records, in that order (optional);
#' the first two columns must be named "longitude" and "latitude"
#' @param occs.grp numeric vector: partition groups for occurrence records (optional)
#' @param bg.z data frame: longitude, latitude, and environmental predictor variable values for background records, in that order (optional);
#' the first two columns must be named "longitude" and "latitude"
#' @param bg.grp numeric vector: partition groups for background records (optional)
#' @param ref.data character: the reference to calculate MESS based on occurrences ("occs") or background ("bg"), with default "occs"
#' @param sim.type character: either "mess" for Multivariate Environmental Similarity Surface, "most_diff" for most different variable,
#' or "most_sim" for most similar variable; uses similarity function from package rmaxent
#' @param categoricals character vector: names of categorical variables in input SpatRaster or data frames to be removed from the analysis;
#' these must be specified as this function was intended for use with continuous data only
#' @param envs.vars character vector: names of a predictor variable to plot similarities for; if left NULL, calculations are done
#' with respect to all variables (optional) 
#' @param bb.buf numeric: distance used to buffer (extend) the mapping extent in map units; for latitude/longitude, this is in degrees (optional)
#' @param occs.testing.z data frame: longitude, latitude, and environmental predictor variable values for fully withheld testing records, 
#' in that order; this is for use only with the "testing" partition option when an ENMevaluation object is not input (optional)
#' @param plot.bg.pts boolean: if TRUE, plot background points when using ref.data = "bg"
#' @param sim.palette character: RColorBrewer palette name to use for plotting discrete variables; if NULL, default is "Set1"
#' @param pts.size numeric: custom point size for ggplot
#' @param gradient.colors character vector: colors used for ggplot2::scale_fill_gradient2
#' @param na.color character: color used for NA values
#' @param return.tbl boolean: if TRUE, return the data frames of similarity values used to make the ggplot instead of the plot itself
#' @param return.ras boolean: if TRUE, return the SpatRaster of similarity values used to make the ggplot instead of the plot itself
#' @param quiet boolean: if TRUE, silence all function messages (but not errors)
#' @details Rasters are plotted showing the environmental similarity estimates for each 
#' partition group. The similarity between environmental values associated with the 
#' validation occurrence or background records per partition group and those associated with 
#' the entire study extent (specified by the extent of the input SpatRaster "envs") are 
#' calculated, and the minimum similarity per grid is returned. For option "mess", higher 
#' negative values indicate greater environmental difference between the validation occurrences 
#' and the study extent, and higher positive values indicate greater similarity. This function 
#' uses the `similarity()` function from the package `rmaxent` 
#' (https://github.com/johnbaums/rmaxent/) to calculate the similarities. Please see the below 
#' reference for details on MESS. 
#' @return A ggplot of environmental similarities between the occurrence or background data 
#' for each partition and all predictor variable values in the extent.
#' @references 
#' Baumgartner J, Wilson P (2021). _rmaxent: Tools for working with Maxent in R_. R package version 0.8.5.9000, <URL: https://github.com/johnbaums/rmaxent>.
#' Elith, J., Kearney, M., and Phillips, S. (2010) The art of modelling range-shifting species. \emph{Methods in Ecology and Evolution}, \bold{1}: 330-342. \doi{doi:10.1111/j.2041-210X.2010.00036.x}
#' 
#' @export

evalplot.envSim.map <- function(e = NULL, envs, occs.z = NULL, bg.z = NULL, occs.grp = NULL, 
                                bg.grp = NULL, ref.data = "occs", sim.type = "mess", 
                                categoricals = NULL, envs.vars = NULL, bb.buf = 0, occs.testing.z = NULL,
                                plot.bg.pts = FALSE, sim.palette = NULL, 
                                pts.size = 1.5, gradient.colors = c("red","white","blue"), na.color = "gray",
                                return.tbl = FALSE, return.ras = FALSE, quiet = FALSE) {
  
  if(return.tbl == TRUE & return.ras == TRUE) {
    stop("*Error: please select only one of return.tbl or return.ras.")
  }
  
  if(is.null(e) & (ref.data == "occs" & any(is.null(occs.z), is.null(occs.grp)))) {
    stop("* Error: If using occurrences as reference group, ensure you input occs.z and occs.grp") 
  }
  if(is.null(e) & (ref.data == "bg" & any(is.null(bg.z), is.null(bg.grp)))) {
    stop("* Error: If using background as reference group, ensure you input bg.z and bg.grp") 
  }
  
  if(!is.numeric(bb.buf)) stop("Please ensure bb.buf is a number.")
  
  pts.plot <- plot.sim.dataPrep(e, envs, occs.z, bg.z, occs.grp, bg.grp, ref.data, categoricals, occs.testing.z, quiet)
  
  if(!is.null(categoricals) & !is.null(envs)) {
    envs <- terra::subset(envs, -which(names(envs) %in% categoricals))
  }
  
  if(!is.null(envs.vars)) {
    if(!quiet) message(paste0("* Similarity values calculated based only on ", 
                              paste(envs.vars, collapse = ", "), "."))
    envs.names <- names(envs)
    envs.rem <- envs.names[-which(envs.names %in% envs.vars)]
    pts.plot <- pts.plot |> dplyr::select(-dplyr::all_of(envs.rem))
    envs <- envs[[envs.vars]]
  }
  
  if(ref.data == "occs") {
    pts.plot <- pts.plot |> dplyr::filter(type == 1) |> dplyr::select(-type)
  }else if (ref.data == "bg") {
    pts.plot <- pts.plot |> dplyr::filter(type == 0) |> dplyr::select(-type)
  }
  
  ras.sim <- list()
  nk <- length(unique(pts.plot$partition))
  if(nk == sum(pts.plot$type)) {
    stop("This plotting function is not available for jackknife (leave-one-out) partitions.")
  }
  
  for(k in 1:nk) {
    test.z <- pts.plot |> dplyr::filter(partition == k) |> 
      dplyr::select(-longitude, -latitude, -partition)
    
    sim <- tryCatch({
      similarity(envs, test.z)
    }, error = function(cond) {
      message('Error: there may be a categorical variable present in the predictor variable data. Please make sure to declare all categorical variables with the "categoricals" argument.')
      # Choose a return value in case of error
      return(NULL)
    })
    sim.sel <- switch(sim.type, mess = sim$similarity_min, 
                      most_diff = sim$mod, most_sim = sim$mos)  
    names(sim.sel) <- paste0("partition", k)
    
    ras.sim[[k]] <- sim.sel
  }
  
  rs.sim <- terra::rast(ras.sim)
  plot.df <- terra::as.data.frame(rs.sim, xy = TRUE, na.rm = FALSE) |>
    tidyr::pivot_longer(cols = 3:dplyr::last_col(), names_to = "ras", values_to = sim.type)
  # add buffer
  plot.df <- plot.df |> dplyr::filter(x > min(pts.plot$longitude) - bb.buf, 
                                       x < max(pts.plot$longitude) + bb.buf,
                                       y > min(pts.plot$latitude) - bb.buf, 
                                       y < max(pts.plot$latitude) + bb.buf)
  
  if(sim.type != "mess") {
    if(is.null(sim.palette)) sim.palette <- "Set1"
    plot.df$ras <- gsub("_var", "", plot.df$ras)
    plot.df[[sim.type]] <- factor(plot.df[[sim.type]])
    title <- paste(switch(sim.type, most_diff = "Most different", 
                          most_sim = "Most similar"), "environmental variable")
  }else{
    title <- "Multivariate environmental similarity" 
  }
  
  plot.text <- paste("\n(Values represent environmental similarity between", 
                     switch(ref.data, occs = "occurrence", bg = "background"),
                     "partitions and all raster cells with values.)")
  
  g <- ggplot2::ggplot() +
    ggplot2::geom_raster(data = plot.df, ggplot2::aes_string(x = "x", y = "y", fill = sim.type))
  if(ref.data == "bg" & plot.bg.pts == FALSE) {
    g <- g + ggplot2::facet_wrap(ggplot2::vars(ras), ncol = 2) +
      ggplot2::ggtitle(paste(title, plot.text, collapse = "\n")) +
      ggplot2::theme_classic()
  }else{
    g <- g + ggplot2::geom_point(data = pts.plot, ggplot2::aes(x = longitude, y = latitude, shape = partition), color = "black", size = pts.size) +
      ggplot2::facet_wrap(ggplot2::vars(ras), ncol = 2) +
      ggplot2::ggtitle(paste(title, plot.text, collapse = "\n")) +
      ggplot2::theme_classic()
  }
  if(sim.type != "mess") {
    g <- g + ggplot2::scale_fill_brewer(palette = sim.palette, na.value = "white", breaks = levels(plot.df[[sim.type]]))
  }else{
    g <- g + ggplot2::scale_fill_gradient2(low = gradient.colors[1], mid = gradient.colors[2], high = gradient.colors[3], na.value = na.color)
  }
  g
  
  if(return.tbl == TRUE) {
    return(plot.df)
  }else if(return.ras == TRUE) {
    return(rs.sim)
  }else{
    return(g)  
  }
}


#' @title ENMevaluation statistics plot
#' @description Plot evaluation statistics over tuning parameter ranges to visualize differences in performance
#' @param e ENMevaluation object
#' @param stats character vector: names of statistics from results table to be plotted; if more than
#' one statistic is specified, the plot will be faceted
#' @param x.var character: variable to be plotted on x-axis
#' @param color.var character: variable used to assign symbology colors 
#' @param dodge numeric: dodge width for points and lines; this improves visibility when there is high overlap (optional)
#' @param error.bars boolean: if TRUE, plot error bars
#' @param facet.labels character vector: custom names for the metric facets
#' @param metric.levels character vector: custom factor levels for metrics; this controls the order that metric statistics are plotted
#' @param return.tbl boolean: if TRUE, return the data frames of results used to make the ggplot instead of the plot itself
#' @details In this plot, the x-axis represents a tuning parameter range, while the y-axis represents the average of a statistic over all partitions.
#' Different colors represent the categories or values of another tuning parameter. 
#' Error bars represent the standard deviation of a statistic around the mean. 
#' Currently, this function can only plot two tuning parameters at a time.
#' @return A ggplot of evaluation statistics. 
#' @export

evalplot.stats <- function(e, stats, x.var, color.var, dodge = NULL, error.bars = TRUE, facet.labels = NULL, metric.levels = NULL, return.tbl = FALSE) {
  exp <- paste(paste0("*", stats), collapse = "|")
  res <- e@results |> 
    tidyr::pivot_longer(cols = auc.train:ncoef, names_to = "metric", values_to = "value") |>
    dplyr::filter(grepl(exp, metric))
  avgs <- res |> 
    dplyr::filter(grepl("avg", metric)) |>
    dplyr::rename(avg = value) |>
    dplyr::mutate(metric = gsub(".avg", "", metric))
  sds <- res |> 
    dplyr::filter(grepl("sd", metric)) |>
    dplyr::rename(sd = value) |>
    dplyr::mutate(metric = gsub(".sd", "", metric))
  join.names <- names(avgs)
  join.names <- join.names[join.names != "avg"]
  res.avgs <- dplyr::left_join(avgs, sds, by = join.names) |>
    dplyr::mutate(lower = avg - sd, upper = avg + sd,
                  metric = factor(metric, levels = stats))
  if(!is.null(facet.labels)) labeller <- ggplot2::as_labeller(facet.labels) else labeller <- NULL
  if(!is.null(metric.levels)) res$metric <- factor(res$metric, levels = metric.levels)
  if(!is.null(metric.levels)) res.avgs$metric <- factor(res.avgs$metric, levels = metric.levels)
  
  if(nrow(res.avgs) > 0) {
    if(is.null(dodge)) dodge <- 0.1
    g <- ggplot2::ggplot(res.avgs, ggplot2::aes_string(x = x.var, y = "avg", color = color.var, group = color.var)) + 
      ggplot2::geom_point(position=ggplot2::position_dodge(width=dodge)) + 
      ggplot2::geom_line(position=ggplot2::position_dodge(width=dodge)) +
      ggplot2::theme_bw()
    if(length(stats) > 1) {
      if(!is.null(labeller)) {
        g <- g + ggplot2::facet_wrap(ggplot2::vars(metric), scales = "free_y", nrow = length(stats), labeller = labeller)
      }else{
        g <- g + ggplot2::facet_wrap(ggplot2::vars(metric), scales = "free_y", nrow = length(stats))  
      }
    }
    if(error.bars == TRUE) {
      g <- g + ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper), width = 0.5, 
                                      position = ggplot2::position_dodge(width=dodge))
    }
    g
    if(return.tbl == TRUE) {
      return(res.avgs)
    }else{
      return(g)  
    }
  }else{
    if(is.null(dodge)) dodge <- 0
    g <- ggplot2::ggplot(res, ggplot2::aes_string(x = x.var, y = "value", color = color.var, group = color.var)) + 
      ggplot2::geom_point(position=ggplot2::position_dodge(width=dodge)) + 
      ggplot2::geom_line(position=ggplot2::position_dodge(width=dodge)) + 
      ggplot2::theme_bw()
    if(!is.null(labeller)) {
      g <- g + ggplot2::facet_wrap(ggplot2::vars(metric), scales = "free_y", nrow = length(stats), labeller = labeller)
    }else{
      g <- g + ggplot2::facet_wrap(ggplot2::vars(metric), scales = "free_y", nrow = length(stats))
    }
    g
    if(return.tbl == TRUE) {
      return(res)
    }else{
      return(g)  
    }
  }
}

#' @title ENMnulls statistics plot
#' @description Plot histogram of evaluation statistics for null ENM simulations
#' @param e.null ENMnull object
#' @param stats character vector: metrics from results table to be plotted; examples are
#' "auc.val" or "or.10p"; if more than one statistic is specified, the histogram plot will be faceted
#' @param plot.type character: either "violin" or "histogram"
#' @param facet.labels named list: custom names for the metric facets, in the form list(old_name = "new_name", ...)
#' @param metric.levels character vector: custom factor levels for metrics; this controls the order that metric statistics are plotted
#' @param return.tbl boolean: if TRUE, return the data frames of null results used to make the ggplot instead of the plot itself
#' @details There are two variations for this plot, but both show null quantiles (0.01, 0.05, 0.5, 0.95, 0.99). 
#' For violin plots, the null distribution is displayed as a vertical shape (i.e., the violin) with horizontal lines showing 
#' the quantiles and the empirical value is plotted as a red point along the vertical axis. 
#' For histogram plots, the null distribution is displayed as a histogram with vertical lines showing the quantiles 
#' and the empirical value as a vertical red line on the distribution.
#' @return A ggplot of null model statistics. 
#' @export

evalplot.nulls <- function(e.null, stats, plot.type, facet.labels = NULL, metric.levels = NULL, return.tbl = FALSE) {
  exp <- paste(paste0("*", stats), collapse = "|")
  null.res <- e.null@null.results |> 
    tidyr::pivot_longer(cols = auc.train:ncoef, names_to = "metric", values_to = "value") |>
    dplyr::filter(grepl(exp, metric)) |>
    dplyr::select(metric, value)
  null.avgs <- null.res |> 
    dplyr::filter(grepl("avg", metric) | metric %in% stats) |>
    dplyr::rename(avg = value) |>
    dplyr::mutate(metric = gsub(".avg", "", metric))
  if(!is.null(metric.levels)) null.avgs$metric <- factor(null.avgs$metric, levels = metric.levels)
  # null.sds <- null.res |> 
  #   dplyr::filter(grepl("sd", metric)) |>
  #   dplyr::rename(sd = value) |>
  #   dplyr::mutate(metric = gsub(".sd", "", metric))
  # null.res.avgs <- dplyr::bind_cols(null.avgs, null.sds |> dplyr::select(sd))
  
  emp.res <- e.null@null.emp.results |> 
    dplyr::slice(1) |>
    tidyr::pivot_longer(cols = stats, names_to = "metric", values_to = "value") |>
    dplyr::select(statistic, metric, value) |>
    tidyr::pivot_wider(names_from = statistic, values_from = value) |>
    dplyr::rename(avg = emp.mean)
  
  if(!is.null(facet.labels)) labeller <- ggplot2::as_labeller(facet.labels) else labeller <- NULL
  
  if(plot.type == "violin") {
    g <- ggplot2::ggplot(null.avgs, ggplot2::aes(x = metric, y = avg)) + 
      ggplot2::geom_violin(draw_quantiles = c(0.01, 0.05, 0.5, 0.95, 0.99)) +
      ggplot2::geom_point(data = emp.res, ggplot2::aes(y = avg), color = "red") +
      ggplot2::theme_bw()  
  }else if(plot.type == "histogram") {
    stats.all <- rbind(null.avgs, emp.res)
    vlines <- null.avgs |> dplyr::group_by(metric) |> 
      dplyr::summarize(`0.01 quantile` = quantile(avg, 0.01),
                       `0.05 quantile` = quantile(avg, 0.05),
                       `0.50 quantile` = quantile(avg, 0.5),
                       `0.95 quantile` = quantile(avg, 0.95),
                       `0.99 quantile` = quantile(avg, 0.99)) |>
      tidyr::pivot_longer(cols = `0.01 quantile`:`0.99 quantile`, names_to = "quantile", values_to = "value")
    vlines <- rbind(vlines, emp.res |> dplyr::mutate(quantile = "empirical value") |> dplyr::rename(value = avg))
    g <- ggplot2::ggplot(mapping = ggplot2::aes(x = avg)) + 
      ggplot2::geom_histogram(data = null.avgs, fill = "gray80") +
      ggplot2::geom_vline(data = vlines, ggplot2::aes(xintercept = value, color = quantile, linetype = quantile, size = quantile)) +
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
      ggplot2::scale_size_manual(values = c(`0.01 quantile` = 1, 
                                            `0.05 quantile` = 1,
                                            `0.50 quantile` = 1,
                                            `0.95 quantile` = 1, 
                                            `0.99 quantile` = 1,
                                            `empirical value` = 0.5)) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.title=ggplot2::element_blank(), 
                     legend.position="bottom")
    if(!is.null(labeller)) {
      g <- g + ggplot2::facet_wrap(ggplot2::vars(metric), scales = "free_x", ncol = 1, labeller = labeller)
    }else{
      g <- g + ggplot2::facet_wrap(ggplot2::vars(metric), scales = "free_x", ncol = 1)
    }
  }
  if(return.tbl == TRUE) {
    return(list(null.avgs = null.avgs, empirical.results = emp.res))
  }else{
    return(g)  
  }
}
