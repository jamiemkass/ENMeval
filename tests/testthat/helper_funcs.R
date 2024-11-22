# testing functions 
test_ENMevaluation <- function(e, alg, parts, tune.args, nparts.occs, nparts.bg, type = "") {
  tune.tbl <- expand.grid(tune.args, stringsAsFactors = FALSE)
  test_that("ENMevaluation object and slots exist", {
    expect_true(!is.null(e))
    expect_true(!is.null(e@algorithm))
    expect_true(!is.null(e@tune.settings))
    expect_true(!is.null(e@partition.method))
    expect_true(!is.null(e@results))
    expect_true(!is.null(e@results.partitions))
    expect_true(!is.null(e@models))
    expect_true(!is.null(e@predictions))
    expect_true(terra::nlyr(e@predictions) > 0)  
    expect_true(!is.null(e@occs))
    expect_true(!is.null(e@occs.grp))
    expect_true(!is.null(e@bg))
    expect_true(!is.null(e@bg.grp))
    expect_true(!is.null(e@overlap))
    expect_equal(length(slotNames(e)), 20)
    expect_equal(slotNames(e),
                 c("algorithm", "tune.settings", "partition.method",
                   "partition.settings", "other.settings", "doClamp",
                   "clamp.directions", "results", "results.partitions",
                   "models", "variable.importance", "predictions", "taxon.name",
                   "occs", "occs.testing", "occs.grp", "bg", "bg.grp",
                   "overlap", "rmm"))
  })  
  
  test_that("Data in ENMevaluation object slots have correct form", {
    # algorithm
    expect_true(e@algorithm == alg)
    # partition method 
    expect_true(e@partition.method == parts)
    # these checks relate to tune.args, which may be NULL
    if(!is.null(tune.args)) {
      # tune.settings 
      expect_true(all(as.data.frame(e@tune.settings[,1:ncol(tune.tbl)]) == tune.tbl))
      # nrow of results
      expect_true(nrow(e@results) == nrow(tune.tbl))
      # tune.args column values are concat of tuning parameters columns
      # expect_true(all(apply(e@results[names(tune.args.ls[[m]])], 1, paste, collapse = "_") == as.character(e@results$tune.args)))
      # number of models
      expect_true(length(e@models) == nrow(tune.tbl))
    }
    # number of rows for occs matches occs.grp
    expect_true(nrow(e@occs) == length(e@occs.grp))
    # number of rows for bg matches bg.grp
    expect_true(nrow(e@bg) == length(e@bg.grp))
    # no overlap is calculated for no tuning or BIOCLIM
    if(length(e@overlap) > 0) {
      # both indices exist for overlap
      expect_true(length(e@overlap) == 2)
      # number of rows of overlap D matches tune.args
      expect_true(nrow(e@overlap$D) == nrow(tune.tbl))
      # number of rows of overlap I matches tune.args
      expect_true(nrow(e@overlap$I) == nrow(tune.tbl))  
    }else{
      # no overlap matrix
      expect_true(length(e@overlap) == 0)
    }
  })
  
  test_that("Records with missing environmental values were removed", {
    expect_true(sum(is.na(e@occs)) == 0)
    expect_true(sum(is.na(e@bg)) == 0)
  })
  
  test_that("Number of partitions is correct", {
    expect_true(length(unique(e@occs.grp)) == nparts.occs)
    expect_true(length(unique(e@bg.grp)) == nparts.bg)
  })
  
  
  test_that("Results table for partitions has correct form", {
    if(parts == "none") {
      expect_true(nrow(e@results.partitions) == 0)
    }else{
      expect_true(nrow(e@results.partitions) == nparts.occs * nrow(tune.tbl))
      if(parts != "testing") {
        expect_true(max(e@results.partitions$fold) == nparts.occs)
      }else{
        expect_true(max(e@results.partitions$fold) == 0)
      }
      # jackknife has NAs for cbi.val
      if(parts == "jackknife" | !requireNamespace("ecospat", quietly = TRUE)) {
        expect_true(sum(is.na(e@results.partitions)) == nrow(e@results.partitions))
      }else{
        expect_true(sum(is.na(e@results.partitions)) == 0)
      }
    }
  })
}

test_clamp <- function(envs, occs.z, bg.z, categoricals) {
  # use occurrences as reference environmental values for clamping
  # restrict to small subset in Amazon to ensure lots of extrapolation for 
  # transfers
  p.z <- occs.z[8:13,-1:-2]
  # remove cats from dataset
  if(!is.null(categoricals)) {
    envs <- envs[[-which(names(envs) %in% categoricals)]]
    p.z <- p.z |> dplyr::select(-categoricals)
  }
  # for sub clamping
  z <- 2:4
  # p.z <- dplyr::bind_rows(occs.z, bg.z)[,-1:-2]
  
  none <- envs
  all <- clamp.vars(orig.vals = envs, ref.vals = p.z)
  left <- clamp.vars(orig.vals = envs, ref.vals = p.z, right = "none",)
  right <- clamp.vars(orig.vals = envs, ref.vals = p.z, left = "none",)
  subboth <- clamp.vars(orig.vals = envs, ref.vals = p.z, left = names(envs)[z], 
                        right = names(envs)[z])
  subleft <- clamp.vars(orig.vals = envs, ref.vals = p.z, right = "none", 
                        left = names(envs)[z])
  subright <- clamp.vars(orig.vals = envs, ref.vals = p.z, left = "none", 
                         right = names(envs)[z])
  # clamps.envs <- list(none=none, all=all, left=left, right=right, subboth=subboth, subleft=subleft, subright=subright)
  
  # enm <- lookup.enm(e@algorithm)
  # m <- e@models[[1]]
  
  # clamp.envs.p <- lapply(clamps.envs, function(x) enm@predict(m, x, list(doClamp = FALSE, pred.type = "cloglog")))
  # combs <- expand.grid(x=names(clamp.envs.p), y=names(clamp.envs.p), stringsAsFactors = FALSE) |> dplyr::filter(x != y)
  
  test_that("Clamped rasters are different from each other", {
    expect_true(all(terra::minmax(none) != terra::minmax(all)))
    expect_true(all(terra::minmax(none)[1,] != terra::minmax(left)[1,]))
    expect_true(all(terra::minmax(none)[2,] == terra::minmax(left)[2,]))
    expect_true(all(terra::minmax(none)[1,] == terra::minmax(right)[1,]))
    expect_true(all(terra::minmax(none)[2,] != terra::minmax(right)[2,]))
    
    expect_true(all(terra::minmax(none)[,z] != terra::minmax(subboth)[,z]))
    expect_true(all(terra::minmax(none)[,-z] == terra::minmax(subboth)[,-z]))
    expect_true(all(terra::minmax(none)[1,z] != terra::minmax(subleft)[1,z]))
    expect_true(all(terra::minmax(none)[2,z] == terra::minmax(subleft)[2,z]))
    expect_true(all(terra::minmax(none)[1,-z] == terra::minmax(subleft)[1,-z]))
    expect_true(all(terra::minmax(none)[2,-z] == terra::minmax(subleft)[2,-z]))
    expect_true(all(terra::minmax(none)[1,z] == terra::minmax(subright)[1,z]))
    expect_true(all(terra::minmax(none)[2,z] != terra::minmax(subright)[2,z]))
    expect_true(all(terra::minmax(none)[1,-z] == terra::minmax(subright)[1,-z]))
    expect_true(all(terra::minmax(none)[2,-z] == terra::minmax(subright)[2,-z]))
    # for(i in 1:nrow(combs)) {
    # if(canExtrapolate == TRUE) {
    # expect_false(all(abs(terra::minmax(clamp.envs.p[[combs[i,1]]] - clamp.envs.p[[combs[i,2]]])) < 1e-7))
    # }else{
    # expect_true(all(abs(terra::minmax(clamp.envs.p[[combs[i,1]]] - clamp.envs.p[[combs[i,2]]])) < 1e-7))
    # }
    # }
  })
  
  # clamps.df <- lapply(clamps.envs, function(x) terra::values(x))
  # clamp.df.p <- lapply(clamps.df, function(x) enm@predict(m, x, list(doClamp = FALSE, pred.type = "cloglog")))
  
  # test_that("Clamped data frames are different from each other", {
  #   for(i in 1:nrow(combs)) {
  #     if(canExtrapolate == TRUE) {
  #       expect_false(isTRUE(all.equal(clamp.df.p[[combs[i,1]]], clamp.df.p[[combs[i,2]]])))
  #     }else{
  #       expect_true(isTRUE(all.equal(clamp.df.p[[combs[i,1]]], clamp.df.p[[combs[i,2]]])))
  #     }
  #   }
  # })
}

test_ENMnulls <- function(e, ns, no.iter, alg, parts, mod.settings, nparts.occs, nparts.bg, n.sims, type = "") {
  mod.settings.tbl <- expand.grid(mod.settings)
  test_that("ENMnulls object and slots exist", {
    expect_true(!is.null(ns))
    expect_true(!is.null(ns@null.algorithm))
    expect_true(!is.null(ns@null.mod.settings))
    expect_true(!is.null(ns@null.partition.method))
    expect_true(!is.null(ns@null.partition.settings))
    expect_true(!is.null(ns@null.other.settings))
    expect_true(!is.null(ns@null.no.iter))
    expect_true(!is.null(ns@null.results))
    expect_true(!is.null(ns@null.results.partitions))
    expect_true(!is.null(ns@null.emp.results))
    expect_true(!is.null(ns@emp.occs))
    expect_true(!is.null(ns@emp.occs.grp))
    expect_true(!is.null(ns@emp.bg))
    expect_true(!is.null(ns@emp.bg.grp))
  })  
  
  test_that("Data in ENMnulls object slots have correct form", {
    # algorithm
    expect_true(ns@null.algorithm == alg)
    # partition method 
    expect_true(ns@null.partition.method == parts)
    # mod.settings 
    if(ncol(mod.settings.tbl) > 1) {
      expect_true(all(ns@null.mod.settings[,1:ncol(mod.settings.tbl)] == mod.settings.tbl))  
    }else{
      expect_true(as.character(ns@null.mod.settings[,1]) == mod.settings.tbl)  
    }
    # no. of iterations
    expect_true(ns@null.no.iter == no.iter)
    # number of rows in results table
    expect_true(nrow(ns@null.results) == no.iter)
    # number of rows in results table for partitions
    if(ns@null.partition.method == "none") {
      expect_true(nrow(ns@null.results.partitions) == 0)
    }else{
      expect_true(nrow(ns@null.results.partitions) == no.iter * nparts.occs)  
    }
    
    # number of rows in empirical vs null results table
    expect_true(nrow(ns@null.emp.results) == 6)
    # check that not all empirical results are NA
    expect_true(sum(apply(ns@null.emp.results[,2:ncol(ns@null.emp.results)], 2, function(x) sum(is.na(x)))) != 
                  nrow(ns@null.emp.results) * (ncol(ns@null.emp.results)-1))
    # check that tables match
    expect_true(all(ns@emp.occs == e@occs))
    expect_true(all(ns@emp.bg == e@bg))
    expect_true(all(ns@emp.occs.grp == e@occs.grp))
  })
  
  test_that("Data in ENMnulls object slots are not NA (except CBI, which can be NA due to low data)", {
    expect_true(all(apply(ns@null.results |> dplyr::select(!starts_with("cbi")), 2, function(x) sum(is.na(x))) == 0))
    if(ns@null.partition.method != "none") expect_true(all(apply(ns@null.results |> dplyr::select(!starts_with("cbi")), 2, function(x) sum(is.na(x))) == 0))
  })
}

#' @title Unit tests for ENMevaluation plotting functions
#' @description All parameters are self-explanatory except the following.
#' Anything with the prefix ".z" is a data frame with longitude, latitude, 
#' and the environmental predictor variable values. Argument "plot.sel" controls
#' whether testing happens for the histogram function, the plotting function, or both
#' (some implementations do not work with one or the other). Argument "bg.sel" controls
#' whether tests should be done with ref.data as "bg" or not (non-spatial implementations
#' cannot be plotted with ref.data as )

test_evalplot.stats <- function(e) {
  
  test_stats <- function(x, stats) {
    if(e@partition.method == "testing") {
      y <- 2
      z <- c("metric", "value")
    }else{
      y <- 5
      z <- c("metric", "avg", "sd", "lower", "upper")
    }
    test_that("Outputs for evalplot.stats have correct form", {
      expect_true(ncol(x) == ncol(e@tune.settings) + y)
      expect_true(nrow(x) == nrow(e@tune.settings) * length(stats))
      expect_true(all(unique(x$metric) == stats))
      n <- (ncol(e@tune.settings)+1):(ncol(x))
      expect_true(all(names(x)[n] == z))
    })
  }
  
  if(e@partition.method == "none") {
    stat1 <- "auc.train"
    stat2 <- c("auc.train", "cbi.train")
  }else{
    stat1 <- "auc.val"
    stat2 <- c("auc.val", "or.10p") 
  }
  
  # defaults
  evalplot.stats(e, stats = stat1, x.var = "rm", color.var = "fc", dodge = FALSE, error.bars = FALSE, facet.labels = NULL, metric.levels = NULL, return.tbl = TRUE) |> test_stats(stat1)
  # two stats
  evalplot.stats(e, stats = stat2, x.var = "rm", color.var = "fc", dodge = FALSE, error.bars = FALSE, facet.labels = NULL, metric.levels = NULL, return.tbl = TRUE) |> test_stats(stat2)
  # dodge
  evalplot.stats(e, stats = stat1, x.var = "rm", color.var = "fc", dodge = TRUE, error.bars = FALSE, facet.labels = NULL, metric.levels = NULL, return.tbl = TRUE) |> test_stats(stat1)
  # error bars
  evalplot.stats(e, stats = stat1, x.var = "rm", color.var = "fc", dodge = FALSE, error.bars = TRUE, facet.labels = NULL, metric.levels = NULL, return.tbl = TRUE) |> test_stats(stat1)
  # facet labels
  evalplot.stats(e, stats = stat1, x.var = "rm", color.var = "fc", dodge = FALSE, error.bars = TRUE, facet.labels = paste0(stat1), metric.levels = NULL, return.tbl = TRUE) |> test_stats(stat1)
  evalplot.stats(e, stats = stat2, x.var = "rm", color.var = "fc", dodge = FALSE, error.bars = TRUE, facet.labels = paste0(stat2), metric.levels = NULL, return.tbl = TRUE) |> test_stats(stat2)
  # metric levels
  evalplot.stats(e, stats = stat1, x.var = "rm", color.var = "fc", dodge = FALSE, error.bars = TRUE, facet.labels = NULL, metric.levels = stat1, return.tbl = TRUE) |> test_stats(stat1)
  evalplot.stats(e, stats = stat2, x.var = "rm", color.var = "fc", dodge = FALSE, error.bars = TRUE, facet.labels = NULL, metric.levels = rev(stat2), return.tbl = TRUE) |> test_stats(stat2)
}

test_evalplot.envSim.hist <- function(e, occs.z, bg.z, occs.grp, bg.grp, bg.sel = 1, occs.testing.z = NULL) {
    test_hist <- function(i) {
      test_that("Outputs for evalplot.envSim.hist have correct form", {
        expect_true(ncol(i) == 2)
        expect_true(names(i)[1] == "partition")
      })
    }
    # with ENMevaluation object
    evalplot.envSim.hist(e = e, ref.data = "occs", return.tbl = TRUE, quiet = TRUE, occs.testing.z = occs.testing.z) |> test_hist()
    if(bg.sel == 1) {
      evalplot.envSim.hist(e = e, ref.data = "bg", return.tbl = TRUE, quiet = TRUE, occs.testing.z = occs.testing.z) |> test_hist()
    }
    evalplot.envSim.hist(e = e, ref.data = "occs", envs.vars = c("bio1", "bio12"), return.tbl = TRUE, quiet = TRUE, occs.testing.z = occs.testing.z) |> test_hist()
    evalplot.envSim.hist(e = e, ref.data = "occs", hist.bins = 50, return.tbl = TRUE, quiet = TRUE, occs.testing.z = occs.testing.z) |> test_hist()
    # with occs and bg data
    evalplot.envSim.hist(occs.z = occs.z, bg.z = bg.z, occs.grp = occs.grp, bg.grp = bg.grp, ref.data = "occs", return.tbl = TRUE, quiet = TRUE, occs.testing.z = occs.testing.z) |> test_hist()
    if(bg.sel == 1) {
      evalplot.envSim.hist(occs.z = occs.z, bg.z = bg.z, occs.grp = occs.grp, bg.grp = bg.grp, ref.data = "bg", return.tbl = TRUE, quiet = TRUE, occs.testing.z = occs.testing.z) |> test_hist() 
    }
    evalplot.envSim.hist(occs.z = occs.z, bg.z = bg.z, occs.grp = occs.grp, bg.grp = bg.grp, ref.data = "occs", envs.vars = c("bio1", "bio12"), return.tbl = TRUE, quiet = TRUE, occs.testing.z = occs.testing.z) |> test_hist()
    evalplot.envSim.hist(occs.z = occs.z, bg.z = bg.z, occs.grp = occs.grp, bg.grp = bg.grp, ref.data = "occs", hist.bins = 50, return.tbl = TRUE, quiet = TRUE, occs.testing.z = occs.testing.z) |> test_hist()
  
}

test_evalplot.envSim.map <- function(e, envs, occs.z, bg.z, occs.grp, bg.grp, bg.sel = 1, occs.testing.z = NULL) { 
  test_map <- function(i) {
    test_that("Outputs for evalplot.envSim.map have correct form", {
      if(inherits(i, "SpatRaster")) {
        if(is.null(occs.testing.z)) {
          expect_true(length(unique(e@occs.grp)) == terra::nlyr(i))
        }else{
          expect_true(2 == terra::nlyr(i))
        }
      }else{
        expect_true(ncol(i) == 4)
        expect_true(names(i)[1] == "x")
        expect_true(names(i)[2] == "y")
        expect_true(names(i)[3] == "ras")
      }
    })
  }
  # with ENMevaluation object
  evalplot.envSim.map(e = e, envs = envs, ref.data = "occs", return.tbl = TRUE, quiet = TRUE, occs.testing.z = occs.testing.z) |> test_map()
  evalplot.envSim.map(e = e, envs = envs, ref.data = "occs", return.ras = TRUE, quiet = TRUE, occs.testing.z = occs.testing.z) |> test_map()
  evalplot.envSim.map(e = e, envs = envs, ref.data = "occs", envs.vars = c("bio1","bio12"), return.tbl = TRUE, quiet = TRUE, occs.testing.z = occs.testing.z) |> test_map()
  evalplot.envSim.map(e = e, envs = envs, ref.data = "occs", envs.vars = c("bio1","bio12"), return.ras = TRUE, quiet = TRUE, occs.testing.z = occs.testing.z) |> test_map()
  # with buffer
  evalplot.envSim.map(e = e, envs = envs, ref.data = "occs", bb.buf = 5, return.tbl = TRUE, quiet = TRUE, occs.testing.z = occs.testing.z) |> test_map()
  evalplot.envSim.map(e = e, envs = envs, ref.data = "occs", bb.buf = 5, return.ras = TRUE, quiet = TRUE, occs.testing.z = occs.testing.z) |> test_map()
  # with occs and bg data
  evalplot.envSim.map(occs.z = occs.z, occs.grp = occs.grp, envs = envs, ref.data = "occs", bb.buf = 5, return.ras = TRUE, quiet = TRUE, occs.testing.z = occs.testing.z) |> test_map()
  if(bg.sel == 1) {
    evalplot.envSim.map(bg.z = bg.z, bg.grp = bg.grp, envs = envs, ref.data = "bg", bb.buf = 5, return.ras = TRUE, quiet = TRUE, occs.testing.z = occs.testing.z) |> test_map()
  }
  
}

test_evalplot.nulls <- function(ns) {
  
  test_nulls <- function(x, stats) {
    test_that("Outputs for evalplot.nulls have correct form", {
      expect_true(length(x) == 2)
      expect_true(all(names(x) == c("null.avgs", "empirical.results")))
      expect_true(ncol(x[[1]]) == 2)
      expect_true(ncol(x[[2]]) == 2)
      expect_true(all(names(x[[1]]) == c("metric", "avg")))
      expect_true(all(names(x[[2]]) == c("metric", "avg")))
      expect_true(nrow(x[[1]]) == ns@null.no.iter * length(stats))
      expect_true(all(unique(x[[1]]$metric) == stats))
      expect_true(all(unique(x[[2]]$metric) == stats))
    })
  }
  
  if(ns@null.partition.method == "none") {
    stat1 <- "auc.train"
    if(requireNamespace("ecospat", quietly = TRUE)) stat2 <- c("auc.train", "cbi.train") else stat2 <- "auc.train"
  }else{
    stat1 <- "auc.val"
    stat2 <- c("auc.val", "or.10p") 
  }
  
  for(i in c("histogram", "violin")) {
    # one metric 
    evalplot.nulls(ns, stat1, plot.type = i, facet.labels = NULL, metric.levels = NULL, return.tbl = TRUE) |> test_nulls(stat1)
    # two metrics 
    evalplot.nulls(ns, stat2, plot.type = i, facet.labels = NULL, metric.levels = NULL, return.tbl = TRUE) |> test_nulls(stat2)
    # one metric labels
    evalplot.nulls(ns, stat1, plot.type = i, facet.labels = paste0(stat1, "2"), metric.levels = NULL, return.tbl = TRUE) |> test_nulls(stat1)
    # two metrics labels
    evalplot.nulls(ns, stat2, plot.type = i, facet.labels = paste0(stat2, "2"), metric.levels = NULL, return.tbl = TRUE) |> test_nulls(stat2)
    # one metric levels
    evalplot.nulls(ns, stat1, plot.type = i, facet.labels = NULL, metric.levels = stat1, return.tbl = TRUE) |> test_nulls(stat1)
    # two metrics levels
    evalplot.nulls(ns, stat2, plot.type = i, facet.labels = NULL, metric.levels = stat2, return.tbl = TRUE) |> test_nulls(stat2)
  }
}
