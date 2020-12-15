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
    if(type == "swd") {
      expect_true(raster::nlayers(e@predictions) == 0)  
    }else{
      expect_true(raster::nlayers(e@predictions) > 0)  
    }
    expect_true(!is.null(e@occs))
    expect_true(!is.null(e@occs.grp))
    expect_true(!is.null(e@bg))
    expect_true(!is.null(e@bg.grp))
    expect_true(!is.null(e@overlap))
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
      if(parts == "jackknife") {
        expect_true(sum(is.na(e@results.partitions)) == nrow(e@results.partitions))
      }else{
        expect_true(sum(is.na(e@results.partitions)) == 0)
      }
    }
  })
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
    expect_true(!is.null(ns@no.iter))
    expect_true(!is.null(ns@null.results))
    expect_true(!is.null(ns@null.results.partitions))
    expect_true(!is.null(ns@emp.vs.null.results))
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
    expect_true(ns@no.iter == no.iter)
    # number of rows in results table
    expect_true(nrow(ns@null.results) == no.iter)
    # number of rows in results table for partitions
    if(ns@null.partition.method == "none") {
      expect_true(nrow(ns@null.results.partitions) == 0)
    }else{
      expect_true(nrow(ns@null.results.partitions) == no.iter * nparts.occs)  
    }
    
    # number of rows in empirical vs null results table
    expect_true(nrow(ns@emp.vs.null.results) == 6)
    # there should only be two NA values for this table: read.sd for auc.train and cbi.train
    if(parts == "jackknife") {
      expect_true(sum(is.na(ns@emp.vs.null.results[2,])) == 3)
      expect_true(sum(is.na(ns@emp.vs.null.results[,6])) == 6) 
    }else if(parts == "testing") {
      expect_true(sum(is.na(ns@emp.vs.null.results[2,])) == 7) 
    }else{
      expect_true(sum(is.na(ns@emp.vs.null.results[2,])) == 2)  
    }
    # check that tables match
    expect_true(all(ns@emp.occs == e@occs))
    expect_true(all(ns@emp.bg == e@bg))
    expect_true(all(ns@emp.occs.grp == e@occs.grp))
  })
}

#' @title Unit tests for ENMevaluation plotting functions
#' @description All parameters are self-explanatory except the following.
#' Anything with the prefix ".z" is a data frame with latitude, longitude,
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
    expect_true(ncol(x) == ncol(e@tune.settings) + y)
    expect_true(nrow(x) == nrow(e@tune.settings) * length(stats))
    expect_true(all(unique(x$metric) == stats))
    n <- (ncol(e@tune.settings)+1):(ncol(x))
    expect_true(all(names(x)[n] == z))
  }
  
  if(e@partition.method == "none") {
    stat1 <- "auc.train"
    stat2 <- c("auc.train", "cbi.train")
  }else{
    stat1 <- "auc.val"
    stat2 <- c("auc.val", "or.10p") 
  }
  
  # defaults
  evalplot.stats(e, stats = stat1, x.var = "rm", color = "fc", dodge = FALSE, error.bars = FALSE, facet.labels = NULL, metric.levels = NULL, return.tbl = TRUE) %>% test_stats(stat1)
  # two stats
  evalplot.stats(e, stats = stat2, x.var = "rm", color = "fc", dodge = FALSE, error.bars = FALSE, facet.labels = NULL, metric.levels = NULL, return.tbl = TRUE) %>% test_stats(stat2)
  # dodge
  evalplot.stats(e, stats = stat1, x.var = "rm", color = "fc", dodge = TRUE, error.bars = FALSE, facet.labels = NULL, metric.levels = NULL, return.tbl = TRUE) %>% test_stats(stat1)
  # error bars
  evalplot.stats(e, stats = stat1, x.var = "rm", color = "fc", dodge = FALSE, error.bars = TRUE, facet.labels = NULL, metric.levels = NULL, return.tbl = TRUE) %>% test_stats(stat1)
  # facet labels
  evalplot.stats(e, stats = stat1, x.var = "rm", color = "fc", dodge = FALSE, error.bars = TRUE, facet.labels = paste0(stat1), metric.levels = NULL, return.tbl = TRUE) %>% test_stats(stat1)
  evalplot.stats(e, stats = stat2, x.var = "rm", color = "fc", dodge = FALSE, error.bars = TRUE, facet.labels = paste0(stat2), metric.levels = NULL, return.tbl = TRUE) %>% test_stats(stat2)
  # metric levels
  evalplot.stats(e, stats = stat1, x.var = "rm", color = "fc", dodge = FALSE, error.bars = TRUE, facet.labels = NULL, metric.levels = stat1, return.tbl = TRUE) %>% test_stats(stat1)
  evalplot.stats(e, stats = stat2, x.var = "rm", color = "fc", dodge = FALSE, error.bars = TRUE, facet.labels = NULL, metric.levels = rev(stat2), return.tbl = TRUE) %>% test_stats(stat2)
}

test_evalplot.envSim.hist <- function(e, occs.z, bg.z, occs.grp, bg.grp, bg.sel = 1, occs.testing.z = NULL) {
  all_tests_hist <- function(sim.type) {
    test_hist <- function(i) {
      expect_true(ncol(i) == 2)
      expect_true(names(i)[1] == "partition")
      expect_true(names(i)[2] == sim.type)
    }
    # with ENMevaluation object
    evalplot.envSim.hist(e = e, ref.data = "occs", sim.type = sim.type, categoricals = "biome", return.tbl = TRUE, quiet = TRUE, occs.testing.z = occs.testing.z) %>% test_hist()
    if(bg.sel == 1) {
      evalplot.envSim.hist(e = e, ref.data = "bg", sim.type = sim.type, categoricals = "biome", return.tbl = TRUE, quiet = TRUE, occs.testing.z = occs.testing.z) %>% test_hist()
    }
    evalplot.envSim.hist(e = e, ref.data = "occs", sim.type = sim.type, categoricals = "biome", envs.var = c("bio1", "bio12"), return.tbl = TRUE, quiet = TRUE, occs.testing.z = occs.testing.z) %>% test_hist()
    evalplot.envSim.hist(e = e, ref.data = "occs", sim.type = sim.type, categoricals = "biome", hist.bins = 50, return.tbl = TRUE, quiet = TRUE, occs.testing.z = occs.testing.z) %>% test_hist()
    # with occs and bg data
    evalplot.envSim.hist(occs.z = occs.z, bg.z = bg.z, occs.grp = occs.grp, bg.grp = bg.grp, ref.data = "occs", sim.type = sim.type, categoricals = "biome", return.tbl = TRUE, quiet = TRUE, occs.testing.z = occs.testing.z) %>% test_hist()
    if(bg.sel == 1) {
      evalplot.envSim.hist(occs.z = occs.z, bg.z = bg.z, occs.grp = occs.grp, bg.grp = bg.grp, ref.data = "bg", sim.type = sim.type, categoricals = "biome", return.tbl = TRUE, quiet = TRUE, occs.testing.z = occs.testing.z) %>% test_hist() 
    }
    evalplot.envSim.hist(occs.z = occs.z, bg.z = bg.z, occs.grp = occs.grp, bg.grp = bg.grp, ref.data = "occs", sim.type = sim.type, categoricals = "biome", envs.var = c("bio1", "bio12"), return.tbl = TRUE, quiet = TRUE, occs.testing.z = occs.testing.z) %>% test_hist()
    evalplot.envSim.hist(occs.z = occs.z, bg.z = bg.z, occs.grp = occs.grp, bg.grp = bg.grp, ref.data = "occs", sim.type = sim.type, categoricals = "biome", hist.bins = 50, return.tbl = TRUE, quiet = TRUE, occs.testing.z = occs.testing.z) %>% test_hist()
  }
  
    all_tests_hist("mess")
    all_tests_hist("most_diff")
    all_tests_hist("most_sim")  
}

test_evalplot.envSim.map <- function(e, envs, occs.z, bg.z, occs.grp, bg.grp, bg.sel = 1, occs.testing.z = NULL) { 
  all_tests_map <- function(sim.type) {
    test_map <- function(i) {
      if(inherits(i, "Raster")) {
        if(is.null(occs.testing.z)) {
          expect_true(length(unique(e@occs.grp)) == raster::nlayers(i))
        }else{
          expect_true(2 == raster::nlayers(i))
        }
      }else{
        expect_true(ncol(i) == 4)
        expect_true(names(i)[1] == "x")
        expect_true(names(i)[2] == "y")
        expect_true(names(i)[3] == "ras")
        expect_true(names(i)[4] == sim.type)  
      }
    }
    # with ENMevaluation object
    evalplot.envSim.map(e = e, envs = envs, ref.data = "occs", sim.type = sim.type, categoricals = "biome", return.tbl = TRUE, quiet = TRUE, occs.testing.z = occs.testing.z) %>% test_map()
    evalplot.envSim.map(e = e, envs = envs, ref.data = "occs", sim.type = sim.type, categoricals = "biome", return.ras = TRUE, quiet = TRUE, occs.testing.z = occs.testing.z) %>% test_map()
    evalplot.envSim.map(e = e, envs = envs, ref.data = "occs", sim.type = sim.type, categoricals = "biome", envs.var = c("bio1","bio12"), return.tbl = TRUE, quiet = TRUE, occs.testing.z = occs.testing.z) %>% test_map()
    evalplot.envSim.map(e = e, envs = envs, ref.data = "occs", sim.type = sim.type, categoricals = "biome", envs.var = c("bio1","bio12"), return.ras = TRUE, quiet = TRUE, occs.testing.z = occs.testing.z) %>% test_map()
    # with buffer
    evalplot.envSim.map(e = e, envs = envs, ref.data = "occs", sim.type = sim.type, categoricals = "biome", bb.buf = 5, return.tbl = TRUE, quiet = TRUE, occs.testing.z = occs.testing.z) %>% test_map()
    evalplot.envSim.map(e = e, envs = envs, ref.data = "occs", sim.type = sim.type, categoricals = "biome", bb.buf = 5, return.ras = TRUE, quiet = TRUE, occs.testing.z = occs.testing.z) %>% test_map()
    # with occs and bg data
    evalplot.envSim.map(occs.z = occs.z, occs.grp = occs.grp, envs = envs, ref.data = "occs", sim.type = sim.type, categoricals = "biome", bb.buf = 5, return.ras = TRUE, quiet = TRUE, occs.testing.z = occs.testing.z) %>% test_map()
    if(bg.sel == 1) {
      evalplot.envSim.map(bg.z = bg.z, bg.grp = bg.grp, envs = envs, ref.data = "bg", sim.type = sim.type, categoricals = "biome", bb.buf = 5, return.ras = TRUE, quiet = TRUE, occs.testing.z = occs.testing.z) %>% test_map()
    }
  }
  
    all_tests_map("mess")
    all_tests_map("most_diff")
    all_tests_map("most_sim")
}

test_evalplot.nulls <- function(ns) {
  
  test_nulls <- function(x, stats) {
    expect_true(length(x) == 2)
    expect_true(all(names(x) == c("null.avgs", "empirical.results")))
    expect_true(ncol(x[[1]]) == 2)
    expect_true(ncol(x[[2]]) == 2)
    expect_true(all(names(x[[1]]) == c("metric", "avg")))
    expect_true(all(names(x[[2]]) == c("metric", "avg")))
    expect_true(nrow(x[[1]]) == ns@no.iter * length(stats))
    expect_true(all(unique(x[[1]]$metric) == stats))
    expect_true(all(unique(x[[2]]$metric) == stats))
  }
  
  if(ns@null.partition.method == "none") {
    stat1 <- "auc.train"
    stat2 <- c("auc.train", "cbi.train")
  }else{
    stat1 <- "auc.val"
    stat2 <- c("auc.val", "or.10p") 
  }
  
  for(i in c("histogram", "violin")) {
    # one metric 
    evalplot.nulls(ns, stat1, plot.type = i, facet.labels = NULL, metric.levels = NULL, return.tbl = TRUE) %>% test_nulls(stat1)
    # two metrics 
    evalplot.nulls(ns, stat2, plot.type = i, facet.labels = NULL, metric.levels = NULL, return.tbl = TRUE) %>% test_nulls(stat2)
    # one metric labels
    evalplot.nulls(ns, stat1, plot.type = i, facet.labels = paste0(stat1, "2"), metric.levels = NULL, return.tbl = TRUE) %>% test_nulls(stat1)
    # two metrics labels
    evalplot.nulls(ns, stat2, plot.type = i, facet.labels = paste0(stat2, "2"), metric.levels = NULL, return.tbl = TRUE) %>% test_nulls(stat2)
    # one metric levels
    evalplot.nulls(ns, stat1, plot.type = i, facet.labels = NULL, metric.levels = stat1, return.tbl = TRUE) %>% test_nulls(stat1)
    # two metrics levels
    evalplot.nulls(ns, stat2, plot.type = i, facet.labels = NULL, metric.levels = stat2, return.tbl = TRUE) %>% test_nulls(stat2)
  }
}
