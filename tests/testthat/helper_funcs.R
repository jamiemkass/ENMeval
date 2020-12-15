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

test_ENMnullSims <- function(e, ns, no.iter, alg, parts, mod.settings, nparts.occs, nparts.bg, n.sims, type = "") {
  mod.settings.tbl <- expand.grid(mod.settings)
  test_that("ENMnullSims object and slots exist", {
    expect_true(!is.null(ns))
    expect_true(!is.null(ns@null.algorithm))
    expect_true(!is.null(ns@null.mod.settings))
    expect_true(!is.null(ns@null.partition.method))
    expect_true(!is.null(ns@null.partition.settings))
    expect_true(!is.null(ns@null.other.settings))
    expect_true(!is.null(ns@no.iter))
    expect_true(!is.null(ns@null.results))
    expect_true(!is.null(ns@null.results.partitions))
    expect_true(!is.null(ns@real.vs.null.results))
    expect_true(!is.null(ns@real.occs))
    expect_true(!is.null(ns@real.occs.grp))
    expect_true(!is.null(ns@real.bg))
    expect_true(!is.null(ns@real.bg.grp))
  })  
  
  test_that("Data in ENMnullSims object slots have correct form", {
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
    
    # number of rows in real vs null results table
    expect_true(nrow(ns@real.vs.null.results) == 6)
    # there should only be two NA values for this table: read.sd for auc.train and cbi.train
    if(parts == "jackknife") {
      expect_true(sum(is.na(ns@real.vs.null.results[2,])) == 3)
      expect_true(sum(is.na(ns@real.vs.null.results[,6])) == 6) 
    }else if(parts == "testing") {
      expect_true(sum(is.na(ns@real.vs.null.results[2,])) == 7) 
    }else{
      expect_true(sum(is.na(ns@real.vs.null.results[2,])) == 2)  
    }
    # check that tables match
    expect_true(all(ns@real.occs == e@occs))
    expect_true(all(ns@real.bg == e@bg))
    expect_true(all(ns@real.occs.grp == e@occs.grp))
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

test_evalPlots <- function(e, envs, occs.z, bg.z, occs.grp, bg.grp, plot.sel = c("hist", "map"), bg.sel = 1, occs.testing.z = NULL) {
  
  test_histPlot <- function(sim.type) {
    test_i <- function(i) {
      expect_true(ncol(i) == 2)
      expect_true(names(i)[1] == "partition")
      expect_true(names(i)[2] == sim.type)
    }
    # with ENMevaluation object
    hist.tbl1 <- evalplot.envSim.hist(e = e, ref.data = "occs", sim.type = sim.type, categoricals = "biome", return.tbl = TRUE, quiet = TRUE, occs.testing.z = occs.testing.z)
    test_i(hist.tbl1)
    if(bg.sel == 1) {
      hist.tbl2 <- evalplot.envSim.hist(e = e, ref.data = "bg", sim.type = sim.type, categoricals = "biome", return.tbl = TRUE, quiet = TRUE, occs.testing.z = occs.testing.z)
      test_i(hist.tbl2)  
    }
    hist.tbl3 <- evalplot.envSim.hist(e = e, ref.data = "occs", sim.type = sim.type, categoricals = "biome", envs.var = c("bio1", "bio12"), return.tbl = TRUE, quiet = TRUE, occs.testing.z = occs.testing.z)
    test_i(hist.tbl3)
    hist.tbl4 <- evalplot.envSim.hist(e = e, ref.data = "occs", sim.type = sim.type, categoricals = "biome", hist.bins = 50, return.tbl = TRUE, quiet = TRUE, occs.testing.z = occs.testing.z)
    test_i(hist.tbl4)
    # with occs and bg data
    hist.tbl5 <- evalplot.envSim.hist(occs.z = occs.z, bg.z = bg.z, occs.grp = occs.grp, bg.grp = bg.grp, ref.data = "occs", sim.type = sim.type, categoricals = "biome", return.tbl = TRUE, quiet = TRUE, occs.testing.z = occs.testing.z)
    test_i(hist.tbl5)
    if(bg.sel == 1) {
      hist.tbl6 <- evalplot.envSim.hist(occs.z = occs.z, bg.z = bg.z, occs.grp = occs.grp, bg.grp = bg.grp, ref.data = "bg", sim.type = sim.type, categoricals = "biome", return.tbl = TRUE, quiet = TRUE, occs.testing.z = occs.testing.z)
      test_i(hist.tbl6)
    }
    hist.tbl7 <- evalplot.envSim.hist(occs.z = occs.z, bg.z = bg.z, occs.grp = occs.grp, bg.grp = bg.grp, ref.data = "occs", sim.type = sim.type, categoricals = "biome", envs.var = c("bio1", "bio12"), return.tbl = TRUE, quiet = TRUE, occs.testing.z = occs.testing.z)
    test_i(hist.tbl7)
    hist.tbl8 <- evalplot.envSim.hist(occs.z = occs.z, bg.z = bg.z, occs.grp = occs.grp, bg.grp = bg.grp, ref.data = "occs", sim.type = sim.type, categoricals = "biome", hist.bins = 50, return.tbl = TRUE, quiet = TRUE, occs.testing.z = occs.testing.z)
    test_i(hist.tbl8)
  }
  
  if("hist" %in% plot.sel) {
    test_histPlot("mess")
    test_histPlot("most_diff")
    test_histPlot("most_sim")  
  }
  
  test_map <- function(sim.type) {
    test_i <- function(i) {
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
    map.tbl1 <- evalplot.envSim.map(e = e, envs = envs, ref.data = "occs", sim.type = sim.type, categoricals = "biome", return.tbl = TRUE, quiet = TRUE, occs.testing.z = occs.testing.z)
    test_i(map.tbl1)
    map.tbl2 <- evalplot.envSim.map(e = e, envs = envs, ref.data = "occs", sim.type = sim.type, categoricals = "biome", return.ras = TRUE, quiet = TRUE, occs.testing.z = occs.testing.z)
    test_i(map.tbl2)
    map.tbl3 <- evalplot.envSim.map(e = e, envs = envs, ref.data = "occs", sim.type = sim.type, categoricals = "biome", envs.var = c("bio1","bio12"), return.tbl = TRUE, quiet = TRUE, occs.testing.z = occs.testing.z)
    test_i(map.tbl3)
    map.tbl4 <- evalplot.envSim.map(e = e, envs = envs, ref.data = "occs", sim.type = sim.type, categoricals = "biome", envs.var = c("bio1","bio12"), return.ras = TRUE, quiet = TRUE, occs.testing.z = occs.testing.z)
    test_i(map.tbl4)
    # with buffer
    map.tbl5 <- evalplot.envSim.map(e = e, envs = envs, ref.data = "occs", sim.type = sim.type, categoricals = "biome", bb.buf = 5, return.tbl = TRUE, quiet = TRUE, occs.testing.z = occs.testing.z)
    test_i(map.tbl5)
    map.tbl6 <- evalplot.envSim.map(e = e, envs = envs, ref.data = "occs", sim.type = sim.type, categoricals = "biome", bb.buf = 5, return.ras = TRUE, quiet = TRUE, occs.testing.z = occs.testing.z)
    test_i(map.tbl6)
    # with occs and bg data
    map.tbl6 <- evalplot.envSim.map(occs.z = occs.z, occs.grp = occs.grp, envs = envs, ref.data = "occs", sim.type = sim.type, categoricals = "biome", bb.buf = 5, return.ras = TRUE, quiet = TRUE, occs.testing.z = occs.testing.z)
    test_i(map.tbl6)
    if(bg.sel == 1) {
      map.tbl6 <- evalplot.envSim.map(bg.z = bg.z, bg.grp = bg.grp, envs = envs, ref.data = "bg", sim.type = sim.type, categoricals = "biome", bb.buf = 5, return.ras = TRUE, quiet = TRUE, occs.testing.z = occs.testing.z)
      test_i(map.tbl6)  
    }
  }
  
  if("map" %in% plot.sel) {
    test_map("mess")
    test_map("most_diff")
    test_map("most_sim")
  }
}
