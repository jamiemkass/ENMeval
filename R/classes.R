#' @importFrom methods setClass setGeneric setMethod setRefClass show
NULL

#' @title ENMevaluation class
#' @description An S4 class that contains the ENMevaluate results.
#' @author Jamie M. Kass, \email{jamie.m.kass@@gmail.com}, Bob Muscarella, \email{bob.muscarella@@gmail.com}
#' @slot algorithm character: algorithm used
#' @slot tune.settings data frame: settings that were tuned
#' @slot partition.method character: partition method used
#' @slot partition.settings list: partition settings used (i.e., value of *k* or aggregation factor)
#' @slot other.settings list: other modeling settings used (i.e., decisions about clamping, AUC diff calculation)
#' @slot results data frame: evaluation summary statistics
#' @slot results.partitions data frame: evaluation k-fold statistics
#' @slot models list: model objects
#' @slot predictions RasterStack: model predictions
#' @slot occs data frame: occurrence coordinates and predictor variable values used for model training
#' @slot occs.grp vector: partition groups for occurrence points
#' @slot bg data frame: background coordinates and predictor variable values used for model training
#' @slot bg.grp vector: partition groups for background points
#' @slot overlap list: matrices of pairwise niche overlap statistics
#' 
#' @export

# class slots match older ENMeval versions
ENMevaluation <- setClass("ENMevaluation",
                          slots = c(algorithm = 'character',
                                    tune.settings = 'data.frame',
                                    partition.method = 'character',
                                    partition.settings = 'list',
                                    other.settings = 'list',
                                    doClamp = 'logical',
                                    clamp.directions = 'list',
                                    results = 'data.frame',
                                    results.partitions = 'data.frame',
                                    models = 'list',
                                    variable.importance = 'list',
                                    predictions = 'RasterStack',
                                    taxon.name = 'character',
                                    occs = 'data.frame',
                                    occs.testing = 'data.frame',
                                    occs.grp = 'factor',
                                    bg = 'data.frame',
                                    bg.grp = 'factor',
                                    overlap = 'list',
                                    rmm = 'list'))

#' @title eval.algorithm generic for ENMevaluation object
#' @param x ENMevaluation object
#' @rdname eval.algorithm
#' @export
setGeneric("eval.algorithm", function(x) standardGeneric("eval.algorithm"))

#' @rdname eval.algorithm
setMethod("eval.algorithm", "ENMevaluation", function(x) x@algorithm)

#' @title eval.tune.settings generic for ENMevaluation object
#' @param x ENMevaluation object
#' @rdname eval.tune.settings
#' @export
setGeneric("eval.tune.settings", function(x) standardGeneric("eval.tune.settings"))

#' @rdname eval.tune.settings
setMethod("eval.tune.settings", "ENMevaluation", function(x) x@tune.settings)

#' @title eval.results generic for ENMevaluation object
#' @param x ENMevaluation object
#' @rdname eval.results
#' @export
setGeneric("eval.results", function(x) standardGeneric("eval.results"))

#' @rdname eval.results
setMethod("eval.results", "ENMevaluation", function(x) x@results)

#' @title eval.results.partitions generic for ENMevaluation object
#' @param x ENMevaluation object
#' @rdname eval.results.partitions
#' @export
setGeneric("eval.results.partitions", function(x) standardGeneric("eval.results.partitions"))

#' @rdname eval.results.partitions
setMethod("eval.results.partitions", "ENMevaluation", function(x) x@results.partitions)

#' @title eval.predictions generic for ENMevaluation object
#' @param x ENMevaluation object
#' @rdname eval.predictions
#' @export
setGeneric("eval.predictions", function(x) standardGeneric("eval.predictions"))

#' @rdname eval.predictions
setMethod("eval.predictions", "ENMevaluation", function(x) x@predictions)

#' @title eval.models generic for ENMevaluation object
#' @param x ENMevaluation object
#' @rdname eval.models
#' @export
setGeneric("eval.models", function(x) standardGeneric("eval.models"))

#' @rdname eval.models
setMethod("eval.models", "ENMevaluation", function(x) x@models)

#' @title eval.partition.method generic for ENMevaluation object
#' @param x ENMevaluation object
#' @rdname eval.partition.method
#' @export
setGeneric("eval.partition.method", function(x) standardGeneric("eval.partition.method"))

#' @rdname eval.partition.method
setMethod("eval.partition.method", "ENMevaluation", function(x) x@partition.method)

#' @title eval.partition.settings generic for ENMevaluation object
#' @param x ENMevaluation object
#' @rdname eval.partition.settings
#' @export
setGeneric("eval.partition.settings", function(x) standardGeneric("eval.partition.settings"))

#' @rdname eval.partition.settings
setMethod("eval.partition.settings", "ENMevaluation", function(x) x@partition.settings)

#' @title eval.other.settings generic for ENMevaluation object
#' @param x ENMevaluation object
#' @rdname eval.other.settings
#' @export
setGeneric("eval.other.settings", function(x) standardGeneric("eval.other.settings"))

#' @rdname eval.other.settings
setMethod("eval.other.settings", "ENMevaluation", function(x) x@other.settings)

#' @title eval.doClamp generic for ENMevaluation object
#' @param x ENMevaluation object
#' @rdname eval.doClamp
#' @export
setGeneric("eval.doClamp", function(x) standardGeneric("eval.doClamp"))

#' @rdname eval.doClamp
setMethod("eval.doClamp", "ENMevaluation", function(x) x@other.settings)

#' @title eval.clamp.directions generic for ENMevaluation object
#' @param x ENMevaluation object
#' @rdname eval.clamp.directions
#' @export
setGeneric("eval.clamp.directions", function(x) standardGeneric("eval.clamp.directions"))

#' @rdname eval.clamp.directions
setMethod("eval.clamp.directions", "ENMevaluation", function(x) x@other.settings)

#' @title eval.taxon.name generic for ENMevaluation object
#' @param x ENMevaluation object
#' @rdname eval.taxon.name
#' @export
setGeneric("eval.taxon.name", function(x) standardGeneric("eval.taxon.name"))

#' @rdname eval.taxon.name
setMethod("eval.taxon.name", "ENMevaluation", function(x) x@taxon.name)

#' @title eval.occs generic for ENMevaluation object
#' @param x ENMevaluation object
#' @rdname eval.occs
#' @export
setGeneric("eval.occs", function(x) standardGeneric("eval.occs"))

#' @rdname eval.occs
setMethod("eval.occs", "ENMevaluation", function(x) x@occs)

#' @title eval.occs.testing generic for ENMevaluation object
#' @param x ENMevaluation object
#' @rdname eval.occs.testing
#' @export
setGeneric("eval.occs.testing", function(x) standardGeneric("eval.occs.testing"))

#' @rdname eval.occs.testing
setMethod("eval.occs.testing", "ENMevaluation", function(x) x@occs.testing)

#' @title eval.occs.grp generic for ENMevaluation object
#' @param x ENMevaluation object
#' @rdname eval.occs.grp
#' @export
setGeneric("eval.occs.grp", function(x) standardGeneric("eval.occs.grp"))

#' @rdname eval.occs.grp
setMethod("eval.occs.grp", "ENMevaluation", function(x) x@occs.grp)

#' @title eval.bg generic for ENMevaluation object
#' @param x ENMevaluation object
#' @rdname eval.bg
#' @export
setGeneric("eval.bg", function(x) standardGeneric("eval.bg"))

#' @rdname eval.bg
setMethod("eval.bg", "ENMevaluation", function(x) x@bg)

#' @title eval.bg.grp generic for ENMevaluation object
#' @param x ENMevaluation object
#' @rdname eval.bg.grp
#' @export
setGeneric("eval.bg.grp", function(x) standardGeneric("eval.bg.grp"))

#' @rdname eval.bg.grp
setMethod("eval.bg.grp", "ENMevaluation", function(x) x@bg.grp)

#' @title eval.overlap generic for ENMevaluation object
#' @param x ENMevaluation object
#' @rdname eval.overlap
#' @export
setGeneric("eval.overlap", function(x) standardGeneric("eval.overlap"))

#' @rdname eval.overlap
setMethod("eval.overlap", "ENMevaluation", function(x) x@overlap)

#' @title eval.rmm generic for ENMevaluation object
#' @param x ENMevaluation object
#' @rdname eval.rmm
#' @export
setGeneric("eval.rmm", function(x) standardGeneric("eval.rmm"))

#' @rdname eval.rmm
setMethod("eval.rmm", "ENMevaluation", function(x) x@rmm)

#' @export
setMethod("show",
          signature = "ENMevaluation",
          definition = function(object) {
            cat("An object of class: ", class(object), "\n")
            cat(" taxon name: ", object@taxon.name, "\n")
            cat(" occurrence/background points: ", nrow(object@occs), '/', nrow(object@bg), "\n")
            cat(" partition method: ", object@partition.method, "\n")
            cat(" partition settings: ", ifelse(length(object@partition.settings) > 0, paste(names(object@partition.settings), unlist(object@partition.settings), sep = " = ", collapse = ", "), "none"), "\n")
            cat(" clamp: ", object@doClamp, "\n")
            cat("clamp directions: ", ifelse(length(object@clamp.directions) > 0, paste(sapply(1:2, function(x) paste0(names(object@clamp.directions[x]), ": ", paste(object@clamp.directions[[x]], collapse = ", "))), collapse = "\n"), "none"), "\n")
            cat(" algorithm: ", object@algorithm, "\n")
            cat(" tune settings (", paste0(names(object@tune.settings)[-ncol(object@tune.settings)], collapse = "_"), "): ", paste0(object@tune.settings[,ncol(object@tune.settings)], collapse = ", "), "\n")
            cat(" overlap: ", !is.null(object@overlap), "\n")
            cat("Refer to ?ENMevaluation for information on slots.", sep = "")
            invisible(NULL)
          })

#' @title ENMdetails class
#' @description An S4 class that details packages, functions, messages associated with a specific species distribution model (SDM) or ecological niche model (ENM). 
#' Objects of this class are generated by \code{ENMdetails()}. For examples, look in the package's R folder for scripts beginning with "enm" -- these are 
#' pre-made ENMdetails object specifications that work with ENMeval out of the box.
#' @author Jamie M. Kass, \email{jamie.m.kass@@gmail.com}
#' @slot name character: name of algorithm
#' @slot fun function: function that runs the algorithm
#' @slot msgs function: prints messages showing the package version number, etc., and those related to the input tuning parameters \code{tune.args}
#' @slot args function: returns the parameters needed to run the algorithm function
#' @slot predict function: specifies how to calculate a model prediction for a Raster* or a data frame
#' @slot ncoefs function: counts the number of non-zero model coefficients
#' @slot varimp function: generates a data frame of variable importance from the model object (if functionality is available)
#' @export

#' @export
ENMdetails <- setClass("ENMdetails",
                       slots = c(name = 'character',
                                 fun = 'function',
                                 errors = 'function',
                                 msgs = 'function',
                                 args = 'function',
                                 predict = 'function',
                                 ncoefs = 'function',
                                 varimp = 'function'))

#' ENMdetails <- function(name, fun, errors, msgs, args, predict, ncoefs, varimp) {
#'   new("ENMdetails", name = name, fun = fun, errors = errors, msgs = msgs, 
#'       args = args, predict = predict, ncoefs = ncoefs, varimp = varimp)
#' }

setGeneric("enm.name", function(x) standardGeneric("enm.name"))
setGeneric("enm.name<-", function(x, value) standardGeneric("enm.name<-"))
#' @export
setMethod("enm.name", "ENMdetails", function(x) x@name)
#' @export
setMethod("enm.name<-", "ENMdetails", function(x, value) {
  x@name <- value
  validObject(x)
  x
})

setGeneric("enm.fun", function(x) standardGeneric("enm.fun"))
setGeneric("enm.fun<-", function(x, value) standardGeneric("enm.fun<-"))
#' @export
setMethod("enm.fun", "ENMdetails", function(x) x@fun)
#' @export
setMethod("enm.fun<-", "ENMdetails", function(x, value) {
  x@fun <- value
  validObject(x)
  x
})

setGeneric("enm.errors", function(x) standardGeneric("enm.errors"))
setGeneric("enm.errors<-", function(x, value) standardGeneric("enm.errors<-"))
#' @export
setMethod("enm.errors", "ENMdetails", function(x) x@errors)
#' @export
setMethod("enm.errors<-", "ENMdetails", function(x, value) {
  x@errors <- value
  validObject(x)
  x
})

setGeneric("enm.msgs", function(x) standardGeneric("enm.msgs"))
setGeneric("enm.msgs<-", function(x, value) standardGeneric("enm.msgs<-"))
#' @export
setMethod("enm.msgs", "ENMdetails", function(x) x@msgs)
#' @export
setMethod("enm.msgs<-", "ENMdetails", function(x, value) {
  x@msgs <- value
  validObject(x)
  x
})

setGeneric("enm.args", function(x) standardGeneric("enm.args"))
setGeneric("enm.args<-", function(x, value) standardGeneric("enm.args<-"))
#' @export
setMethod("enm.args", "ENMdetails", function(x) x@args)
#' @export
setMethod("enm.args<-", "ENMdetails", function(x, value) {
  x@args <- value
  validObject(x)
  x
})

setGeneric("enm.predict", function(x) standardGeneric("enm.predict"))
setGeneric("enm.predict<-", function(x, value) standardGeneric("enm.predict<-"))
#' @export
setMethod("enm.predict", "ENMdetails", function(x) x@predict)
#' @export
setMethod("enm.predict<-", "ENMdetails", function(x, value) {
  x@predict <- value
  validObject(x)
  x
})

setGeneric("enm.ncoefs", function(x) standardGeneric("enm.ncoefs"))
setGeneric("enm.ncoefs<-", function(x, value) standardGeneric("enm.ncoefs<-"))
#' @export
setMethod("enm.ncoefs", "ENMdetails", function(x) x@ncoefs)
#' @export
setMethod("enm.ncoefs<-", "ENMdetails", function(x, value) {
  x@ncoefs <- value
  validObject(x)
  x
})

setGeneric("enm.varimp", function(x) standardGeneric("enm.varimp"))
setGeneric("enm.varimp<-", function(x, value) standardGeneric("enm.varimp<-"))
#' @export
setMethod("enm.varimp", "ENMdetails", function(x) x@varimp)
#' @export
setMethod("enm.varimp<-", "ENMdetails", function(x, value) {
  x@varimp <- value
  validObject(x)
  x
})

#' @export
setMethod("show",
          signature = "ENMdetails",
          definition = function(object) {
            cat("An object of class: ", class(object), "\n")
            cat(" Name: ", object@name, "\n")
            cat("Refer to ?ENMdetails for information on slots, and to the vignette for how to construct a custom object.", sep = "")
            invisible(NULL)
          })


#' @title ENMnull class
#' @description An S4 class that contains the ENMnulls results. 
#' @author Jamie M. Kass, \email{jamie.m.kass@@gmail.com}, Corentin Bohl, \email{corentinbohl@@gmail.com}
#' @slot null.algorithm character: algorithm used
#' @slot null.mod.settings data frame: model settings used
#' @slot null.partition.method character: partition method used
#' @slot null.partition.settings list: partition settings used (i.e., value of *k* or aggregation factor)
#' @slot null.other.settings list: other modeling settings used (i.e., decisions about clamping, AUC diff calculation)
#' @slot no.iter numeric: number of null model iterations
#' @slot null.results data frame: evaluation summary statistics for null models
#' @slot null.results.partitions data frame: evaluation k-fold statistics for null models
#' @slot emp.vs.null.results data frame: evaluation summary statistics for the empirical model, means for all null models, z-scores, and p-values
#' @slot emp.occs data frame: occurrence coordinates and predictor variable values used for model training (empirical model)
#' @slot emp.occs.grp vector: partition groups for occurrence points (empirical model)
#' @slot emp.bg data frame: background coordinates and predictor variable values used for model training (empirical model)
#' @slot emp.bg.grp vector: partition groups for background points (empirical model)
#' @export

# class slots match older ENMeval versions
ENMnull <- setClass("ENMnull",
                    slots = c(null.algorithm = 'character',
                              null.mod.settings = 'data.frame',
                              null.partition.method = 'character',
                              null.partition.settings = 'list',
                              null.other.settings = 'list',
                              no.iter = 'numeric',
                              null.results = 'data.frame',
                              null.results.partitions = 'data.frame',
                              emp.vs.null.results = 'data.frame',
                              emp.occs = 'data.frame',
                              emp.occs.grp = 'factor',
                              emp.bg = 'data.frame',
                              emp.bg.grp = 'factor'))

setGeneric("null.algorithm", function(x) standardGeneric("null.algorithm"))

#' @export
setMethod("null.algorithm", "ENMnull", function(x) x@null.algorithm)

setGeneric("null.mod.settings", function(x) standardGeneric("null.mod.settings"))
#' @export
setMethod("null.mod.settings", "ENMnull", function(x) x@null.mod.settings)

setGeneric("null.partition.method", function(x) standardGeneric("null.partition.method"))
#' @export
setMethod("null.partition.method", "ENMnull", function(x) x@null.partition.method)

setGeneric("null.partition.settings", function(x) standardGeneric("null.partition.settings"))
#' @export
setMethod("null.partition.settings", "ENMnull", function(x) x@null.partition.settings)

setGeneric("null.other.settings", function(x) standardGeneric("null.other.settings"))
#' @export
setMethod("null.other.settings", "ENMnull", function(x) x@null.other.settings)

setGeneric("no.iter", function(x) standardGeneric("no.iter"))
#' @export
setMethod("no.iter", "ENMnull", function(x) x@no.iter)

setGeneric("null.results", function(x) standardGeneric("null.results"))
#' @export
setMethod("null.results", "ENMnull", function(x) x@null.results)

setGeneric("null.results.partitions", function(x) standardGeneric("null.results.partitions"))
#' @export
setMethod("null.results.partitions", "ENMnull", function(x) x@null.results.partitions)

setGeneric("emp.vs.null.results", function(x) standardGeneric("emp.vs.null.results"))
#' @export
setMethod("emp.vs.null.results", "ENMnull", function(x) x@emp.vs.null.results)

setGeneric("emp.occs", function(x) standardGeneric("emp.occs"))
#' @export
setMethod("emp.occs", "ENMnull", function(x) x@emp.occs)

setGeneric("emp.occs.grp", function(x) standardGeneric("emp.occs.grp"))
#' @export
setMethod("emp.occs.grp", "ENMnull", function(x) x@emp.occs.grp)

setGeneric("emp.bg", function(x) standardGeneric("emp.bg"))
#' @export
setMethod("emp.bg", "ENMnull", function(x) x@emp.bg)

setGeneric("emp.bg.grp", function(x) standardGeneric("emp.bg.grp"))
#' @export
setMethod("emp.bg.grp", "ENMnull", function(x) x@emp.bg.grp)


#' @export
setMethod("show",
          signature = "ENMnull",
          definition = function(object) {
            cat("An object of class: ", class(object), "\n")
            cat(" no. iterations: ", object@no.iter, "\n")
            cat(" empirical occurrence/background points: ", nrow(object@emp.occs), '/', nrow(object@emp.bg), "\n")
            cat(" partition method: ", object@null.partition.method, "\n")
            cat(" partition settings: ", paste(names(object@null.partition.settings), unlist(object@null.partition.settings), sep = " = ", collapse = ", "), "\n")
            cat(" other settings: ", paste(names(object@null.other.settings), unlist(object@null.other.settings), sep = " = ", collapse = ", "), "\n")
            cat(" algorithm: ", object@null.algorithm, "\n")
            cat(" model settings: \n")
            print(object@null.mod.settings[,-ncol(object@null.mod.settings)], row.names = FALSE)
            cat("Refer to ?ENMnull for information on slots.", sep = "")
            invisible(NULL)
          })

