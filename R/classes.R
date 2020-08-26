#' @importFrom methods setClass setGeneric setMethod setRefClass
NULL

#' @title ENMevaluation class
#' @description An S4 class that contains the ENMevaluate results.
#' @author Jamie M. Kass, \email{jamie.m.kass@@gmail.com}, Bob Muscarella, \email{bob.muscarella@@gmail.com}
#' @slot algorithm character of algorithm used
#' @slot tune.settings data.frame of settings that were tuned
#' @slot partition.method character of partition method used
#' @slot partition.settings list of partition settings used (i.e., value of *k* or aggregation factor)
#' @slot other.settings list of other modeling settings used (i.e., decisions about clamping, AUC diff calculation)
#' @slot results data.frame of evaluation summary statistics
#' @slot results.grp data.frame of evaluation k-fold statistics
#' @slot models list of model objects
#' @slot predictions RasterStack of model predictions
#' @slot occs data.frame of occurrence coordinates and predictor variable values used for model training
#' @slot occs.grp vector of partition groups for occurrence points
#' @slot bg data.frame of background coordinates and predictor variable values used for model training
#' @slot bg.grp vector of partition groups for background points
#' @slot overlap list of matrices of pairwise niche overlap statistics
#' @export

# class slots match older ENMeval versions
ENMevaluation <- setClass("ENMevaluation",
                          slots = c(algorithm = 'character',
                                  tune.settings = 'data.frame',
                                  partition.method = 'character',
                                  partition.settings = 'list',
                                  other.settings = 'list',
                                  results = 'data.frame',
                                  results.grp = 'data.frame',
                                  models = 'list',
                                  predictions = 'RasterStack',
                                  taxon.name = 'character',
                                  occs = 'data.frame',
                                  occs.grp = 'factor',
                                  bg = 'data.frame',
                                  bg.grp = 'factor',
                                  overlap = 'list',
                                  rmm = 'list'))

setGeneric("eval.algorithm", function(x) standardGeneric("eval.algorithm"))
#' @export
setMethod("eval.algorithm", "ENMevaluation", function(x) x@algorithm)

setGeneric("eval.tune.settings", function(x) standardGeneric("eval.tune.settings"))
#' @export
setMethod("eval.tune.settings", "ENMevaluation", function(x) x@tune.settings)

setGeneric("eval.results", function(x) standardGeneric("eval.results"))
#' @export
setMethod("eval.results", "ENMevaluation", function(x) x@results)

setGeneric("eval.results.grp", function(x) standardGeneric("eval.results.grp"))
#' @export
setMethod("eval.results.grp", "ENMevaluation", function(x) x@results.grp)

setGeneric("eval.predictions", function(x) standardGeneric("eval.predictions"))
#' @export
setMethod("eval.predictions", "ENMevaluation", function(x) x@predictions)

setGeneric("eval.models", function(x) standardGeneric("eval.models"))
#' @export
setMethod("eval.models", "ENMevaluation", function(x) x@models)

setGeneric("eval.partition.method", function(x) standardGeneric("eval.partition.method"))
#' @export
setMethod("eval.partition.method", "ENMevaluation", function(x) x@partition.method)

setGeneric("eval.partition.settings", function(x) standardGeneric("eval.partition.settings"))
#' @export
setMethod("eval.partition.settings", "ENMevaluation", function(x) x@partition.settings)

setGeneric("eval.other.settings", function(x) standardGeneric("eval.other.settings"))
#' @export
setMethod("eval.other.settings", "ENMevaluation", function(x) x@other.settings)

setGeneric("eval.taxon.name", function(x) standardGeneric("eval.taxon.name"))
#' @export
setMethod("eval.taxon.name", "ENMevaluation", function(x) x@taxon.name)

setGeneric("eval.occs", function(x) standardGeneric("eval.occs"))
#' @export
setMethod("eval.occs", "ENMevaluation", function(x) x@occs)

setGeneric("eval.occs.grp", function(x) standardGeneric("eval.occs.grp"))
#' @export
setMethod("eval.occs.grp", "ENMevaluation", function(x) x@occs.grp)

setGeneric("eval.bg", function(x) standardGeneric("eval.bg"))
#' @export
setMethod("eval.bg", "ENMevaluation", function(x) x@bg)

setGeneric("eval.bg.grp", function(x) standardGeneric("eval.bg.grp"))
#' @export
setMethod("eval.bg.grp", "ENMevaluation", function(x) x@bg.grp)

setGeneric("eval.overlap", function(x) standardGeneric("eval.overlap"))
#' @export
setMethod("eval.overlap", "ENMevaluation", function(x) x@overlap)

setGeneric("eval.rmm", function(x) standardGeneric("eval.rmm"))
#' @export
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
#' @slot name character of name of model
#' @slot fun function that runs the model
#' @slot pkgs vector of package names needed to run the model function
#' @slot msgs function that prints messages showing the package version number, etc., and those related to the input tuning parameters \code{tune.args}
#' @slot args function specifying the parameters needed to run the model function
#' @slot evaluate function specifying how to calculate evaluation statistics
#' @slot train function specifying how to calculate training evaluation statistics
#' @slot validate function specifying how to calculate validation evaluation statistics
#' @slot predict function specifying how to calculate a model prediction for a Raster* or a data frame
#' @slot nparams function specifying how to measure the number of model coefficients
#' @export

#' @export
ENMdetails <- setClass("ENMdetails",
                       slots = c(name = 'character',
                                 fun = 'function',
                                 pkgs = 'character',
                                 msgs = 'function',
                                 args = 'function',
                                 evaluate = 'function',
                                 train = 'function',
                                 validate = 'function',
                                 predict = 'function',
                                 nparams = 'function'))
#' @export
ENMdetails <- function(name, fun, pkgs, msgs, args, evaluate, train, validate, predict, nparams) {
  new("ENMdetails", name = name, fun = fun, pkgs = pkgs, msgs = msgs, args = args, evaluate = evaluate,
      train = train, validate = validate, predict = predict, nparams = nparams)
}

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

setGeneric("enm.pkgs", function(x) standardGeneric("enm.pkgs"))
setGeneric("enm.pkgs<-", function(x, value) standardGeneric("enm.pkgs<-"))
#' @export
setMethod("enm.pkgs", "ENMdetails", function(x) x@pkgs)
#' @export
setMethod("enm.pkgs<-", "ENMdetails", function(x, value) {
  x@pkgs <- value
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

setGeneric("enm.evaluate", function(x) standardGeneric("enm.evaluate"))
setGeneric("enm.evaluate<-", function(x, value) standardGeneric("enm.evaluate<-"))
#' @export
setMethod("enm.evaluate", "ENMdetails", function(x) x@train)
#' @export
setMethod("enm.evaluate<-", "ENMdetails", function(x, value) {
  x@train <- value
  validObject(x)
  x
})

setGeneric("enm.train", function(x) standardGeneric("enm.train"))
setGeneric("enm.train<-", function(x, value) standardGeneric("enm.train<-"))
#' @export
setMethod("enm.train", "ENMdetails", function(x) x@train)
#' @export
setMethod("enm.train<-", "ENMdetails", function(x, value) {
  x@train <- value
  validObject(x)
  x
})

setGeneric("enm.validate", function(x) standardGeneric("enm.validate"))
setGeneric("enm.validate<-", function(x, value) standardGeneric("enm.validate<-"))
#' @export
setMethod("enm.validate", "ENMdetails", function(x) x@validate)
#' @export
setMethod("enm.validate<-", "ENMdetails", function(x, value) {
  x@validate <- value
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

setGeneric("enm.nparam", function(x) standardGeneric("enm.nparam"))
setGeneric("enm.nparam<-", function(x, value) standardGeneric("enm.nparam<-"))
#' @export
setMethod("enm.nparam", "ENMdetails", function(x) x@nparam)
#' @export
setMethod("enm.nparam<-", "ENMdetails", function(x, value) {
  x@nparam <- value
  validObject(x)
  x
})

#' @export
setMethod("show",
          signature = "ENMdetails",
          definition = function(object) {
            cat("An object of class: ", class(object), "\n")
            cat(" Name: ", object@name, "\n")
            cat(" Packages: ", paste(object@pkgs, sep = ", "), "\n")
            cat("Refer to ?ENMdetails for information on slots, and to the vignette for how to construct a custom object.", sep = "")
            invisible(NULL)
          })


#' @title ENMnull class
#' @description An S4 class that contains the ENMnullSims results. 
#' @author Jamie M. Kass, \email{jamie.m.kass@@gmail.com}, Corentin Bohl, \email{corentinbohl@@gmail.com}
#' @slot null.algorithm character of algorithm used
#' @slot null.mod.settings data.frame of model settings used
#' @slot null.partition.method character of partition method used
#' @slot null.partition.settings list of partition settings used (i.e., value of *k* or aggregation factor)
#' @slot null.other.settings list of other modeling settings used (i.e., decisions about clamping, AUC diff calculation)
#' @slot no.iter numeric of number of null model iterations
#' @slot null.results data.frame of evaluation summary statistics for null models
#' @slot null.results.grp data.frame of evaluation k-fold statistics for null models
#' @slot real.vs.null.results data.frame of evaluation summary statistics for the real model, means for all null models, z-scores, and p-values
#' @slot real.occs data.frame of occurrence coordinates and predictor variable values used for model training (real model)
#' @slot real.occs.grp vector of partition groups for occurrence points
#' @slot real.bg data.frame of background coordinates and predictor variable values used for model training
#' @slot real.bg.grp vector of partition groups for background points
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
                                    null.results.grp = 'data.frame',
                                    real.vs.null.results = 'data.frame',
                                    real.occs = 'data.frame',
                                    real.occs.grp = 'factor',
                                    real.bg = 'data.frame',
                                    real.bg.grp = 'factor'))

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

setGeneric("null.results.grp", function(x) standardGeneric("null.results.grp"))
#' @export
setMethod("null.results.grp", "ENMnull", function(x) x@null.results.grp)

setGeneric("real.vs.null.results", function(x) standardGeneric("real.vs.null.results"))
#' @export
setMethod("real.vs.null.results", "ENMnull", function(x) x@real.vs.null.results)

setGeneric("real.occs", function(x) standardGeneric("real.occs"))
#' @export
setMethod("real.occs", "ENMnull", function(x) x@real.occs)

setGeneric("real.occs.grp", function(x) standardGeneric("real.occs.grp"))
#' @export
setMethod("real.occs.grp", "ENMnull", function(x) x@real.occs.grp)

setGeneric("real.bg", function(x) standardGeneric("real.bg"))
#' @export
setMethod("real.bg", "ENMnull", function(x) x@real.bg)

setGeneric("real.bg.grp", function(x) standardGeneric("real.bg.grp"))
#' @export
setMethod("real.bg.grp", "ENMnull", function(x) x@real.bg.grp)


#' @export
setMethod("show",
          signature = "ENMnull",
          definition = function(object) {
            cat("An object of class: ", class(object), "\n")
            cat(" no. iterations: ", object@no.iter, "\n")
            cat(" real occurrence/background points: ", nrow(object@real.occs), '/', nrow(object@real.bg), "\n")
            cat(" partition method: ", object@null.partition.method, "\n")
            cat(" partition settings: ", paste(names(object@null.partition.settings), unlist(object@null.partition.settings), sep = " = ", collapse = ", "), "\n")
            cat(" other settings: ", paste(names(object@null.other.settings), unlist(object@null.other.settings), sep = " = ", collapse = ", "), "\n")
            cat(" algorithm: ", object@null.algorithm, "\n")
            cat(" model settings: \n")
            print(object@null.mod.settings[,-ncol(object@null.mod.settings)], row.names = FALSE)
            cat("Refer to ?ENMnull for information on slots.", sep = "")
            invisible(NULL)
          })

