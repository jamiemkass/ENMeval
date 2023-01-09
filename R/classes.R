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
#' @slot doClamp logical: whether or not clamping was used 
#' @slot clamp.directions list: the clamping directions specified 
#' @slot results data frame: evaluation summary statistics
#' @slot results.partitions data frame: evaluation k-fold statistics
#' @slot models list: model objects
#' @slot variable.importance list: variable importance data frames (when available)
#' @slot predictions RasterStack: model predictions
#' @slot taxon.name character: the name of the focal taxon (optional)
#' @slot occs data frame: occurrence coordinates and predictor variable values used for model training
#' @slot occs.testing data frame: when provided, the coordinates of the fully-withheld testing records
#' @slot occs.grp vector: partition groups for occurrence points
#' @slot bg data frame: background coordinates and predictor variable values used for model training
#' @slot bg.grp vector: partition groups for background points
#' @slot overlap list: matrices of pairwise niche overlap statistics
#' @slot rmm list: the rangeModelMetadata objects for each model
#' 
#' @details The following are brief descriptions of the columns in the results table, which prints
#' when accessing `e@results` or `results(e)` if `e` is the ENMevaluation object. Those columns
#' that represent evaluations of validation data (__.val.__) end in either "avg" (average of the
#' metric across the models trained on withheld data during cross-validation) or "sd" (standard
#' deviation of the metric across these models).\cr*
#' fc = feature class\cr*
#' rm = regularization multiplier\cr*
#' tune.args = combination of arguments that define the complexity settings used for tuning (i.e., fc and rm for Maxent)\cr*
#' auc.train = AUC calculated on the full dataset\cr*
#' cbi.train = Continuous Boyce Index calculated on the full dataset\cr*
#' auc.val = average/sd AUC calculated on the validation datasets (the data withheld during cross-validation)\cr*
#' auc.diff = average/sd difference between auc.train and auc.val\cr*
#' or.mtp = average/sd omission rate with threshold as the minimum suitability value across occurrence records\cr*
#' or.10p = average/sd omission rate with threshold as the minimum suitability value across occurrence records after removing the lowest 10%\cr*
#' cbi.val = average/sd Continuous Boyce Index calculated on the validation datasets (the data withheld during cross-validation)\cr*
#' AICc = AIC corrected for small sample sizes\cr*
#' delta.AICc = highest AICc value across all models minus this model's AICc value, where lower values mean higher performance and 0 is the highest performing model\cr*
#' w.AIC = AIC weights, calculated by exp( -0.5 * delta.AIC), where higher values mean higher performance\cr*
#' ncoef = number of non-zero beta values (model coefficients)
#' 
#' @references 
#' For references on performance metrics, see the following:
#' 
#' In general for ENMeval: 
#' 
#' Muscarella, R., Galante, P. J., Soley-Guardia, M., Boria, R. A., Kass, J. M., Uriarte, M., & Anderson, R. P. (2014). ENMeval: An R package for conducting spatially independent evaluations and estimating optimal model complexity for Maxent ecological niche models. \emph{Methods in Ecology and Evolution}, \bold{5}: 1198-1205. \doi{10.1111/2041-210X.12261}
#' 
#' \emph{AUC}
#' 
#' Fielding, A. H., & Bell, J. F. (1997). A review of methods for the assessment of prediction errors in conservation presence/absence models. \emph{Environmental Conservation}, \bold{24}: 38-49. \doi{10.1017/S0376892997000088}
#' 
#' Jiménez‐Valverde, A. (2012). Insights into the area under the receiver operating characteristic curve (AUC) as a discrimination measure in species distribution modelling. \emph{Global Ecology and Biogeography}, \bold{21}: 498-507. \doi{10.1111/j.1466-8238.2011.00683.x}
#' 
#' \emph{AUC diff}
#' 
#' Warren, D. L., Glor, R. E., Turelli, M. & Funk, D. (2008) Environmental niche equivalency versus conservatism: quantitative approaches to niche evolution. \emph{Evolution}, \bold{62}: 2868-2883. \doi{10.1111/j.1558-5646.2008.00482.x}
#' 
#' Radosavljevic, A., & Anderson, R. P. (2014). Making better Maxent models of species distributions: complexity, overfitting and evaluation. \emph{Journal of Biogeography}, \bold{41}(4), 629-643. \doi{10.1111/jbi.12227} 
#' 
#' \emph{Omission rates}
#' 
#' Radosavljevic, A., & Anderson, R. P. (2014). Making better Maxent models of species distributions: complexity, overfitting and evaluation. \emph{Journal of Biogeography}, \bold{41}(4), 629-643. \doi{10.1111/jbi.12227}
#' 
#' \emph{Continuous Boyce Index}
#' 
#' Hirzel, A. H., Le Lay, G., Helfer, V., Randin, C., & Guisan, A. (2006). Evaluating the ability of habitat suitability models to predict species presences. \emph{Ecological Modelling}, \bold{199}: 142-152. \doi{10.1016/j.ecolmodel.2006.05.017}
#' 
#' @rdname ENMevaluation
#' @export ENMevaluation

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

#' @title eval.models generic for ENMevaluation object
#' @param x ENMevaluation object
#' @rdname eval.models
#' @export
setGeneric("eval.models", function(x) standardGeneric("eval.models"))

#' @rdname eval.models
setMethod("eval.models", "ENMevaluation", function(x) x@models)

#' @title eval.variable.importance (variable importance) generic for ENMevaluation object
#' @param x ENMevaluation object
#' @rdname eval.variable.importance
#' @export
setGeneric("eval.variable.importance", function(x) standardGeneric("eval.variable.importance"))

#' @rdname eval.models
setMethod("eval.variable.importance", "ENMevaluation", function(x) x@variable.importance)

#' @title eval.predictions generic for ENMevaluation object
#' @param x ENMevaluation object
#' @rdname eval.predictions
#' @export
setGeneric("eval.predictions", function(x) standardGeneric("eval.predictions"))

#' @rdname eval.predictions
setMethod("eval.predictions", "ENMevaluation", function(x) x@predictions)

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

#' @param object ENMevaluation object
#' @rdname ENMevaluation
setMethod("show",
          signature = "ENMevaluation",
          definition = function(object) {
            cat("An object of class: ", class(object), "\n")
            if(nchar(object@taxon.name)>0) cat(" taxon name: ", object@taxon.name, "\n")
            cat(" occurrence/background points: ", nrow(object@occs), '/', nrow(object@bg), "\n")
            cat(" partition method: ", object@partition.method, "\n")
            cat(" partition settings: ", ifelse(length(object@partition.settings) > 0, paste(names(object@partition.settings), unlist(object@partition.settings), sep = " = ", collapse = ", "), "none"), "\n")
            clamp.dir.spacing <- "\n         "
            if(object@doClamp == FALSE) cat(" clamp: ", object@doClamp, "\n")
            if(object@doClamp == TRUE) cat(" clamp: ", paste(sapply(1:2, function(x) paste0(names(object@clamp.directions[x]), ": ", paste(object@clamp.directions[[x]], collapse = ", "))), collapse = clamp.dir.spacing), "\n")
            cat(" categoricals: ", paste(object@other.settings$categoricals, collapse = ", "), "\n")
            cat(" algorithm: ", object@algorithm, "\n")
            tune.args.df <- object@tune.settings[,-which(names(object@tune.settings)=="tune.args"), drop = FALSE]
            tune.spacing <- "\n                 "
            cat(" tune settings: ", paste0(names(tune.args.df), ": ", apply(tune.args.df, 2, function(x) paste(unique(x), collapse = ",")), collapse = tune.spacing), "\n")
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
#' @slot errors function: returns errors chosen by the user to prevent any malfunction in the analysis.
#' The available arguments are: occs, envs, bg, tune.args, partitions, algorithm, partition.settings, other.settings, 
#' categoricals, doClamp, clamp.directions.
#' @slot msgs function: prints messages showing the package version number, etc., and those related to the input tuning parameters \code{tune.args}.
#' The available arguments are: tune.args, other.settings.
#' @slot args function: returns the parameters needed to run the algorithm function.
#' The available arguments are: occs.z, bg.z, tune.tbl.i, other.settings (where x.z is a data.frame of the envs values at
#' coordinates of x, and tune.tbl.i is a single set of tuning parameters).
#' @slot predict function: specifies how to calculate a model prediction for a Raster* or a data frame.
#' The available arguments are: mod, envs, other.settings.
#' @slot ncoefs function: counts the number of non-zero model coefficients.
#' The available arguments are: mod.
#' @slot variable.importance function: generates a data frame of variable importance from the model object (if functionality is available).
#' The available arguments are: mod.
#' @rdname ENMdetails
#' @export ENMdetails

ENMdetails <- setClass("ENMdetails",
                       slots = c(name = 'character',
                                 fun = 'function',
                                 errors = 'function',
                                 msgs = 'function',
                                 args = 'function',
                                 predict = 'function',
                                 ncoefs = 'function',
                                 variable.importance = 'function'))

#' @title eval.name generic for ENMdetails object
#' @param x ENMdetails object
#' @param value input value
#' @rdname enm.name
#' @export
setGeneric("enm.name", function(x) standardGeneric("enm.name"))

#' @rdname enm.name
#' @export
setGeneric("enm.name<-", function(x, value) standardGeneric("enm.name<-"))

#' @rdname enm.name
setMethod("enm.name", "ENMdetails", function(x) x@name)

#' @rdname enm.name
setMethod("enm.name<-", "ENMdetails", function(x, value) {
  x@name <- value
  validObject(x)
  x
})

#' @title enm.fun generic for ENMdetails object
#' @param x ENMdetails object
#' @param value input value
#' @rdname enm.fun
#' @export
setGeneric("enm.fun", function(x) standardGeneric("enm.fun"))

#' @rdname enm.fun
#' @export
setGeneric("enm.fun<-", function(x, value) standardGeneric("enm.fun<-"))

#' @rdname enm.fun
setMethod("enm.fun", "ENMdetails", function(x) x@fun)

#' @rdname enm.fun
setMethod("enm.fun<-", "ENMdetails", function(x, value) {
  x@fun <- value
  validObject(x)
  x
})

#' @title enm.errors generic for ENMdetails object
#' @param x ENMdetails object
#' @param value input value
#' @rdname enm.errors
#' @export
setGeneric("enm.errors", function(x) standardGeneric("enm.errors"))

#' @rdname enm.errors
#' @export
setGeneric("enm.errors<-", function(x, value) standardGeneric("enm.errors<-"))

#' @rdname enm.errors
setMethod("enm.errors", "ENMdetails", function(x) x@errors)

#' @rdname enm.errors
setMethod("enm.errors<-", "ENMdetails", function(x, value) {
  x@errors <- value
  validObject(x)
  x
})

#' @title enm.msgs generic for ENMdetails object
#' @param x ENMdetails object
#' @param value input value
#' @rdname enm.msgs
#' @export
setGeneric("enm.msgs", function(x) standardGeneric("enm.msgs"))

#' @rdname enm.msgs
#' @export
setGeneric("enm.msgs<-", function(x, value) standardGeneric("enm.msgs<-"))

#' @rdname enm.msgs
setMethod("enm.msgs", "ENMdetails", function(x) x@msgs)

#' @rdname enm.msgs
setMethod("enm.msgs<-", "ENMdetails", function(x, value) {
  x@msgs <- value
  validObject(x)
  x
})

#' @title enm.args generic for ENMdetails object
#' @param x ENMdetails object
#' @param value input value
#' @rdname enm.args
#' @export
setGeneric("enm.args", function(x) standardGeneric("enm.args"))

#' @rdname enm.args
#' @export
setGeneric("enm.args<-", function(x, value) standardGeneric("enm.args<-"))

#' @rdname enm.args
setMethod("enm.args", "ENMdetails", function(x) x@args)

#' @rdname enm.args
setMethod("enm.args<-", "ENMdetails", function(x, value) {
  x@args <- value
  validObject(x)
  x
})

#' @title enm.predict generic for ENMdetails object
#' @param x ENMdetails object
#' @param value input value
#' @rdname enm.predict
#' @export
setGeneric("enm.predict", function(x) standardGeneric("enm.predict"))

#' @rdname enm.predict
#' @export
setGeneric("enm.predict<-", function(x, value) standardGeneric("enm.predict<-"))

#' @rdname enm.predict
setMethod("enm.predict", "ENMdetails", function(x) x@predict)

#' @rdname enm.predict
setMethod("enm.predict<-", "ENMdetails", function(x, value) {
  x@predict <- value
  validObject(x)
  x
})

#' @title enm.ncoefs generic for ENMdetails object
#' @param x ENMdetails object
#' @param value input value
#' @rdname enm.ncoefs
#' @export
setGeneric("enm.ncoefs", function(x) standardGeneric("enm.ncoefs"))

#' @rdname enm.ncoefs
#' @export
setGeneric("enm.ncoefs<-", function(x, value) standardGeneric("enm.ncoefs<-"))

#' @rdname enm.ncoefs
setMethod("enm.ncoefs", "ENMdetails", function(x) x@ncoefs)

#' @rdname enm.ncoefs
setMethod("enm.ncoefs<-", "ENMdetails", function(x, value) {
  x@ncoefs <- value
  validObject(x)
  x
})

#' @title enm.variable.importance generic for ENMdetails object
#' @param x ENMdetails object
#' @param value input value
#' @rdname enm.variable.importance
#' @export
setGeneric("enm.variable.importance", function(x) standardGeneric("enm.variable.importance"))

#' @rdname enm.variable.importance
#' @export
setGeneric("enm.variable.importance<-", function(x, value) standardGeneric("enm.variable.importance<-"))

#' @rdname enm.variable.importance
setMethod("enm.variable.importance", "ENMdetails", function(x) x@variable.importance)

#' @rdname enm.variable.importance
setMethod("enm.variable.importance<-", "ENMdetails", function(x, value) {
  x@variable.importance <- value
  validObject(x)
  x
})

#' @param object ENMdetails object
#' @rdname ENMdetails
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
#' @slot null.doClamp logical: whether to clamp model predictions or not
#' @slot null.other.settings list: other modeling settings used (i.e., decisions about clamping, AUC diff calculation)
#' @slot null.no.iter numeric: number of null model iterations
#' @slot null.results data frame: evaluation summary statistics for null models
#' @slot null.results.partitions data frame: evaluation k-fold statistics for null models
#' @slot null.emp.results data frame: evaluation summary statistics for the empirical model, means for all null models, z-scores, and p-values
#' @slot emp.occs data frame: occurrence coordinates and predictor variable values used for model training (empirical model)
#' @slot emp.occs.grp vector: partition groups for occurrence points (empirical model)
#' @slot emp.bg data frame: background coordinates and predictor variable values used for model training (empirical model)
#' @slot emp.bg.grp vector: partition groups for background points (empirical model)
#' @rdname ENMnull
#' @export ENMnull

# class slots match older ENMeval versions
ENMnull <- setClass("ENMnull",
                    slots = c(null.algorithm = 'character',
                              null.mod.settings = 'data.frame',
                              null.partition.method = 'character',
                              null.partition.settings = 'list',
                              null.doClamp = 'logical',
                              null.other.settings = 'list',
                              null.no.iter = 'numeric',
                              null.results = 'data.frame',
                              null.results.partitions = 'data.frame',
                              null.emp.results = 'data.frame',
                              emp.occs = 'data.frame',
                              emp.occs.grp = 'factor',
                              emp.bg = 'data.frame',
                              emp.bg.grp = 'factor'))

#' @title null.algorithm generic for ENMnull object
#' @param x ENMnull object
#' @rdname null.algorithm
#' @export
setGeneric("null.algorithm", function(x) standardGeneric("null.algorithm"))

#' @rdname null.algorithm
setMethod("null.algorithm", "ENMnull", function(x) x@null.algorithm)

#' @title null.mod.settings generic for ENMnull object
#' @param x ENMnull object
#' @rdname null.mod.settings
#' @export
setGeneric("null.mod.settings", function(x) standardGeneric("null.mod.settings"))

#' @rdname null.mod.settings
setMethod("null.mod.settings", "ENMnull", function(x) x@null.mod.settings)

#' @title null.partition.method generic for ENMnull object
#' @param x ENMnull object
#' @rdname null.partition.method
#' @export
setGeneric("null.partition.method", function(x) standardGeneric("null.partition.method"))

#' @rdname null.partition.method
setMethod("null.partition.method", "ENMnull", function(x) x@null.partition.method)

#' @title null.partition.settings generic for ENMnull object
#' @param x ENMnull object
#' @rdname null.partition.settings
#' @export
setGeneric("null.partition.settings", function(x) standardGeneric("null.partition.settings"))

#' @rdname null.partition.settings
setMethod("null.partition.settings", "ENMnull", function(x) x@null.partition.settings)

#' @title null.doClamp generic for ENMnull object
#' @param x ENMnull object
#' @rdname null.doClamp
#' @export
setGeneric("null.doClamp", function(x) standardGeneric("null.doClamp"))

#' @rdname null.doClamp
setMethod("null.doClamp", "ENMnull", function(x) x@null.doClamp)

#' @title null.other.settings generic for ENMnull object
#' @param x ENMnull object
#' @rdname null.other.settings
#' @export
setGeneric("null.other.settings", function(x) standardGeneric("null.other.settings"))

#' @rdname null.other.settings
setMethod("null.other.settings", "ENMnull", function(x) x@null.other.settings)

#' @title null.no.iter generic for ENMnull object
#' @param x ENMnull object
#' @rdname null.no.iter
#' @export
setGeneric("null.no.iter", function(x) standardGeneric("null.no.iter"))

#' @rdname null.no.iter
setMethod("null.no.iter", "ENMnull", function(x) x@null.no.iter)

#' @title null.results generic for ENMnull object
#' @param x ENMnull object
#' @rdname null.results
#' @export
setGeneric("null.results", function(x) standardGeneric("null.results"))

#' @rdname null.results
setMethod("null.results", "ENMnull", function(x) x@null.results)

#' @title null.results.partitions generic for ENMnull object
#' @param x ENMnull object
#' @rdname null.results.partitions
#' @export
setGeneric("null.results.partitions", function(x) standardGeneric("null.results.partitions"))

#' @rdname null.results.partitions
setMethod("null.results.partitions", "ENMnull", function(x) x@null.results.partitions)

#' @title null.emp.results generic for ENMnull object
#' @param x ENMnull object
#' @rdname null.emp.results
#' @export
setGeneric("null.emp.results", function(x) standardGeneric("null.emp.results"))

#' @rdname null.emp.results
setMethod("null.emp.results", "ENMnull", function(x) x@null.emp.results)

#' @title emp.occs generic for ENMnull object
#' @param x ENMnull object
#' @rdname emp.occs
#' @export
setGeneric("emp.occs", function(x) standardGeneric("emp.occs"))

#' @rdname emp.occs
setMethod("emp.occs", "ENMnull", function(x) x@emp.occs)

#' @title emp.occs.grp generic for ENMnull object
#' @param x ENMnull object
#' @rdname emp.occs.grp
#' @export
setGeneric("emp.occs.grp", function(x) standardGeneric("emp.occs.grp"))

#' @rdname emp.occs.grp
setMethod("emp.occs.grp", "ENMnull", function(x) x@emp.occs.grp)

#' @title emp.bg generic for ENMnull object
#' @param x ENMnull object
#' @rdname emp.bg
#' @export
setGeneric("emp.bg", function(x) standardGeneric("emp.bg"))

#' @rdname emp.bg
setMethod("emp.bg", "ENMnull", function(x) x@emp.bg)

#' @title emp.bg.grp generic for ENMnull object
#' @param x ENMnull object
#' @rdname emp.bg.grp
#' @export
setGeneric("emp.bg.grp", function(x) standardGeneric("emp.bg.grp"))

#' @rdname emp.bg.grp
setMethod("emp.bg.grp", "ENMnull", function(x) x@emp.bg.grp)

#' @param object ENMnull object
#' @rdname ENMnull
setMethod("show",
          signature = "ENMnull",
          definition = function(object) {
            cat("An object of class: ", class(object), "\n")
            cat(" no. iterations: ", object@null.no.iter, "\n")
            cat(" empirical occurrence/background points: ", nrow(object@emp.occs), '/', nrow(object@emp.bg), "\n")
            cat(" partition method: ", object@null.partition.method, "\n")
            cat(" partition settings: ", paste(names(object@null.partition.settings), unlist(object@null.partition.settings), sep = " = ", collapse = ", "), "\n")
            clamp.dir.spacing <- "\n         "
            if(object@null.doClamp == FALSE) cat(" clamp: ", object@null.doClamp, "\n")
            if(object@null.doClamp == TRUE) cat(" clamp: ", paste(sapply(1:2, function(x) paste0(names(object@null.other.settings$clamp.directions[x]), ": ", paste(object@null.other.settings$clamp.directions[[x]], collapse = ", "))), collapse = clamp.dir.spacing), "\n")
            cat(" categoricals: ", paste(object@null.other.settings$categoricals, collapse = ", "), "\n")
            cat(" algorithm: ", object@null.algorithm, "\n")
            # cat(" model settings: \n")
            tune.args.df <- object@null.mod.settings[,-which(names(object@null.mod.settings)=="tune.args"), drop = FALSE]
            tune.spacing <- "\n                  "
            cat(" model settings: ", paste0(names(tune.args.df), ": ", apply(tune.args.df, 2, function(x) paste(unique(x), collapse = ",")), collapse = tune.spacing), "\n")
            # print(object@null.mod.settings[,-ncol(object@null.mod.settings)], row.names = FALSE)
            cat("Refer to ?ENMnull for information on slots.", sep = "")
            invisible(NULL)
          })

