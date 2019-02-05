#' @export
ENMevaluation <- setClass("ENMevaluation",
                          slots=c(algorithm='character',
                                  stats='data.frame',
                                  kstats='data.frame',
                                  predictions='RasterStack',
                                  models='list',
                                  partition.method='character',
                                  occs='data.frame',
                                  occs.folds='numeric',
                                  bg='data.frame',
                                  bg.folds='numeric',
                                  overlap='matrix'))
#' @export
setMethod("show",
		  signature="ENMevaluation",
		  definition=function(object) {
		  	cat("An object of class: ", class(object), "\n", sep="")
		  	cat(" ", "Occurrence/background points: ",
		  		paste(nrow(object@occs), '/', nrow(object@bg), sep=''), '\n',
		  		" ",  "Partition method: ", object@partition.method, '\n',
		  		" ",  "@algorithm             : character of algorithm used", '\n',
		  		" ",  "@stats                 : data.frame of evaluation summary statistics", '\n',
		  		" ",  "@kstats                : data.frame of evaluation k-fold statistics", '\n',
		  		" ",  "@predictions           : RasterStack of model predictions", '\n',
          " ",  "@models                : list of model objects", '\n',
		  		" ",  "@partition.method      : character of partition method used", '\n',
		  		" ",  "@occs                  : data.frame of occurrence coordinates", '\n',
		  		" ",  "@occs.folds            : vector of bins for occurrence points", '\n',
		  		" ",  "@bg                    : data.frame of background coordinates", '\n',
		  		" ",  "@bg.folds              : vector of bins for background points", '\n',
		  		" ",  if (nrow(object@overlap) > 0){"@overlap     : matrix of pairwise niche overlap metric"}, '\n', sep="")
		  	invisible(NULL)
		  })
