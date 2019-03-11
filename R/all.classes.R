#' @export

# class slots match older ENMeval versions
ENMevaluation <- setClass("ENMevaluation",
                          slots=c(algorithm='character',
                                  tuned.settings = 'data.frame',
                                  results='data.frame',
                                  results.grp='data.frame',
                                  predictions='RasterStack',
                                  models='list',
                                  partition.method='character',
                                  occ.pts='data.frame',
                                  occ.grp='numeric',
                                  bg.pts='data.frame',
                                  bg.grp='numeric',
                                  overlap='list'))
#' @export
setMethod("show",
		  signature="ENMevaluation",
		  definition=function(object) {
		  	cat("An object of class: ", class(object), "\n", sep="")
		  	cat(" ", "Occurrence/background points: ",
		  		paste(nrow(object@occ.pts), '/', nrow(object@bg), sep=''), '\n',
		  		" ",  "Partition method: ", object@partition.method, '\n',
		  		" ",  "@algorithm             : character of algorithm used", '\n',
		  		" ",  "@tuned.settings        : data.frame of settings that were tuned", '\n',
		  		" ",  "@results               : data.frame of evaluation summary statistics", '\n',
		  		" ",  "@results.grp           : data.frame of evaluation k-fold statistics", '\n',
		  		" ",  "@predictions           : RasterStack of model predictions", '\n',
          " ",  "@models                : list of model objects", '\n',
		  		" ",  "@partition.method      : character of partition method used", '\n',
		  		" ",  "@occ.pts               : data.frame of occurrence coordinates", '\n',
		  		" ",  "@occ.grp               : vector of bins for occurrence points", '\n',
		  		" ",  "@bg.pts                : data.frame of background coordinates", '\n',
		  		" ",  "@bg.grp                : vector of bins for background points", '\n',
		  		" ",  if (length(object@overlap) > 0){"@overlap               : list of matrices of pairwise niche overlap statistic"}, '\n', sep="")
		  	invisible(NULL)
		  })
