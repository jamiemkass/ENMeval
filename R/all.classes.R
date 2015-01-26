ENMevaluation <- setClass("ENMevaluation",
					   slots=c(results='data.frame', 
					   predictions='RasterStack',
					   partition.method='character',
					   occ.pts='data.frame',
					   occ.grp='numeric',
					   bg.pts='data.frame',
					   bg.grp='numeric',
					   overlap='matrix'))

setMethod("show",
		  signature="ENMevaluation",
		  definition=function(object) {
		  	cat("An object of class: ", class(object), "\n", sep="")
		  	cat(" ", "Occurrence/background points: ", 
		  		paste(nrow(object@occ.pts), '/', nrow(object@bg.pts), sep=''), '\n',
		  		" ",  "Partition method: ", object@partition.method, '\n',
		  	    " ",  "Feature classes: ", paste(as.character(unique(object@results[,2])), collapse=', '), '\n',
		  		" ",  "Regularization multipliers: ", paste(unique(object@results[,3]), collapse=', '), '\n',
		  		" ",  "@results     : data.frame of evaluation results", '\n',
		  		" ",  "@predictions : RasterStack of model predictions", '\n',
		  		" ",  "@occ.pts     : data.frame of occurrence coordinates", '\n',
		  		" ",  "@occ.grp     : vector of bins for occurrence points", '\n',
		  		" ",  "@bg.pts      : data.frame of background coordinates", '\n',
		  		" ",  "@bg.grp      : vector of bins for background points", '\n',
		  		" ",  if(nrow(object@overlap) > 0){"@overlap     : matrix of pairwise niche overlap metric"}, '\n', sep="")
		  	invisible(NULL)
		  })
