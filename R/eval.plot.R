#####################################################
#########	 PLOT OPTIMIZATION CRITERIA #############
#####################################################

#' @export

eval.plot <- function(results, value="delta.AICc", variance=NULL, legend=TRUE, legend.position="topright") {
	res <- results
	if(!class(res)=='data.frame') {
			stop("check input - eval.plot requires results table as data.frame", call.=F)
		}
	if(sum(value %in% colnames(res))==0) {
			stop("value not in results table provided", call.=F)
		}
	fc <- length(unique(res$features))
	col <- rainbow(fc)
	rm <- length(unique(res$rm))
	xlab <- "Regularization Multiplier"

	y <- res[,value]
	ylim <- c(min(y, na.rm=TRUE), max(y, na.rm=TRUE))

	if(!is.null(variance)){
		if(!variance %in% colnames(res)){
			stop("Check name of variance column")
			}
		v <- res[,variance]
		ylim <- c(min(y-v), max(y+v))
		}

	par(mar=c(4.5,4.5,1,1))
	plot(res$rm, y, col='white', ylim=ylim, ylab=value, xlab=xlab, axes=F, cex.lab=1.5)
	if(value=="delta.AICc") abline(h=2, lty=3)
	axis(1, at= unique(res$rm))
	axis(2)
	box()
	for (j in 1:length(unique(res$features))){
		s <- ((fc*rm)-fc+j)
		points(res$rm[seq(j, s, fc)], y[seq(j, s, fc)], type="l", col=col[j])
		if(!is.null(variance)){
			arrows(res$rm[seq(j, s, fc)], 
				y[seq(j, s, fc)] + v[seq(j, s, fc)], 
				res$rm[seq(j, s, fc)], 
				y[seq(j, s, fc)] - v[seq(j, s, fc)],
				code=3, length=.05, angle=90, col=col[j])
		}
	}
	points(res$rm, y, bg=col, pch=21)

	if(legend==TRUE){
	legend(legend.position, legend=unique(res$features), pt.bg=col, pch=21, bg='white')
	}
}
