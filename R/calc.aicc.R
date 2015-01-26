#############################################################
######### CALCULATE AICc  ###################################
#############################################################

calc.aicc <- function(nparam, occ, predictive.maps){
	AIC.valid <- nparam < nrow(occ)
	vals <- extract(predictive.maps, occ)
	probsum <- cellStats(predictive.maps, sum)
# The log-likelihood was incorrectly calculated (see next line) in ENMeval v.1.0.0 when working with >1 model at once.
# 	LL <- colSums(log(vals/probsum), na.rm=T)
# The corrected calculation (since v.0.1.1) is:
	LL <- colSums(log(t(t(vals)/probsum)), na.rm=T)
	AICc <- (2*nparam - 2*LL) + (2*(nparam)*(nparam+1)/(nrow(occ)-nparam-1))
	AICc[AIC.valid==FALSE] <- NA
	AICc[is.infinite(AICc)] <- NA
	if(sum(is.na(AICc))==length(AICc)){
		warning("AICc not valid - returning NA's.")
		res <- data.frame(cbind(AICc, delta.AICc=NA, w.AIC=NA, nparam=NA))
	} else {
		delta.AICc <- (AICc - min(AICc, na.rm=TRUE))
		w.AIC <- (exp(-0.5*delta.AICc))/(sum(exp(-0.5*delta.AICc), na.rm=TRUE))
		res <- data.frame(AICc, delta.AICc, w.AIC, nparam)
		rownames(res) <- NULL
	}
    rownames(res) <- NULL
    return(res)
}
