############################################################################
#########	 CALCULATE SCHOENERS D-STATISTIC (NICHE OVERLAP)  #############
############################################################################

calc.niche.overlap <- function(predictive.maps, stat="D", maxent.args){
	overlap <- matrix(nrow=nlayers(predictive.maps), ncol=nlayers(predictive.maps))
	pb <- txtProgressBar(0, nlayers(predictive.maps)-1, style=3)
	for (i in 1:(nlayers(predictive.maps)-1)){
		setTxtProgressBar(pb, i)
		for (j in (i+1):nlayers(predictive.maps)){
			overlap[j, i] <- nicheOverlap(predictive.maps[[i]], predictive.maps[[j]], stat=stat)
		}
	}
colnames(overlap) <- names(predictive.maps)
rownames(overlap) <- names(predictive.maps)
return(overlap)
}
