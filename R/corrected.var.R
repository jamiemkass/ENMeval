########################################################
######### CALCULATE CORRECTED VARIANCE  ################
########################################################

corrected.var <- function(x, nk){
	rowSums((x - rowMeans(x))^2) * ((nk-1)/nk)
}
