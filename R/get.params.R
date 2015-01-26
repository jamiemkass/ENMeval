#################################################################
######### GET NUMBER OF PARAMETERS FROM FULL MODEL  #############
#################################################################

get.params <- function(model){
	lambdas <- model@lambdas[1:(length(model@lambdas)-4)]
	countNonZeroParams <- function(x) {
	if (strsplit(x, split=", ")[[1]][2] != '0.0') 1
	}
	no.params <- sum(unlist(sapply(lambdas, countNonZeroParams)))
	return(no.params)
}
