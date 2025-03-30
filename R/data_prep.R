# make a list of categorical variable levels
catLevs <- function(occs.z, envs, categoricals) {
  # make categorical levels list
  cat.levs <- list()
  for(i in 1:length(categoricals)) {
    if(!is.null(envs)) {
      cat.levs[[i]] <- terra::levels(envs[[categoricals[i]]])[[1]][,2]
    }else{
      cat.levs[[i]] <- levels(occs.z[, categoricals[i]])
    }  
  }
  return(cat.levs)
}

# # converts categorical variable fields to factors and coerces them to numbers for maxent.jar
# numFactors <- function(x, d, i) {
#   x[, categoricals[i]] <- factor(as.numeric(x[, categoricals[i]]), levels = levels(d[, categoricals[i]]))
#   return(x)
# }
# 
# # converts categorical variable fields to factors
# regFactors <- function(x, d, i) {
#   x[, categoricals[i]] <- factor(x[, categoricals[i]], levels = levels(d[, categoricals[i]]))
#   return(x)
# }

  

  
  