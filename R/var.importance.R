var.importance <- function(mod) {
  if(!'MaxEnt' %in% class(mod)){
    stop('Sorry, variable importance cannot currently be calculated with Maxnet models (only Maxent)')
    } else {
    res <- mod@results
    pc <- res[grepl('contribution', rownames(res)),]
    pi <- res[grepl('permutation', rownames(res)),]
    varnames <- sapply(strsplit(names(pc), '.contribution'), function(x) x[1])
    df <- data.frame(variable=varnames, percent.contribution=pc, permutation.importance=pi, row.names=NULL)
    return(df)
  }
}
