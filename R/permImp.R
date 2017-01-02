permImp <- function(m) {
  res <- m@results
  pi <- res[grepl('permutation', rownames(res)),]
  pi <- sort(pi, decreasing = TRUE)
  varnames <- sapply(strsplit(names(pi), '.permutation'), function(x) x[1])
  pi.df <- data.frame(variable=varnames, permImp=pi, row.names=NULL)
  pi.df$variable <- factor(pi.df$variable, levels = varnames)
  return(pi.df)
}

