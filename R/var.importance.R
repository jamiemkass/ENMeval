var.importance <- function(mod) {
  res <- mod@results
  pc <- res[grepl('contribution', rownames(res)),]
  pi <- res[grepl('permutation', rownames(res)),]
  varnames <- sapply(strsplit(names(pc), '.contribution'), function(x) x[1])
  df <- data.frame(variable=varnames, percent.contribution=pc, permutation.importance=pi, row.names=NULL)
  return(df)
}
