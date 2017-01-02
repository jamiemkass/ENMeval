pntContr <- function(m) {
  res <- m@results
  pc <- res[grepl('contribution', rownames(res)),]
  pc <- sort(pc, decreasing = TRUE)
  varnames <- sapply(strsplit(names(pc), '.contribution'), function(x) x[1])
  pc.df <- data.frame(variable=varnames, pntContr=pc, row.names=NULL)
  pc.df$variable <- factor(pc.df$variable, levels = varnames)
  return(pc.df)
}