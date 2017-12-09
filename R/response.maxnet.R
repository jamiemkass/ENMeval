response.maxnet <- function(mod, v, type) {
  nr <- 100
  mm <- mod$samplemeans
  m <- data.frame(matrix(unlist(mm), nr, length(mm), byrow=T))
  colnames(m) <- names(mm)
  min <- mod$varmin[v]
  max <- mod$varmax[v]
  m[,v] <- seq(min - 0.1*(max-min), max+0.1*(max-min), length=100)
  preds <- predict(mod, m, type=type)
  resp.tbl <- cbind(m[,v], preds)
  return(resp.tbl)
}