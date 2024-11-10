#' @title Parse Maxent lambdas information
#'
#' @description NOTICE: This function was borrowed from the rmaxent package 
#' written by John Baumgartner (https://github.com/johnbaums/rmaxent/).
#  It is included here with John's permission to make ENMeval CRAN-compatible 
#' (dependencies on Github-only packages are not allowed for CRAN).
#' 
#' Parse Maxent .lambdas files to extract the types, weights, minima and maxima 
#' of features, as well as the fitted model's entropy and other values required
#' for predicting to new data.
#'
#' @param lambdas Either a `MaxEnt` fitted model object (fitted with the 
#'   `maxent` function in the `dismo` package), or a file path to a 
#'   Maxent .lambdas file.
#' @return A list (of class `lambdas`) with five elements: 
#' * `lambdas`: a `data.frame` describing the features used in
#'   a Maxent model, including their weights (lambdas), maxima, minima, and 
#'   type;
#' * `linearPredictorNormalizer`: a constant that ensures the
#'   linear predictor (the sum of clamped features multiplied by their 
#'   respective feature weights) is always negative (for numerical stability);
#' * `densityNormalizer`: a scaling constant that ensures Maxent's 
#'   raw output sums to 1 over background points;
#' * `numBackgroundPoints`: the number of background points used in
#'   model training; and
#' * `entropy`: the entropy of the fitted model.
#' @references 
#' * Wilson, P. W. (2009) [_Guidelines for computing MaxEnt model output values from a lambdas file_](http://gis.humboldt.edu/OLM/Courses/GSP_570/Learning\%20Modules/10\%20BlueSpray_Maxent_Uncertinaty/MaxEnt\%20lambda\%20files.pdf).
#' * _Maxent software for species habitat modeling, version 3.3.3k_ help file (software freely available [here](https://www.cs.princeton.edu/~schapire/maxent/)).
#' @importFrom methods is
#' @importFrom utils count.fields
#' @importFrom stats setNames
#' @export
#' @examples
#' # Below we use the predicts::MaxEnt example to fit a model:
#' library(predicts)
#' occs <- read.csv(file.path(system.file(package="predicts"),
#' "/ex/bradypus.csv"))[,2:3]
#' predictors <- rast(file.path(system.file(package='predicts'), '/ex/bio.tif'))
#' me <- MaxEnt(predictors, occs)
#' # ... and then parse the lambdas information:
#' lam <- parse_lambdas(me)
#' lam
#' str(lam, 1)

parse_lambdas <- function(lambdas) {
  if(methods::is(lambdas, 'MaxEnt_model')) {
    lambdas <- lambdas@lambdas
  } else {
    lambdas <- readLines(lambdas)
  }
  con <- textConnection(lambdas)
  n <- utils::count.fields(con, ',', quote='')
  close(con)
  meta <- stats::setNames(lapply(strsplit(lambdas[n==2], ', '), 
                                 function(x) as.numeric(x[2])),
                          sapply(strsplit(lambdas[n==2], ', '), '[[', 1))
  lambdas <- stats::setNames(data.frame(do.call(
    rbind, strsplit(lambdas[n==4], ', ')), stringsAsFactors=FALSE),
    c('feature', 'lambda', 'min', 'max'))
  lambdas[, -1] <- lapply(lambdas[, -1], as.numeric)
  lambdas$feature <- sub('=', '==', lambdas$feature)
  lambdas$feature <- sub('<', '<=', lambdas$feature)
  lambdas$type <- factor(sapply(lambdas$feature, function(x) {
    switch(gsub("\\w|\\.|-|\\(|\\)", "", x),
           "==" = 'categorical',
           "<=" = "threshold",
           "^" = "quadratic",
           "*" = "product", 
           "`" = "reverse_hinge",
           "'" = 'forward_hinge',
           'linear')
  }))
  vars <- gsub("\\^2|\\(.*<=|\\((.*)==.*|`|\\'|\\)", "\\1", lambdas$feature)
  lambdas$var <- sub('\\*', ',', vars)
  l <- c(list(lambdas=lambdas[, c(1, 6, 2:5)]), meta)
  class(l) <- 'lambdas'
  l
}
