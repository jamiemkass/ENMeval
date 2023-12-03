.onAttach <- function(libname, pkgname) {
  packageStartupMessage("This is ", pkgname, " version ", packageVersion(pkgname), ". \nFor worked examples, please consult the vignette: <https://jamiemkass.github.io/ENMeval/articles/ENMeval-2.0-vignette.html>.")
}

.onLoad <- function(lib, pkg) {
  # Make sure the dbplyr methods are loaded
  loadNamespace("predicts")
}