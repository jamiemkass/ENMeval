[![CRAN version](http://www.r-pkg.org/badges/version/wallace)](https://CRAN.R-project.org/package=wallace) [![downloads](http://cranlogs.r-pkg.org/badges/grand-total/wallace?color=orange)](http://cranlogs.r-pkg.org/badges/grand-total/wallace?color=orange)

# ENMeval
R package for automated runs and evaluations of ecological niche models.

[`ENMeval`](https://cran.r-project.org/package=ENMeval) is an R package that performs automated runs and evaluations of ecological niche models, and currently implements Maxent using the either (now by default) the 'maxnet' algoritm developed by [Phillips *et al.* (2017)](https://onlinelibrary.wiley.com/doi/full/10.1111/ecog.03049) using the ['maxnet' R package](https://cran.r-project.org/package=maxnet) or [the original java program (http://biodiversityinformatics.amnh.org/open_source/maxent/). `ENMeval` was made for those who want to "tune" their models to maximize predictive ability and avoid overfitting, or in other words, optimize model complexity to balance goodness-of-fit and predictive ability. The primary function, `ENMevaluate`, does all the heavy lifting and returns several items including a table of evaluation statistics and, for each setting combination, a model object and a raster layer showing the model prediction across the study extent. There are also options for calculating niche overlap between predictions, running in parallel to speed up computation, and more. For a more detailed description of the package, check out the open-access publication:

[Muscarella, R., Galante, P. J., Soley-Guardia, M., Boria, R. A., Kass, J. M., Uriarte, M. and Anderson, R. P. (2014), ENMeval: An R package for conducting spatially independent evaluations and estimating optimal model complexity for Maxent ecological niche models. Methods in Ecology and Evolution, 5: 1198–1205.](http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12261/full)

Also see the vignette for examples of implementation.

Note that as of version 0.3.0, the default implementation uses the ['maxnet' R package](https://cran.r-project.org/package=maxnet).  The output from this differs from that of the original java program and so some features are not compatible (e.g., variable importance, the old html output).
