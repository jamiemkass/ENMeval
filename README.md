[![CRAN version](https://www.r-pkg.org/badges/version/ENMeval)](https://CRAN.R-project.org/package=ENMeval) [![downloads](https://cranlogs.r-pkg.org:443/badges/grand-total/ENMeval?color=orange)](https://cranlogs.r-pkg.org:443/badges/grand-total/ENMeval?color=orange) [![Build Status](https://travis-ci.com/jamiemkass/ENMeval.svg?branch=master)](https://travis-ci.com/jamiemkass/ENMeval)

# ENMeval version 1.9.0
R package for automated tuning and evaluations of ecological niche models

[`ENMeval`](https://cran.r-project.org/package=ENMeval) is an R package that performs automated tuning and evaluations of ecological niche models. Version 1.9.0 represents an extensive restructure and expansion of version 0.3.1, and has many new features, including customizable specification of algorithms besides Maxent using the new **ENMdetails** object, comprehensive metadata output, null model evaluations, new visualization tools, and more flexibility for different analyses and data types. Most of these new features were created in response to user requests -- thank you for your input!

`ENMeval` 1.9.0 includes the functionality to specify any algorithm of choice, but comes out of the box with two implementations of Maxent: maxnet [(Phillips *et al.* 2017)](https://onlinelibrary.wiley.com/doi/full/10.1111/ecog.03049) from the ['maxnet' R package](https://cran.r-project.org/package=maxnet) and maxent.jar [(Phillips *et al.* 2006)](https://doi.org/10.1016/j.ecolmodel.2005.03.026), available [here](http://biodiversityinformatics.amnh.org/open_source/maxent/). 

Model tuning refers to the process of building models with varying complexity settings, then choosing optimal settings based on some criteria. As it is difficult to predict in advance what level of complexity best fits your data and results in the most ecologically realistic response for your species, model tuning and evaluations are essential for ENM studies. This process helps researchers maximize predictive ability and avoid overfitting with models that are too complex. 

The main function  `ENMevaluate()` returns several items including a table of evaluation results, lists of models and raster predictions for each complexity setting combination, a metadata object with details about the analysis, and other related data. The package also has tools for null model simulations to determine significance and effect sizes of performance metrics, tools for calculating niche overlap between predictions, a suite of visualization tools, parallel processing functionality, a vignette with a complete analysis walk-through, and more. 

For a more detailed description of the original package version, check out the open-access publication:

[Muscarella, R., Galante, P. J., Soley-Guardia, M., Boria, R. A., Kass, J. M., Uriarte, M. and Anderson, R. P. (2014), ENMeval: An R package for conducting spatially independent evaluations and estimating optimal model complexity for Maxent ecological niche models. Methods in Ecology and Evolution, 5: 1198–1205.](http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12261/full)

* Note that as of version 0.3.0, the default implementation uses the ['maxnet' R package](https://cran.r-project.org/package=maxnet). The output from this differs from that of the original java program and so some features are not compatible (e.g., variable importance, the old html output).
