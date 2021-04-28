[![CRAN version](https://www.r-pkg.org/badges/version/ENMeval)](https://CRAN.R-project.org/package=ENMeval) [![downloads](https://cranlogs.r-pkg.org:443/badges/grand-total/ENMeval?color=orange)](https://cranlogs.r-pkg.org:443/badges/grand-total/ENMeval?color=orange) [![Build Status](https://travis-ci.com/jamiemkass/ENMeval.svg?branch=master)](https://travis-ci.com/jamiemkass/ENMeval) [![DOI](https://zenodo.org/badge/29864043.svg)](https://zenodo.org/badge/latestdoi/29864043)

# ENMeval version 2.0.0
R package for automated tuning and evaluations of ecological niche models

[`ENMeval`](https://cran.r-project.org/package=ENMeval) is an R package that performs automated tuning and evaluations of ecological niche models and species distribution models. Version 2.0.0 represents an extensive restructure and expansion of version 0.3.1, and has many new features, including customizable specification of algorithms besides Maxent using the new **ENMdetails** object, comprehensive metadata output, null model evaluations, new visualization tools, a completely updated and extensive vignette, and more flexibility for different analyses and data types. Many of these new features were created in response to user requests -- thank you for your input!

`ENMeval` 2.0.0 includes the functionality to specify any algorithm of choice, but comes out of the box with two implementations of Maxent: maxnet [(Phillips *et al.* 2017)](https://onlinelibrary.wiley.com/doi/full/10.1111/ecog.03049) from the [maxnet R package](https://cran.r-project.org/package=maxnet) and the Java software maxent.jar [(Phillips *et al.* 2006)](https://doi.org/10.1016/j.ecolmodel.2005.03.026), available [here](http://biodiversityinformatics.amnh.org/open_source/maxent/), as well as BIOCLIM implemented with the [dismo R package](https://cran.r-project.org/package=dismo).

Model tuning refers to the process of building models with varying complexity settings, then choosing optimal settings based on some criteria. As it is difficult to predict in advance what level of complexity best fits your data and results in the most ecologically realistic response for your species, model tuning and evaluations are essential for ENM studies. This process helps researchers maximize predictive ability and avoid overfitting with models that are too complex. 

The main function `ENMevaluate()` returns several items including a table of evaluation results, lists of models and raster predictions for each complexity setting combination, a metadata object with details about the analysis, and other related data. The package also has tools for null model simulations to determine significance and effect sizes of performance metrics, tools for calculating niche overlap between predictions, a suite of visualization tools, parallel processing functionality, a vignette with a complete analysis walk-through, and more. 

For a more detailed description of version 2.0.0, please reference the new open-access publication in Methods in Ecology and Evolution:

Kass, J. M., Muscarella, R., Galante, P. J., Bohl, C., Pinilla-Buitrago, G. E., Boria, R. A., Soley‐Guardia, M., & Anderson, R. P. (2021). ENMeval 2.0: redesigned for customizable and reproducible modeling of species’ niches and distributions. Methods in Ecology and Evolution, Accepted.

For the original package version, please reference the new open-access publication in Methods in Ecology and Evolution:

[Muscarella, R., Galante, P. J., Soley-Guardia, M., Boria, R. A., Kass, J. M., Uriarte, M. and Anderson, R. P. (2014), ENMeval: An R package for conducting spatially independent evaluations and estimating optimal model complexity for Maxent ecological niche models. Methods in Ecology and Evolution, 5: 1198–1205.](http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12261/full)

NOTES:

1. Please make sure to use the most recent version of [maxent.jar](https://biodiversityinformatics.amnh.org/open_source/maxent/) (currently 3.4.4), as recent bug fixes were made.

2. Note that as of version 0.3.0, the default implementation uses the ['maxnet' R package](https://cran.r-project.org/package=maxnet). The output from this differs from that of the original Java program and so some features are not compatible (e.g., variable importance, html output).
