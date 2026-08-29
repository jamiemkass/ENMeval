# Iterate tuning of ENMs

Internal functions to tune and summarize results for ecological niche
models (ENMs) iteratively across a range of user-specified tuning
settings. See
[ENMevaluate](https://jamiemkass.github.io/ENMeval/reference/ENMevaluate.md)
for descriptions of shared arguments. Function `tune.parallel()` tunes
ENMs with parallelization. Function `cv.enm()` calculates training and
validation evaluation statistics for one set of specified tuning
parameters.

Validation CBI is calculated here with background values, not raster
data, in order to standardize the methodology for both training and
validation data for spatial partitions, as ENMeval does not mask rasters
to partition areas and hence does not have partitioned raster data.
Further, predictions for occurrence and background localities are
combined as input for the parameter "fit" in `ecospat::ecospat_boyce()`
because the interval is determined from "fit" only, and if test
occurrences all have higher predictions than the background, the
interval will be cut short.

## Usage

``` r
tune.train(
  enm,
  occs.z,
  bg.z,
  mod.full,
  tune.tbl.i,
  other.settings,
  partitions,
  quiet
)

tune.validate(
  enm,
  occs.train.z,
  occs.val.z,
  bg.train.z,
  bg.val.z,
  mod.k,
  nk,
  tune.tbl.i,
  other.settings,
  partitions,
  user.eval,
  quiet
)

tune(
  d,
  enm,
  partitions,
  tune.tbl,
  doClamp,
  other.settings,
  partition.settings,
  user.val.grps,
  occs.testing.z,
  numCores,
  parallel,
  user.eval,
  algorithm,
  updateProgress,
  quiet
)

cv.enm(
  d,
  enm,
  partitions,
  tune.tbl.i,
  doClamp,
  other.settings,
  partition.settings,
  user.val.grps,
  occs.testing.z,
  user.eval,
  algorithm,
  quiet
)
```

## Arguments

- enm:

  [ENMdetails](https://jamiemkass.github.io/ENMeval/reference/ENMdetails.md)
  object

- occs.z:

  data.frame: the envs values for the coordinates at the full dataset
  occurrence records

- bg.z:

  data.frame: the envs values for the coordinates at the full dataset
  background records

- mod.full:

  model object: the model trained on the full dataset

- tune.tbl.i:

  vector: single set of tuning parameters

- other.settings:

  named list: used to specify extra settings for the analysis. All of
  these settings have internal defaults, so if they are not specified
  the analysis will be run with default settings. See Details for
  descriptions of these settings, including how to specify arguments for
  maxent.jar.

- partitions:

  character: name of partitioning technique (see
  [`?partitions`](https://jamiemkass.github.io/ENMeval/reference/partitions.md))

- quiet:

  boolean: if TRUE, silence all function messages (but not errors).

- occs.train.z:

  data.frame: the envs values for the coordinates at the training
  occurrence records

- occs.val.z:

  data.frame: the envs values for the coordinates at the validation
  occurrence records

- bg.train.z:

  data.frame: the envs values for the coordinates at the training
  background records

- bg.val.z:

  data.frame: the envs values for the coordinates at the validation
  background records

- mod.k:

  model object: the model trained on the training dataset that becomes
  evaluated on the validation data

- nk:

  numeric: the number of folds (i.e., partitions) – will be equal to
  `kfolds` for random partitions

- user.eval:

  function: custom function for specifying performance metrics not
  included in ENMeval. The function must first be defined and then input
  as the argument `user.eval`. This function should have a single
  argument called `vars`, which is a list that includes different data
  that can be used to calculate the metric. See Details below and the
  vignette for a worked example.

- d:

  data frame: data frame from
  [`ENMevaluate()`](https://jamiemkass.github.io/ENMeval/reference/ENMevaluate.md)
  with occurrence and background coordinates (or coordinates plus
  predictor variable values) and partition group values

- tune.tbl:

  data frame: all combinations of tuning parameters

- doClamp:

  boolean: if TRUE (default), model prediction extrapolations will be
  restricted to the upper and lower bounds of the predictor variables.
  Clamping avoids extreme predictions for environment values outside the
  range of the training data. If free extrapolation is a study aim, this
  should be set to FALSE, but for most applications leaving this at the
  default of TRUE is advisable to avoid unrealistic predictions. When
  predictor variables are input, they are clamped internally before
  making model predictions when clamping is on. When no predictor
  variables are input and data frames of coordinates and variable values
  are used instead (SWD format), validation data is clamped before
  making model predictions when clamping is on.

- partition.settings:

  named list: used to specify certain settings for partitioning schema.
  See Details and ?partitions for descriptions of these settings.

- user.val.grps:

  matrix / data frame: user-defined validation record coordinates and
  predictor variable values. This is used internally by
  [`ENMnulls()`](https://jamiemkass.github.io/ENMeval/reference/ENMnulls.md)
  to force each null model to evaluate with empirical validation data,
  and does not have any current use when running
  [`ENMevaluate()`](https://jamiemkass.github.io/ENMeval/reference/ENMevaluate.md)
  independently.

- occs.testing.z:

  data.frame: when fully withheld testing data is provided, the envs
  values for the coordinates at the testing occurrence records

- numCores:

  numeric: number of cores to use for parallel processing. If NULL, all
  available cores will be used.

- parallel:

  boolean: if TRUE, run with parallel processing.

- algorithm:

  character: name of the algorithm used to build models. Currently one
  of "maxnet", "maxent.jar", or "bioclim", else the name from a custom
  ENMdetails implementation.

- updateProgress:

  boolean: if TRUE, use shiny progress bar. This is only for use in
  shiny apps.
