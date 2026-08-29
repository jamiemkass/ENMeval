# ENMnull class

An S4 class that contains the ENMnulls results.

## Usage

``` r
# S4 method for class 'ENMnull'
show(object)
```

## Arguments

- object:

  ENMnull object

## Slots

- `null.algorithm`:

  character: algorithm used

- `null.mod.settings`:

  data frame: model settings used

- `null.partition.method`:

  character: partition method used

- `null.partition.settings`:

  list: partition settings used (i.e., value of \*k\* or aggregation
  factor)

- `null.doClamp`:

  logical: whether to clamp model predictions or not

- `null.other.settings`:

  list: other modeling settings used (i.e., decisions about clamping,
  AUC diff calculation)

- `null.no.iter`:

  numeric: number of null model iterations

- `null.results`:

  data frame: evaluation summary statistics for null models

- `null.results.partitions`:

  data frame: evaluation k-fold statistics for null models

- `null.emp.results`:

  data frame: evaluation summary statistics for the empirical model,
  means for all null models, z-scores, and p-values

- `emp.occs`:

  data frame: occurrence coordinates and predictor variable values used
  for model training (empirical model)

- `emp.occs.grp`:

  vector: partition groups for occurrence points (empirical model)

- `emp.bg`:

  data frame: background coordinates and predictor variable values used
  for model training (empirical model)

- `emp.bg.grp`:

  vector: partition groups for background points (empirical model)

## Author

Jamie M. Kass, <jamie.m.kass@gmail.com>, Corentin Bohl,
<corentinbohl@gmail.com>
