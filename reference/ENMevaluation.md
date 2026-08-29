# ENMevaluation class

An S4 class that contains the ENMevaluate results.

## Usage

``` r
# S4 method for class 'ENMevaluation'
show(object)
```

## Arguments

- object:

  ENMevaluation object

## Details

The following are brief descriptions of the columns in the results
table, which prints when accessing \`e@results\` or \`results(e)\` if
\`e\` is the ENMevaluation object. Those columns that represent
evaluations of validation data (\_\_.val.\_\_) end in either "avg"
(average of the metric across the models trained on withheld data during
cross-validation) or "sd" (standard deviation of the metric across these
models).  
\* fc = feature class  
\* rm = regularization multiplier  
\* tune.args = combination of arguments that define the complexity
settings used for tuning (i.e., fc and rm for Maxent)  
\* auc.train = AUC calculated on the full dataset  
\* cbi.train = Continuous Boyce Index calculated on the full dataset  
\* auc.val = average/sd AUC calculated on the validation datasets (the
data withheld during cross-validation)  
\* auc.diff = average/sd difference between auc.train and auc.val  
\* or.mtp = average/sd omission rate with threshold as the minimum
suitability value across occurrence records  
\* or.10p = average/sd omission rate with threshold as the minimum
suitability value across occurrence records after removing the lowest 10
cbi.val = average/sd Continuous Boyce Index calculated on the validation
datasets (the data withheld during cross-validation)  
\* AICc = AIC corrected for small sample sizes  
\* delta.AICc = highest AICc value across all models minus this model's
AICc value, where lower values mean higher performance and 0 is the
highest performing model  
\* w.AIC = AIC weights, calculated by exp( -0.5 \* delta.AIC), where
higher values mean higher performance  
\* ncoef = number of non-zero beta values (model coefficients)

## Slots

- `algorithm`:

  character: algorithm used

- `tune.settings`:

  data frame: settings that were tuned

- `partition.method`:

  character: partition method used

- `partition.settings`:

  list: partition settings used (i.e., value of \*k\* or aggregation
  factor)

- `other.settings`:

  list: other modeling settings used (i.e., decisions about clamping,
  AUC diff calculation)

- `doClamp`:

  logical: whether or not clamping was used

- `clamp.directions`:

  list: the clamping directions specified

- `results`:

  data frame: evaluation summary statistics

- `results.partitions`:

  data frame: evaluation k-fold statistics

- `models`:

  list: model objects

- `variable.importance`:

  list: variable importance data frames (when available)

- `predictions`:

  SpatRaster: model predictions

- `taxon.name`:

  character: the name of the focal taxon (optional)

- `occs`:

  data frame: occurrence coordinates and predictor variable values used
  for model training

- `occs.testing`:

  data frame: when provided, the coordinates of the fully-withheld
  testing records

- `occs.grp`:

  vector: partition groups for occurrence points

- `bg`:

  data frame: background coordinates and predictor variable values used
  for model training

- `bg.grp`:

  vector: partition groups for background points

- `overlap`:

  list: matrices of pairwise niche overlap statistics

- `rmm`:

  list: the rangeModelMetadata objects for each model

## References

For references on performance metrics, see the following:

In general for ENMeval:

Muscarella, R., Galante, P. J., Soley-Guardia, M., Boria, R. A., Kass,
J. M., Uriarte, M., & Anderson, R. P. (2014). ENMeval: An R package for
conducting spatially independent evaluations and estimating optimal
model complexity for Maxent ecological niche models. *Methods in Ecology
and Evolution*, **5**: 1198-1205.
[doi:10.1111/2041-210X.12261](https://doi.org/10.1111/2041-210X.12261)

*AUC*

Fielding, A. H., & Bell, J. F. (1997). A review of methods for the
assessment of prediction errors in conservation presence/absence models.
*Environmental Conservation*, **24**: 38-49.
[doi:10.1017/S0376892997000088](https://doi.org/10.1017/S0376892997000088)

Jimenez-Valverde, A. (2012). Insights into the area under the receiver
operating characteristic curve (AUC) as a discrimination measure in
species distribution modelling. *Global Ecology and Biogeography*,
**21**: 498-507.
[doi:10.1111/j.1466-8238.2011.00683.x](https://doi.org/10.1111/j.1466-8238.2011.00683.x)

*AUC diff*

Warren, D. L., Glor, R. E., Turelli, M. & Funk, D. (2008) Environmental
niche equivalency versus conservatism: quantitative approaches to niche
evolution. *Evolution*, **62**: 2868-2883.
[doi:10.1111/j.1558-5646.2008.00482.x](https://doi.org/10.1111/j.1558-5646.2008.00482.x)

Radosavljevic, A., & Anderson, R. P. (2014). Making better Maxent models
of species distributions: complexity, overfitting and evaluation.
*Journal of Biogeography*, **41**(4), 629-643.
[doi:10.1111/jbi.12227](https://doi.org/10.1111/jbi.12227)

*Omission rates*

Radosavljevic, A., & Anderson, R. P. (2014). Making better Maxent models
of species distributions: complexity, overfitting and evaluation.
*Journal of Biogeography*, **41**(4), 629-643.
[doi:10.1111/jbi.12227](https://doi.org/10.1111/jbi.12227)

*Continuous Boyce Index*

Hirzel, A. H., Le Lay, G., Helfer, V., Randin, C., & Guisan, A. (2006).
Evaluating the ability of habitat suitability models to predict species
presences. *Ecological Modelling*, **199**: 142-152.
[doi:10.1016/j.ecolmodel.2006.05.017](https://doi.org/10.1016/j.ecolmodel.2006.05.017)

## Author

Jamie M. Kass, <jamie.m.kass@gmail.com>, Bob Muscarella,
<bob.muscarella@gmail.com>
