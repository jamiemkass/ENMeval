# Calculate AICc from Maxent model prediction

This function calculates AICc for Maxent models based on Warren and
Seifert (2011).

## Usage

``` r
aic.maxent(p.occs, ncoefs, p = NULL)
```

## Arguments

- p.occs:

  data frame: raw (maxent.jar) or exponential (maxnet) predictions for
  the occurrence localities based on one or more models

- ncoefs:

  numeric: number of non-zero model coefficients

- p:

  SpatRaster: raw (maxent.jar) or exponential (maxnet) model
  predictions; if NULL, AICc will be calculated based on the background
  points, which already have predictions that sum to 1 and thus need no
  correction – this assumes that the background points represent a good
  sample of the study extent

## Value

data frame with three columns: `AICc` is the Akaike Information
Criterion corrected for small sample sizes calculated as: \$\$ (2 \* K -
2 \* logLikelihood) + (2 \* K) \* (K+1) / (n - K - 1)\$\$ where *K* is
the number of non-zero coefficients in the model and *n* is the number
of occurrence localities. The *logLikelihood* is calculated as: \$\$
sum(log(vals))\$\$ where *vals* is a vector of Maxent raw/exponential
values at occurrence localities and the sum of these values across the
study extent is equal to 1. `delta.AICc` is the difference between the
AICc of a given model and the AICc of the model with the lowest AICc.
`w.AICc` is the Akaike weight (calculated as the relative likelihood of
a model (exp(-0.5 \* `delta.AICc`)) divided by the sum of the likelihood
values of all models included in a run. These can be used for model
averaging (Burnham and Anderson 2002).

## Details

As motivated by Warren and Seifert (2011) and implemented in ENMTools
(Warren *et al.* 2010), this function calculates the small sample size
version of Akaike Information Criterion for ENMs (Akaike 1974). We use
AICc (instead of AIC) regardless of sample size based on the
recommendation of Burnham and Anderson (1998, 2004). The number of
coefficients is determined by counting the number of non-zero
coefficients in the `maxent` lambda file (`m@lambdas` for maxent.jar and
`m$betas` for maxnet. See Warren *et al.* (2014) for limitations of this
approach, namely that the number of non-zero coefficients is an estimate
of the true degrees of freedom. For Maxent ENMs, AICc is calculated by
first standardizing the raw predictions such that all cells in the study
extent sum to 1, then extracting the occurrence record predictions. The
predictions of the study extent may not sum to 1 if the background does
not cover every grid cell – as the background predictions sum to 1 by
definition, extra predictions for grid cells not in the training data
will add to this sum. When no raster data is provided, the raw
predictions of the occurrence records are used to calculate AICc without
standardization, with the assumption that the background records have
adequately represented the occurrence records. The standardization is
not necessary here because the background predictions sum to 1 already,
and the occurrence data is a subset of the background. This will not be
true if the background does not adequately represent the occurrence
records, in which case the occurrences are not a subset of the
background and the raster approach should be used instead. The
likelihood of the data for a given model is then calculated by taking
the product of the raw occurrence predictions (Warren and Seifert 2011),
or the sum of their logs, as is implemented here.

## Note

Returns all `NA`s if the number of non-zero coefficients is larger than
the number of observations (occurrence localities).

## References

Akaike, H. (1974) A new look at the statistical model identification.
*IEEE Transactions on Automatic Control*, **19**: 716-723.
[doi:10.1109/TAC.1974.1100705](https://doi.org/10.1109/TAC.1974.1100705)

Burnham, K. P. and Anderson, D. R. (1998) Model selection and multimodel
inference: a practical information-theoretic approach. Springer, New
York.

Burnham, K. P. and Anderson, D. R. (2004) Multimodel inference:
understanding AIC and BIC in model selection. *Sociological Methods and
Research*, **33**: 261-304.
[doi:10.1177/0049124104268644](https://doi.org/10.1177/0049124104268644)

Warren, D. L., Glor, R. E, and Turelli, M. (2010) ENMTools: a toolbox
for comparative studies of environmental niche models. *Ecography*,
**33**: 607-611.
[doi:10.1111/j.1600-0587.2009.06142.x](https://doi.org/10.1111/j.1600-0587.2009.06142.x)

Warren, D. L., & Seifert, S. N. (2011). Ecological niche modeling in
Maxent: the importance of model complexity and the performance of model
selection criteria. *Ecological Applications*, **21**: 335-342.
[doi:10.1890/10-1171.1](https://doi.org/10.1890/10-1171.1)

Warren, D. L., Wright, A. N., Seifert, S. N., and Shaffer, H. B. (2014).
Incorporating model complexity and sampling bias into ecological niche
models of climate change risks faced by 90 California vertebrate species
of concern. *Diversity and Distributions*, **20**: 334-343.
[doi:10.1111/ddi.12160](https://doi.org/10.1111/ddi.12160)

## See also

`MaxEnt` for the predicts package.
