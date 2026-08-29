# Corrected variance function

Calculate variance corrected for non-independence of *k*-fold iterations

## Usage

``` r
corrected.var(x, nk)
```

## Arguments

- x:

  numeric vector: input values

- nk:

  numeric: number of *k*-fold iterations

## Value

A numeric value of the corrected variance.

## Details

\`corrected.var\` calculates variance corrected for non-independence of
*k*-fold iterations. See Appendix of Shcheglovitova & Anderson (2013)
and other references (Miller 1974; Parr 1985; Shao and Wu 1989) for
additional details. This function calculates variance that is corrected
for the non-independence of *k* cross-validation iterations. Following
Shao and Wu (1989):

\$\$Sum Of Squares \* ((n-1)/n)\$\$

where *n* = the number of *k*-fold iterations.

## References

Miller, R. G. (1974) The jackknife - a review. *Biometrika*, **61**:
1-15. [doi:10.1093/biomet/61.1.1](https://doi.org/10.1093/biomet/61.1.1)

Parr, W. C. (1985) Jackknifing differentiable statistical functionals.
*Journal of the Royal Statistics Society, Series B*, **47**: 56-66.
[doi:10.1111/j.2517-6161.1985.tb01330.x](https://doi.org/10.1111/j.2517-6161.1985.tb01330.x)

Shao J. and Wu, C. F. J. (1989) A general theory for jackknife variance
estimation. *Annals of Statistics*, **17**: 1176-1197.
[doi:10.1214/aos/1176347263](https://doi.org/10.1214/aos/1176347263)

Shcheglovitova, M. and Anderson, R. P. (2013) Estimating optimal
complexity for ecological niche models: a jackknife approach for species
with small sample sizes. *Ecological Modelling*, **269**: 9-17.
[doi:10.1016/j.ecolmodel.2013.08.011](https://doi.org/10.1016/j.ecolmodel.2013.08.011)

## Author

Robert Muscarella \<bob.muscarella@gmail.com\>
