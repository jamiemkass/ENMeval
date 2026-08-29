# Calculate Similarity of ENMs in Geographic Space

Compute pairwise "niche overlap" in geographic space for Maxent
predictions. The value ranges from 0 (no overlap) to 1 (identical
predictions). The function uses the `nicheOverlap` function of the dismo
package (Hijmans *et al.* 2011).

## Usage

``` r
calc.niche.overlap(predictors, overlapStat, quiet = FALSE)
```

## Arguments

- predictors:

  SpatRaster: at least 2 Maxent raster predictions

- overlapStat:

  character: either "D" or "I", the statistic calculated by the
  `nicheOverlap` function of the dismo package (default: "D"), which we
  updated for terra as no correlate currently exists in the new predicts
  package

- quiet:

  boolean: if TRUE, silence all function messages (but not errors)

## Value

A matrix with the lower triangle giving values of pairwise "niche
overlap" in geographic space. Row and column names correspond to the
results table output by
[`ENMevaluate()`](https://jamiemkass.github.io/ENMeval/reference/ENMevaluate.md).

## Details

"D" refers to Schoeners *D* (Schoener 1968), while "I" refers to the *I*
similarity statistic from Warren *et al.* (2008).

## References

Hijmans, R. J., Phillips, S., Leathwick, J. & Elith, J. (2011) dismo
package for R. Available online at:
<https://cran.r-project.org/package=dismo>.

Schoener, T. W. (1968) The *Anolis* lizards of Bimini: resource
partitioning in a complex fauna. *Ecology*, **49**: 704-726.
[doi:10.2307/1935534](https://doi.org/10.2307/1935534)

Warren, D. L., Glor, R. E., Turelli, M. & Funk, D. (2008) Environmental
niche equivalency versus conservatism: quantitative approaches to niche
evolution. *Evolution*, **62**: 2868-2883.
[doi:10.1111/j.1558-5646.2008.00482.x](https://doi.org/10.1111/j.1558-5646.2008.00482.x)

## See also

\`nicheOverlap\` in the dismo package

## Author

Based on dismo::`nicheOverlap`, which is based on SDMTools::`Istat`,
updated for terra package Robert Muscarella \<bob.muscarella@gmail.com\>
