# Build metadata object from ENMeval results

Builds a `rangeModelMetadata` object from the output of `ENMevaluate`.
See Merow *et al.* (2019) for more details on the nature of the metadata
and the `rangeModelMetadata` package. To improve reproducibility of the
study, this metadata object can be used as supplemental information for
a manuscript, shared with collaborators, etc.

## Usage

``` r
buildRMM(e, envs, rmm = NULL)
```

## Arguments

- e:

  ENMevaluation object

- envs:

  SpatRaster: environmental predictor variables used in analysis; needed
  to pull information on the predictor variables not included in the
  ENMevaluation object

- rmm:

  rangeModelMetadata object: if included, fields are appended to this
  RMM object as opposed to returning a new RMM object

## References

Merow, C., Maitner, B. S., Owens, H. L., Kass, J. M., Enquist, B. J.,
Jetz, W., & Guralnick, R. (2019). Species' range model metadata
standards: RMMS. *Global Ecology and Biogeography*, **28**: 1912-1924.
[doi:10.1111/geb.12993](https://doi.org/10.1111/geb.12993)
