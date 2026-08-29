# Save ENMevaluation object

Save an ENMevaluation object as an .rds file. This is necessary to use
instead of
[`saveRDS()`](https://rspatial.github.io/terra/reference/serialize.html)
because terra SpatRasters require
[`wrap()`](https://rspatial.github.io/terra/reference/wrap.html) before
saving to preserve the connections to the raster data. This convenience
function does that for you.

## Usage

``` r
saveENMevaluation(e, filename)
```

## Arguments

- e:

  ENMevaluation object

- filename:

  character: path to the file to create with .rds extension
