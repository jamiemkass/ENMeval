# Load ENMevaluation object

Load an ENMevaluation object as an .rds file. This is necessary to use
instead of
[`readRDS()`](https://rspatial.github.io/terra/reference/serialize.html)
because wrapped terra SpatRasters require
[`unwrap()`](https://rspatial.github.io/terra/reference/wrap.html) after
loading for the raster data. This convenience function does that for
you.

## Usage

``` r
loadENMevaluation(filename)
```

## Arguments

- filename:

  character: path to the .rds file to load
