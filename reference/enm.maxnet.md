# ENMdetails maxnet

This is the ENMdetails implementation for maxnet, the R version of the
Maxent algorithm. The configuration for running the model now includes
addsamplestobackground = TRUE, which explicitly adds presences to the
background for model training, though as the current version of maxnet
has this set to TRUE as default, behavior between ENMeval versions
should not differ.

## Usage

``` r
enm.maxnet
```
