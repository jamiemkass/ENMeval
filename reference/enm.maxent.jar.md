# ENMdetails maxent.jar

This is the ENMdetails implementation for maxent.jar, the Java version
of the Maxent algorithm. The configuration for running the model differs
slightly from that in previous versions of ENMeval (0.3.0 and before) in
that this version (\>=2.0.0) uses the default of adding presences to the
background for model training, while previous versions had turned this
off. Specifically, previous versions ran maxent() with
"noaddsamplestobackground" in the "args" vector argument, while this
version does not.

## Usage

``` r
enm.maxent.jar
```
