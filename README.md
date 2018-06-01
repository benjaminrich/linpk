# linpk

[![Travis-CI Build Status](https://travis-ci.org/benjaminrich/linpk.svg?branch=master)](https://travis-ci.org/benjaminrich/linpk)

An R package for generating concentration-time profiles from linear pharmacokinetic (PK) systems.

## Installation

To install from CRAN:

``` r
install.packages("linpk")
```

To install the latest development version directly from GitHub:

``` r
require(devtools)
devtools::install_github("benjaminrich/linpk")
```


For an introduction to the package, with usage examples, see the [vignette](https://benjaminrich.github.io/linpk/vignettes/linpk-intro.html).

There is a shiny app that provides a demo of the package capabilities, and also generates code that can be placed in an R script. To run it, paste the following lines in an R console:

``` r
require(shiny)
shiny::runGitHub("linpk", "benjaminrich", subdir="inst/demo")
```

