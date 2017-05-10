# linpk

An R package for generating concentration-time profiles from linear pharmacokinetic (PK) systems.

To install directly from GitHub:

``` r
require(devtools)
devtools::install_github("benjaminrich/linpk")
```

To try a demo shiny app:

``` r
# Make sure the required packages are installed
require(shiny)
require(shinyjs)
require(shinyAce)
require(dygraphs)
require(linpk)  # See above for installation
shiny::runGitHub("linpk", "benjaminrich", subdir="inst/demo-app")
```

