# mvsusieR

[![Appveyor Build status](https://ci.appveyor.com/api/projects/status/fhrp1e868f40skp1?svg=true](https://ci.appveyor.com/project/pcarbo/mvsusieR)]
[![CircleCI](https://circleci.com/gh/stephenslab/mvsusieR/tree/master.svg?style=svg)](https://app.circleci.com/pipelines/github/stephenslab/mvsusieR?branch=master)

## Setup

To automatically retrieve and install the latest version of the R
package from this repository, run

```r
remotes::install_github("stephenslab/mvsusieR")
```

## Quick Start

[Here](https://stephenslab.github.io/mvsusieR/articles/prediction.html) is
a quick document to show mvsusieR in action.  For more documentation and
examples please visit: https://stephenslab.github.io/mvsusieR

## Developer notes

+ When any changes are made to `roxygen2` markup, simply run 
`devtools::document()` to update package `NAMESPACE`
and documentation files.

+ To install and test the package, run the following commands
in the shell:

    ```bash
    VERSION=`grep Version DESCRIPTION | awk '{print $2}'`
    R CMD build --resave-data --no-build-vignettes mvsusieR
    R CMD INSTALL mvsusieR_$VERSION.tar.gz
    R CMD check --as-cran --ignore-vignettes mvsusieR_$VERSION.tar.gz
    ```

+ Run `pkgdown::build_site(lazy = TRUE,examples = FALSE)` to build the
  website.
