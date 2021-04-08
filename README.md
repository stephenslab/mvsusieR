# mmbr

## Setup

To automatically retrieve and install `mmbr` from this repository,

```R
devtools::install_github("stephenslab/mmbr")
```

## Quick Start

[Here](https://stephenslab.github.io/mmbr/articles/prediction.html) is
a quick document to show `mmbr` in action.  For more documentation and
examples please visit: https://stephenslab.github.io/mmbr

## Developer notes

+ When any changes are made to `roxygen2` markup, simply run 
`devtools::document()` to update package `NAMESPACE`
and documentation files.

+ To install and test the package, run the following commands
in the shell:

    ```bash
    VERSION=`grep Version DESCRIPTION | awk '{print $2}'`
    R CMD build --resave-data --no-build-vignettes mmbr
    R CMD INSTALL mmbr_$VERSION.tar.gz
    R CMD check --as-cran --ignore-vignettes mmbr_$VERSION.tar.gz
    ```

+ Run `pkgdown::build_site(lazy=T,examples=F)` to build the website.
