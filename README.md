# mvsusieR

[![Continuous Integration](https://github.com/stephenslab/mvsusieR/actions/workflows/ci.yml/badge.svg)](https://github.com/stephenslab/mvsusieR/actions/workflows/ci.yml)

Implements a multivariate generalization of the "Sum of Single
Effects" (SuSiE) model for variable selection in multivariate linear
regression (Y = XB + E, where Y is N x R). The methods are
particularly well-suited to settings where some of the X variables
are highly correlated, the true effects are sparse, and effects may
be shared across multiple conditions (traits). A major motivation is
multi-trait genetic fine-mapping, but the methods are applicable more
generally.

Go [here][pkgdown-site] for documentation and tutorials.

## Setup

To install the latest version of the mvsusieR package from GitHub:

```R
# install.packages("remotes")
remotes::install_github("stephenslab/mvsusieR")
```

This command should automatically install all required packages if
they are not installed already.

Note that mvsusieR depends on [mashr][mashr], which requires the GNU
Scientific Library (GSL). You may want to install mashr first before
attempting to install mvsusieR.

## Citing this work

If you find the `mvsusieR` package useful for your work, please cite:

> Zou, Y., Carbonetto, P., Xie, D., Wang, G. & Stephens, M. (2026).
> Fast and flexible joint fine-mapping of multiple traits via the Sum
> of Single Effects model. *Nature Genetics* **58**, 454-462.
> https://doi.org/10.1038/s41588-025-02486-7

## Developer notes

+ When any changes are made to `roxygen2` markup, run
`devtools::document()` to update package `NAMESPACE`
and documentation files.

+ To install and test the package:

    ```bash
    VERSION=`grep Version DESCRIPTION | awk '{print $2}'`
    R CMD build --resave-data --no-build-vignettes mvsusieR
    R CMD INSTALL mvsusieR_$VERSION.tar.gz
    R CMD check --as-cran --ignore-vignettes mvsusieR_$VERSION.tar.gz
    ```

+ Run `pkgdown::build_site()` to build the pkgdown site.

[pkgdown-site]: https://stephenslab.github.io/mvsusieR/
[mashr]: https://github.com/stephenslab/mashr
