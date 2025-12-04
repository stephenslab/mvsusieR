# mvsusieR

[![Continuous Integration](https://github.com/stephenslab/mvsusieR/actions/workflows/ci.yml/badge.svg)](https://github.com/stephenslab/mvsusieR/actions/workflows/ci.yml)

Implements a multivariate generalization of the "Sum of Single
Effects" (SuSiE) model for variable selection in multivariate linear
regression. Go [here][pkgdown-site] for documentation and tutorials.
  
## Setup

To install the latest version of the mvsusieR package from GitHub, run
the following code in R:

```R
install.packages("remotes")
library(remotes)
install_github("stephenslab/mvsusieR")
```

This command should automatically install all required packages if
they are not installed already.

Note that mvsusieR uses [mashr][mashr], which in turn requires the GNU
Scientific libraries (GSL). You may want to install mashr first before
attempting to install mvsusieR.

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

+ Run `pkgdown::build_site()` to build the pkgdown site.


+ To format R codes in the `R` folder,

   ```bash
   for i in `ls R/*.R`; do bash inst/misc/format_r_code.sh $i; done
   ```

[pkgdown-site]: https://stephenslab.github.io/mvsusieR/
[prediction-vignette]: https://stephenslab.github.io/mvsusieR/articles/prediction.html
[mashr]: https://github.com/stephenslab/mashr
