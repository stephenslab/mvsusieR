default:
	Rscript -e "install.packages('./', type='source', repos=NULL)" && Rscript -e "library(mvsusieR); testthat::test_examples('.')"

document:
	Rscript -e "devtools::document()"
test:
	Rscript -e "devtools::test()"
html:
	mkdir -p prototypes && cp -a inst/prototypes/*.ipynb prototypes
	./release
	rm -r prototypes
pkg_down:
	R --slave -e "pkgdown::build_site(lazy = TRUE,examples = FALSE)"
