default:
	Rscript -e "install.packages('./', type='source', repos=NULL)" && Rscript -e "library(mmbr); testthat::test_examples('.')"

document:
	Rscript -e "devtools::document()"
test:
	Rscript -e "devtools::test()"
html:
	mkdir -p prototypes && cp -a inst/prototypes/*.ipynb prototypes
	cp docs/index.html docs/.index.html
	./release
	mv docs/.index.html docs/index.html
	rm -r prototypes
pkg_down:
	R --slave -e "pkgdown::build_site(lazy=T,examples=F)"
