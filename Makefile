all:
	mkdir -p prototypes && cp -a inst/prototypes/*.ipynb prototypes
	cp docs/index.html docs/.index.html
	./release
	mv docs/.index.html docs/index.html
	rm -r prototypes

