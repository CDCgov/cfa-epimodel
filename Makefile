ifndef ENGINE
ENGINE := podman
endif

ifndef IMAGE
IMAGE := cfa-epimodel
endif

help:
	@echo "Makefile for the $(IMAGE) project."
	@echo ""
	@echo "Usage: make [target]"
	@echo ""
	@echo "Available targets:"
	@echo "  help              : Displays this message."
	@echo "  build             : Builds the R package under epimodelcfa."
	@echo "  install           : After building the package, installs the package."
	@echo "  clean             : Deletes epimodelcfa_*tar*gz files."
	@echo "  docs              : Runs devtools::document() within the package."
	@echo "  check             : Runs R CMD check on the package."
	@echo "  bic               : Builds, installs, and cleans the package."
	@echo ""

build: clean
	R CMD build .

install:
	R CMD INSTALL --preclean --clean -l epimodelcfa_*.tar.gz --dependencies=TRUE

clean:
	rm -f epimodelcfa_*tar.gz

docs:
	Rscript -e 'devtools::document()'

check:
	R CMD check . --no-manual

bic: build install clean

.PHONY: help install clean docs check build bic
