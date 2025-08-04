ifndef ENGINE
ENGINE := podman
endif

ifndef IMAGE
IMAGE := cfa-epimodel
endif

lib_path := /home/RUU7/r_libs

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
	@echo ""

build: clean
	R CMD build .

install:
	R CMD INSTALL --preclean --clean -l $(lib_path) epimodelcfa_*.tar.gz --dependencies=TRUE

clean:
	rm -f epimodelcfa_*tar.gz

docs:
	Rscript -e 'devtools::document()'

check:
	R CMD check .

.PHONY: container_build container_run r_package_build r_package_install docs \
	install
