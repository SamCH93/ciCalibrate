## ----------------------------------------------------------------------------
## Document/Build/Check/Install R package (depends on devtools and roxygen2)
## Author: Samuel Pawel
## adapted from Manuela Ott, Sebastian Meyer, Florian Gerber
## ----------------------------------------------------------------------------

PACKAGE = ciCalibrate
VERSION = 0.42
TAR = $(PACKAGE)_$(VERSION).tar.gz

all: build

description:
	sed -i -r -- 's/^Version:.*/Version: '$(VERSION)'/g' DESCRIPTION ;      
	sed -i -r -- 's/^Date:.*/Date: '`date +'%F'`'/g' DESCRIPTION ;
	R -e 'roxygen2::roxygenize()'

document: description
	R -e 'devtools::document()'

manual: document
	R -e 'devtools::build_manual()'

$(TAR): manual
	R -e 'devtools::build()'

build: $(TAR)

install: $(TAR)
	R -e 'devtools::install(build = FALSE)'

check:
	R -e 'devtools::check(cran = FALSE)

cran:
	R -e 'devtools::check(cran = TRUE, remote = TRUE)'

test:
	R -e 'devtools::test()'

.PHONY: all document manual build install check cran description test
