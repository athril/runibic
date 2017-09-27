# Makefile for generating R packages.
# Author: 2011 Andrew Redd
# Modifications: 2017 Patryk Orzechowski
# Roxygen uses the roxygen2 package, and will run automatically on check and all.


PKG_VERSION=$(shell grep -i ^Version: DESCRIPTION | sed 's/^.*: //g')
PKG_NAME=$(shell grep -i ^Package: DESCRIPTION | sed 's/^.*: //g')
CXX_STD = CXX11
PKG_CXXFLAGS = -std=c++11 -fopenmp


R_FILES := $(wildcard R/*.R)
SRC_FILES := $(wildcard src/*) $(addprefix src/, $(COPY_SRC))
PKG_FILES := DESCRIPTION NAMESPACE $(R_FILES) $(SRC_FILES)

.PHONY: tarball install check clean build list


all:build install

tarball: $(PKG_NAME)_$(PKG_VERSION).tar.gz
$(PKG_NAME)_$(PKG_VERSION).tar.gz: $(PKG_FILES)
	R CMD build .

check: $(PKG_NAME)_$(PKG_VERSION).tar.gz
	R CMD check $(PKG_NAME)_$(PKG_VERSION).tar.gz

build: $(PKG_NAME)_$(PKG_VERSION).tar.gz
	R CMD INSTALL --build $(PKG_NAME)_$(PKG_VERSION).tar.gz

install: $(PKG_NAME)_$(PKG_VERSION).tar.gz
	R CMD INSTALL $(PKG_NAME)_$(PKG_VERSION).tar.gz

NAMESPACE: $(R_FILES)
	Rscript -e "library(roxygen2);roxygenize('.')"

clean:
	-rm -f $(PKG_NAME)_*.tar.gz
	-rm -rf $(PKG_NAME).Rcheck
	-rm -rf  man
	-rm -rf  NAMESPACE
	-rm -f src/*.o
	-rm -f src/*.so

list:
	@echo "R files:"
	@echo $(R_FILES)
	@echo "Source files:"
	@echo $(SRC_FILES)