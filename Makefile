
DEPS=Biopython
PYSOURCES=$(wildcard *.py)

## help			: print the help message
help: Makefile
		@sed -n 's/^##//p' $<

## install-deps	: install the dependencies to run dn develop the project
install-deps: install-dependencies

install-dependencies: 
	pip2 install --upgrade $(DEPS) || pip install --upgrade $(DEPS)

