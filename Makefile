# Parts of the Makefile were adapted from the Khmer project one.
# Copyright (c) 2010-2014, Michigan State University. All rights reserved.

# project dependencies
DEPS=Biopython

# pattern to find Python source files
SOURCES=$(wildcard *.py)

# default parameters used in Python scripts
output=output
threshold=1.0

# checks for required variables
guard-%:
	@ if [ "${${*}}" = "" ]; then \
		echo "Required variable '$*' not set"; \
		exit 1; \
	fi

## help:			: print the help message
help: Makefile
	@sed -n 's/^##//p' $<

## install-deps		: install the dependencies to run dn develop the project
install-deps: install-dependencies

install-dependencies: 
	pip2 install --upgrade $(DEPS) || pip install --upgrade $(DEPS)

## clean			: clean temporary build files
clean:
	rm filter/*.pyc || true
	./setup.py clean --all || true
	rm coverage-debug || true
	rm -Rf .coverage || true
	rm -Rf $(output) || true

## test			: run the test suite
test:
	./setup.py nosetests

## pep8			: check Python code style
pep8: $(SOURCES)
	pep8 --exclude=_version.py  --show-source --show-pep8 setup.py filter/ \
		parsers/ ./ || true

pep8_report.txt: $(SOURCES)
	pep8 --exclude=_version.py setup.py filter/ parsers/ ./ \
		> pep8_report.txt || true

diff_pep8_report: pep8_report.txt
	diff-quality --violations=pep8 pep8_report.txt

## pep257      		: check Python code style
pep257: $(SOURCES)
	#pep257 --ignore=D100,D101,D102,D103 \
	pep257 \
		setup.py filter/ parsers/ || true

pep257_report.txt: $(SOURCES) $(wildcard tests/*.py)
	pep257 setup.py filter/ parsers/ \
		> pep257_report.txt 2>&1 || true

diff_pep257_report: pep257_report.txt
	diff-quality --violations=pep8 pep257_report.txt

## --- running the pipeline

## filter_by_expect	: run the filterByExpect.py script
filter_by_expect: guard-infile
	mkdir -p $(output)
	./filterByExpect.py --infile=$(infile) --output=$(output) --threshold=$(threshold)
