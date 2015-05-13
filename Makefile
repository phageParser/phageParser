# The MIT License (MIT)
#
# Copyright (c) 2014 Sidhartha Goyal, Shirish Goyal
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

DEPS=Biopython
SOURCES=$(wildcard *.py)

## help		: print the help message
help: Makefile
		@sed -n 's/^##//p' $<

## install-deps	: install the dependencies to run dn develop the project
install-deps: install-dependencies

install-dependencies: 
	pip2 install --upgrade $(DEPS) || pip install --upgrade $(DEPS)

## clean		: clean temporary build files
clean:
	rm filter/*.pyc || true
	./setup.py clean --all || true
	rm coverage-debug || true
	rm -Rf .coverage || true

## test 		: run the test suite
test:
	./setup.py nosetests

## pep8        : check Python code style
pep8: $(SOURCES)
	pep8 --exclude=_version.py  --show-source --show-pep8 setup.py filter/ \
		parsers/ ./ || true

pep8_report.txt: $(SOURCES)
	pep8 --exclude=_version.py setup.py filter/ parsers/ ./ \
		> pep8_report.txt || true

diff_pep8_report: pep8_report.txt
	diff-quality --violations=pep8 pep8_report.txt

## pep257      : check Python code style
pep257: $(SOURCES)
	#pep257 --ignore=D100,D101,D102,D103 \
	pep257 \
		setup.py filter/ parsers/ || true

pep257_report.txt: $(SOURCES) $(wildcard tests/*.py)
	pep257 setup.py filter/ parsers/ \
		> pep257_report.txt 2>&1 || true

diff_pep257_report: pep257_report.txt
	diff-quality --violations=pep8 pep257_report.txt

