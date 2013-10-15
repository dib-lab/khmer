all:
	python setup.py build

install:
	python setup.py install

clean:
	python setup.py clean --all
	cd lib && make clean
	cd tests && rm -rf khmertest_*

doc: FORCE
	python setup.py build_sphinx --fresh-env

lib:
	cd lib && \
	$(MAKE)

test: all
	python setup.py nosetests

FORCE:
