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
	@echo ''
	@echo '--> docs in build/sphinx/html <--'
	@echo ''

lib:
	cd lib && \
	$(MAKE)

test: all
	python setup.py nosetests

FORCE:
