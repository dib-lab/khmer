all:
	python setup.py build

clean:
	python setup.py clean --all
	cd lib && make clean
	cd tests && rm -rf khmertest_*

doc: FORCE
	python setup.py build_sphinx

lib:
	cd lib && \
	$(MAKE)

test: all
	python setup.py nosetests

FORCE:
