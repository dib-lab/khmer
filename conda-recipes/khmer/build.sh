#!/bin/sh

# When building from the zip file, this line should be commented out.
#export SRC_DIR=~/src/khmer
#cd $SRC_DIR

$PYTHON setup.py install

# No idea why the following is necessary... scripts/ don't appear to turn up
# anywhere other than the $PREFIX/bin stubs that setuptools installs.
for _f in $(ls $SRC_DIR/scripts/*.py); do
    cp -v $_f $PREFIX/bin
done
