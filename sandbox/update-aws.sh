#! /bin/bash
cd /root/screed
git pull origin master
python setup.py install

cd /root/khmer
make clean
git pull origin refactor
make clean test
