#! /bin/bash
rm -fr subsets merge1 groups *.part *.paired
mkdir subsets
mkdir merge1
mkdir groups

cd subsets
python ../../scripts/do-th-subset-save.py ../part-test.fa
cd ..
python ../scripts/do-subset-merge.py 1 subsets merge1

python ../scripts/do-th-subset-load.py part-test.fa merge1/*.pmap

python ../scripts/combine-pe.py part-test.fa.part

cd groups
python ../../scripts/extract-partitions.py ../*.part.paired part-g part-test.dist

cd ..

