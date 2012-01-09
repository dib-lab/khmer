# Packages that need to be installed:
#    Python 2.6
#    numpy
#    khmer (github.com/ctb/khmer)
#    screed (github.com/ctb/screed)
#    R
#    graphviz
#    octave

all: bitsperkmer graphvisual clustsize numparts seqerror diam perc

bitsperkmer:
	python bitsperk.py > results/bitsperk.txt 

seqerror:
	python seqerror.py data/ecoli.fa data/reads.fa > results/seqerr.txt

graphvisual: 
	python gengraphs.py
	neato -Tpdf -oresults/001.pdf 0.01.dot
	neato -Tpdf -oresults/005.pdf 0.05.dot
	neato -Tpdf -oresults/010.pdf 0.1.dot
	neato -Tpdf -oresults/015.pdf 0.15.dot

clustsize:
	python clustsize.py > clustdata.txt
	Rscript clustsize.R clustdata.txt
	mv Rplots.pdf results/clustsize.pdf

diam:
	python diam.py > diamdata.txt
	Rscript diam.R diamdata.txt
	mv Rplots.pdf results/diam.pdf

numparts:
	rm -f var.txt
	rm -f novar.txt 
	./partition.sh > partout.txt
	./output.sh
	mv Rplots.pdf results/part.pdf

perc:
	perc/do.sh > percout.txt
	mv perc/S1/S1.pdf results/perc.pdf

clean:
	rm -f *.dot
	rm -rf 0.0*
	rm -rf 0.1*
	rm -f *.ht
	rm -f *.txt
	rm -f genome.fa
	rm -f results/*
	rm -f perc/S1/*.txt
	rm -f perc/scripts/*.txt
	rm -f perc/S1/ppt
	rm -f perc/scripts/ppt
	rm -f perc/src/ppt
