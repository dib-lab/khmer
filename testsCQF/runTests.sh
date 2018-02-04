## initialize the parameters needed to generate the test data
K=20                  ## k-mer size
NoKmers=1000000       ## number of unique kmers used for the Unique and Zipifan datasets
NoUnseenKmers=10000   ## number of unique kmers used to generate the Unseen dataset
outputPrefix="cqfTest"

###DATASet Generation
echo "Generate test dataSets"
python3 generateSeq.py $NoKmers $K  $NoUnseenKmers $outputPrefix
python3 plotGoldHist.py $outputPrefix
cat $outputPrefix.dat |uniq > $outputPrefix.uniq.dat

##Compile cqf tests
echo ""
echo "Compile CQF tests"
make


### CQF construction test
echo ""
echo "CQF Construction Test"
./cqfConstruction

### CQF UNIT TEST
echo ""
echo "CQF Unit Test"
ln -s ../tests ./
py.test test_CQF.py

### Load Factor Test
echo ""
echo "Load Factor Test"
echo "M\tMaximum Number of Unique Kmers" > $outputPrefix.loadfactor1-10.res.tsv
seq 1 10  | parallel --gnu -k "python3 testLoadFactorCQF.py $outputPrefix.uniq.dat  8192  {} 2>> LoadFactor1.log |tail -n1" >> $outputPrefix.loadfactor1-10.res.tsv
python plotLoadingFactor.py  $outputPrefix.loadfactor1-10.res.tsv normal $outputPrefix.loadfactor1-10.res.png

echo "M\tMaximum Number of Unique Kmers" > $outputPrefix.loadfactor_log1-16.res.tsv
seq 1 16| awk '{print (2^$1)-1}' | parallel --gnu -k "python3 testLoadFactorCQF.py $outputPrefix.uniq.dat  8192  {} 2>> LoadFactor2.log |tail -n1" >> $outputPrefix.loadfactor_log1-16.res.tsv
python plotLoadingFactor.py $outputPrefix.loadfactor_log1-16.res.tsv log $outputPrefix.loadfactor_log1-16.res.png

### Accuracy Test
echo ""
echo "Accuracy Test"
##Exp1
echo "Experiment 1 Quotient Filter Vs Bloom filter"
python testSketchesAccuracy.py $outputPrefix.uniq.dat $outputPrefix.none.dat 2>> $outputPrefix.log > $outputPrefix.bloom_vs_cqf.tsv
python plotQFvsBloom.py $outputPrefix.bloom_vs_cqf.tsv $outputPrefix.bloom_vs_cqf.png
##Exp2
echo "Experiment 2 (CQF Vs Count-min sketch)"
parallel --gnu  -k "python3 testPerfomance.py $outputPrefix {1} {2} 2>> $outputPrefix.log" :::     23 24 25 26 27 28  :::  --cqf --cm > $outputPrefix.result
python3 plotPerformanceBoxPlot.py $outputPrefix

### Merge and possible resizing
echo ""
echo "Merge Test Same Size"
./mergeTest_SameSize
echo ""
echo "Merge Test Different Size"
./mergeTest_DifferentSize

echo ""
echo "Resize Test"
./mergeTest_Resize
