NoKmers=1000000
K=20
NoUnseenKmers=10000
outputPrefix=new

###DATASet Generation
echo "Generate DataSet"
python generateSeq.py $NoKmers $K  $NoUnseenKmers $outputPrefix
python3 plotGoldHist.py $outputPrefix
cat $outputPrefix.dat |uniq > $outputPrefix.uniq.dat

##Compile cqf tests
echo ""
echo "Compile CQF tests"
make

### CQF UNIT TEST
echo ""
echo "Unit Test"
py.test test_CQF.py

### Load Factor Test
echo ""
echo "Load Factor"
echo "M\tMaximum Number of Unique Kmers" > $outputPrefix.loadfactor.res.tsv
seq 1 16| awk '{print (2^$1)-1}' | parallel --gnu -k "python3 testLoadFactorCQF.py $outputPrefix.uniq.dat  8192  {} 2>> $outputPrefix.log |tail -n1" >> $outputPrefix.loadfactor.res.tsv
seq 1 10  | parallel --gnu -k "python3 testLoadFactorCQF.py $outputPrefix.uniq.dat  8192  {}  |tail -n1" >> $outputPrefix.loadfactor1-10.res.tsv


### Accuracy Test
echo ""
echo "Accuracy Test"
##Exp1
python testSketchesAccuracy.py $outputPrefix.uniq.dat $outputPrefix.none.dat 2>> $outputPrefix.log > $outputPrefix.bloom_vs_cqf.tsv

##Exp2
parallel --gnu  -k "python3 testPerfomance.py $outputPrefix {1} {2} 2>> $outputPrefix.log" :::     23 24 25 26 27 28  :::  --cqf --cm > $outputPrefix.result
python3 plotPerformanceBoxPlot.py $outputPrefix

### Merge Test
echo ""
echo "Merge Test Same Size"
./mergeTest_SameSize
echo ""
echo "Merge Test Different Size"
./mergeTest_DifferentSize

#### Resize Test
echo ""
echo "Change reminder size"
./mergeTest_Qbits
echo ""
echo "Resize Test"
./mergeTest_Resize
