NoKmers=1000000
K=20
NoUnseenKmers=10000
outputPrefix=data1000000

###DATASet Generation
python generateSeq.py $NoKmers $K  $NoUnseenKmers $outputPrefix
python3 plotGoldHist.py $outputPrefix
cat $outputPrefix |uniq > $outputPrefix.uniq.dat

##Compile cqf tests
make

### CQF UNIT TEST
py.test test_CQF.py

### Load Factor Test
echo "M\tMaximum Number of Unique Kmers" > $outputPrefix.loadfactor.res.tsv
seq 1 19| awk '{print (2^$1)-1}' |
 parallel --gnu -k "python3 testLoadFactorCQF.py $outputPrefix.uniq.dat  8192  {} |tail -n1"
  >> $outputPrefix.loadfactor.res.tsv


### Accuracy Test

##Exp1
python testSketchesAccuracy.py $outputPrefix.uniq.dat $outputPrefix.none.dat > $outputPrefix.bloom_vs_cqf.tsv

##Exp2
parallel --gnu  -k "python3 testPerfomance.py $outputPrefix {1} {2}" :::     23 24 25 26 27 28  :::  --cqf --cm > $outputPrefix.result
python3 plotPerformanceBoxPlot.py $outputPrefix

### Merge Test
./mergeTest_SameSize
./mergeTest_DifferentSize

#### Resize Test
./mergeTest_Qbits
./mergeTest_Resize
