NoKmers=1000000
K=20
NoUnseenKmers=10000
outputPrefix=data1000000

python generateSeq.py $NoKmers $K  $NoUnseenKmers $outputPrefix

python3 plotGoldHist.py $outputPrefix

parallel --gnu  -k "python3 testPerfomance.py $outputPrefix {1} {2}" :::    21 22 23 24 25  :::  --cqf --cm > $outputPrefix.result 

python3 plotPerformanceBoxPlot.py $outputPrefix
