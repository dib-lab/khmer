readsFile=$1
percent=$2
sketchSize=$3
outputPrefix=$4


countMinSketchSize=$(echo 100| awk '{print $1/4}')


sketchFile=$(mktemp mmapExperiment.XXX.kh)
sketchCQFFile=$(mktemp mmapExperiment.XXX.kh)

distFile=$(mktemp /tmp/mmapExperiment.XXX.dist)
distFileMMap=$(mktemp /tmp/mmapExperiment.XXX.dist)


noTotalLines=$(wc -l $readsFile|cut -f1 -d ' ' )
subset=$(echo $noTotalLines $percent | awk '{print ($1/4) * $2 * 4}')



subsetFile=$readsFile

head -n $subset < $readsFile > $subsetFile 

echo "Counting Normal"

python3 scripts/load-into-counting.py   -x $sketchSize $sketchCQFFile $subsetFile


echo "Counting CQF"
python3 scripts/load-into-counting.py  -c  -x $sketchSize $sketchFile $subsetFile


echo "Normal"

/usr/bin/time -v python3 scripts/abundance-dist.py   -s $sketchFile $subsetFile $distFile  2> $outputPrefix.res.$percent.$sketchSize


echo "CQF"


/usr/bin/time -v python3 scripts/abundance-dist.py  -m -s $sketchCQFFile $subsetFile $distFileMMap  2> $outputPrefix.res.mmap.$percent.$sketchSize  2> $outputPrefix.res.mmap.$percent.$sketchSize



diff $distFile $distFileMMap >  $outputPrefix.diff.$percent.$sketchSize

rm  -f $subsetFile $sketchFile $distFile $distFileMMap

