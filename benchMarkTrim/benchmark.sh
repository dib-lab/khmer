reads=$1
refIndex=$2
outPrefix=$3

bowtie2 -p 4 -x $refIndex -U $reads |samtools view -bS - > $outPrefix.bam

./sam-scan.py $refIndex.fa <(samtools view $outPrefix.bam) -o $outPrefix.bam.pos

./summarize-pos-file.py $outPrefix.bam.pos $reads  > $outPrefix.report
