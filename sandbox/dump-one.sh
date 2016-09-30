#!/usr/bin/env bash
set -eo pipefail

n=$1
refr=$2

seqid="chr$n"
workdir=PotentiallyInterestingReads/$seqid
mkdir -p $workdir
for sample in NA19238 NA19239 NA19240
do
    outprefix=${workdir}/${sample}.${seqid}.nb
    dump-boring-reads.py --seqid=$seqid $refr ${sample}.bam \
        2> ${outprefix}.log \
        | gzip -c > ${outprefix}.fq.gz &
done
wait
