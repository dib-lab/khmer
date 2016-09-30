#!/usr/bin/env bash
set -eo pipefail

numjobs=$1

for n in {1..22} X Y
do
    seqid=chr${n}
    for sample in NA19238 NA19239 NA19240
    do
        echo PotentiallyInterestingReads/${seqid}/${sample}.${seqid}.nb
    done
done | parallel --gnu --ungroup --jobs $numjobs load-into-counting.py -M 12G --ksize=31 {}.countgraph {}.fq.gz > {}.count.log 2>&1
