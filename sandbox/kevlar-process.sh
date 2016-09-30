#!/usr/bin/env bash
set -eo pipefail

REFR=/share/hormozdiarilab/Data/ReferenceGenomes/GRCh38_full_analysis_set_plus_decoy_hla.fa

parallel --gnu --jobs 8 --ungroup dump-one.sh {} $REFR ::: {1..22} X Y
count-one.sh 16
parallel --gnu --jobs 4 --ungroup novel-one.sh {} ::: {1..22} X Y
