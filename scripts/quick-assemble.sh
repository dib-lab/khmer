#! /bin/bash
python ~/dev/khmer/scripts/strip-partition.py $1 | velveth $1.ass 33 -short -
velvetg $1.ass
