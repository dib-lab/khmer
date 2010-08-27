#! /bin/bash
READ_FILE=data/25k.fa

export PYTHONPATH=python/
env/bin/python ./scripts/do-th-subset-save.py $READ_FILE
