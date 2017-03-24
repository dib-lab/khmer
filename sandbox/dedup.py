#! /usr/bin/env python
import screed
import sys
import khmer.utils

names = set()

for record in screed.open(sys.argv[1]):
    if record.name not in names:
        khmer.utils.write_record(record, sys.stdout)
        names.add(record.name)
