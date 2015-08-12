#!/usr/bin/env python

from __future__ import print_function
from __future__ import unicode_literals

from nameparser import HumanName
import codecs
import textwrap

authors = []

with codecs.open('authors.csv', 'r', encoding='utf-8') as namefile:
    for line in namefile:
        name, address = line.split(',')
        authors.append((HumanName(name), address))

authors = sorted(authors, key=lambda author: author[0].last)
authors.insert(0, (HumanName("Michael R. Crusoe"), "crusoe@ucdavis.edu"))
authors.append((HumanName("C. Titus Brown"), "titus@idyll.org"))

# print(authors)

bibtex = u'   author = \"'

for tup in authors:
    name = tup[0]
    name.string_format = "{last}, {first} {middle} and"
    bibtex += unicode(name) + " "

bibtex = bibtex[:-5] + '"'  # remove last 'and' and close the quote

print(u'  @article{khmer2015,')
for line in textwrap.wrap(unicode(bibtex), 77):
    print('  ' + line)

print(
    '''     title = "The khmer software package: enabling efficient nucleotide
  sequence analysis",
     year = "2015",
     month = "08",
     publisher = "F1000",
     url = "http://dx.doi.org/10.12688/f1000research.6924.1"
  }''')
