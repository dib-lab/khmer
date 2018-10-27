#! /usr/bin/env python3
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2015, The Regents of the University of California.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#
#     * Neither the name of the Michigan State University nor the names
#       of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written
#       permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Contact: khmer-project@idyll.org

from nameparser import HumanName
import codecs
import textwrap


authors = []

with codecs.open('authors.csv', 'r', encoding='utf-8') as namefile:
    for line in namefile:
        name, address = line.split(',')
        authors.append((HumanName(name), address))

authors = sorted(authors, key=lambda author: author[0].last)
authors.insert(0, (HumanName("Daniel Standage"), "daniel.standage@gmail.com"))
authors.append((HumanName("C. Titus Brown"), "titus@idyll.org"))

# print(authors)

bibtex = '   author = \"'

for tup in authors:
    name = tup[0]
    name.string_format = "{last}, {first} {middle} and"
    bibtex += str(name) + " "

bibtex = bibtex[:-5] + '"'  # remove last 'and' and close the quote

print('  @article{khmer2015,')
for line in textwrap.wrap(str(bibtex), 77):
    print('  ' + line)

print(
    '''     title = "The khmer software package: enabling efficient nucleotide
  sequence analysis",
     year = "2015",
     month = "08",
     publisher = "F1000",
     url = "https://doi.org/10.12688/f1000research.6924.1"
  }''')

doclist = u':Authors: '

for tup in authors:
    name = tup[0]
    name.string_format = "{first} {middle} {last}"
    doclist += str(name) + ", "

doclist = doclist[:-2]

for line in textwrap.wrap(str(doclist), 71):
    print('        ' + line)
