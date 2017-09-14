# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2016, The Regents of the University of California.
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
"""Common methods for trimming short reads on k-mer abundance."""
import screed


def trim_record(countgraph, record, cutoff, variable_coverage=False,
                normalize_to=None):
    name = record.name
    seq = record.sequence
    seqN = record.cleaned_seq

    if variable_coverage:  # only trim when sequence has high enough C
        if not countgraph.median_at_least(seqN, normalize_to):
            return record, False                 # return unmodified

    _, trim_at = countgraph.trim_on_abundance(seqN, cutoff)

    # too short? eliminate read
    if trim_at < countgraph.ksize():
        return None, True

    # would we trim? if not, return unmodified.
    if trim_at == len(seq):
        return record, False

    # construct new record
    trim_seq = seq[:trim_at]
    if hasattr(record, 'quality'):
        trim_qual = record.quality[:trim_at]
        trim_rec = screed.Record(name=name, sequence=trim_seq,
                                 quality=trim_qual)
    else:
        trim_rec = screed.Record(name=name, sequence=trim_seq)

    return trim_rec, True
