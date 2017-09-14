# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2011-2015, Michigan State University.
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
# pylint: disable=missing-docstring,invalid-name

from khmer.thread_utils import ThreadedSequenceProcessor, SequenceGroup
from io import StringIO
from screed.fasta import fasta_iter
from screed.fastq import fastq_iter
import screed

# stdlib queue module was renamed on Python 3
try:
    import queue
except ImportError:
    import Queue as queue


def load_records(stringio_fp):
    records = list(fasta_iter(StringIO(stringio_fp.getvalue())))
    return records


def load_records_fastq(stringio_fp):
    records = list(fastq_iter(StringIO(stringio_fp.getvalue())))
    return records


def load_records_d(stringio_fp):
    return dict([(r.name, r.sequence)
                 for r in load_records(stringio_fp)])

# simple processing function: keep all sequences.


def idem(record):
    return record['name'], record['sequence']


# keep every *other* sequence
odd_counter = 0


def every_other(record):
    global odd_counter  # pylint: disable=global-statement
    odd_counter += 1
    if odd_counter % 2 == 1:
        return None, None

    return record['name'], record['sequence']

#


def test_basic():
    tsp = ThreadedSequenceProcessor(idem, 1, 1, verbose=False)

    inseqs = [screed.Record(name='a', sequence='AAA'),
              screed.Record(name='b', sequence='TTT'), ]
    outfp = StringIO()

    tsp.start(inseqs, outfp)

    x = load_records_d(outfp)
    assert len(x) == 2, x
    assert x['a'] == 'AAA'
    assert x['b'] == 'TTT'


def test_basic_fastq_like():
    tsp = ThreadedSequenceProcessor(idem, 1, 1, verbose=False)

    inseqs = [screed.Record(name='a', sequence='AAA', quality='###'),
              screed.Record(name='b', sequence='TTT', quality='###'), ]
    outfp = StringIO()

    tsp.start(inseqs, outfp)

    x = load_records_fastq(outfp)
    for i in x:
        assert i['quality'] == '###'


def test_odd():
    tsp = ThreadedSequenceProcessor(every_other, 1, 1, verbose=False)

    inseqs = [screed.Record(name='a', sequence='AAA'),
              screed.Record(name='b', sequence='TTT'), ]
    outfp = StringIO()

    tsp.start(inseqs, outfp)

    x = load_records_d(outfp)
    assert len(x) == 1, x
    assert x['b'] == 'TTT'


def test_basic_2thread():
    tsp = ThreadedSequenceProcessor(idem, 2, 1, verbose=False)

    inseqs = [screed.Record(name='a', sequence='AAA'),
              screed.Record(name='b', sequence='TTT'), ]
    outfp = StringIO()

    tsp.start(inseqs, outfp)

    x = load_records_d(outfp)
    assert len(x) == 2, x
    assert x['a'] == 'AAA'
    assert x['b'] == 'TTT'


def test_paired_2thread():
    class TSP_TestPairedProcess(ThreadedSequenceProcessor):
        # write a new do_process function that ensures paired ends are kept.

        def do_process(self):
            inq = self.inqueue

            while not self.done or not inq.empty():
                try:
                    g = inq.get(True, 1)
                except queue.Empty:
                    continue

                assert len(g.seqlist) == 2
                first_rec = g.seqlist[0]
                second_rec = g.seqlist[1]

                assert first_rec['name'][:-1] == second_rec['name'][:-1]
                assert first_rec['name'][-1] == '1'
                assert second_rec['name'][-1] == '2'

                keep = []
                for record in g.seqlist:
                    name, sequence = self.process_fn(record)
                    if name:
                        keep.append((name, sequence, None))

                self.outqueue.put(SequenceGroup(0, keep))

            # end of thread; exit, decrement worker count.
            self.worker_count -= 1

    #

    tsp = TSP_TestPairedProcess(idem, 1, 1, verbose=False)

    inseqs = [screed.Record(name='a/1', sequence='AAA'),
              screed.Record(name='a/2', sequence='TTT'), ]
    outfp = StringIO()

    tsp.start(inseqs, outfp)

    x = load_records_d(outfp)
    assert len(x) == 2, x
    assert x['a/1'] == 'AAA'
    assert x['a/2'] == 'TTT'


def test_paired_2thread_more_seq():
    class TSP_TestPairedProcess(ThreadedSequenceProcessor):
        # write a new do_process function that ensures paired ends are kept.

        def do_process(self):
            inq = self.inqueue

            while not self.done or not inq.empty():
                try:
                    g = inq.get(True, 1)
                except queue.Empty:
                    continue

                if len(g.seqlist) == 2:
                    first_rec = g.seqlist[0]
                    second_rec = g.seqlist[1]

                    assert first_rec['name'][:-1] == second_rec['name'][:-1]
                    assert first_rec['name'][-1] == '1'
                    assert second_rec['name'][-1] == '2'

                keep = []
                for record in g.seqlist:
                    name, sequence = self.process_fn(record)
                    if name:
                        keep.append((name, sequence, None))

                self.outqueue.put(SequenceGroup(0, keep))

            # end of thread; exit, decrement worker count.
            self.worker_count -= 1

    #

    tsp = TSP_TestPairedProcess(idem, 1, 1, verbose=False)

    inseqs = [screed.Record(name='b/1', sequence='AAA'),
              screed.Record(name='a/1', sequence='AAA'),
              screed.Record(name='a/2', sequence='TTT'),
              screed.Record(name='c/2', sequence='AAA'), ]
    outfp = StringIO()

    tsp.start(inseqs, outfp)

    x = load_records_d(outfp)
    assert len(x) == 4, x
    assert x['a/1'] == 'AAA'
    assert x['a/2'] == 'TTT'
    assert x['b/1'] == 'AAA'
    assert x['c/2'] == 'AAA'
