import sys
from khmer.thread_utils import ThreadedSequenceProcessor, SequenceGroup
from io import StringIO
from screed.fasta import fasta_iter
from screed.fastq import fastq_iter
import queue
from nose.plugins.attrib import attr


def load_records(stringio_fp):
    records = list(fasta_iter(StringIO(stringio_fp.getvalue())))
    return records


def load_records_fastq(stringio_fp):
    records = list(fastq_iter(StringIO(stringio_fp.getvalue())))
    return records


def load_records_d(stringio_fp):
    return dict([(r['name'], r['sequence'])
                 for r in load_records(stringio_fp)])

# simple processing function: keep all sequences.


def idem(record):
    return record['name'], record['sequence']

# keep every *other* sequence
odd_counter = 0


def every_other(record):
    global odd_counter
    odd_counter += 1
    if odd_counter % 2 == 1:
        return None, None

    return record['name'], record['sequence']

#


def test_basic():
    tsp = ThreadedSequenceProcessor(idem, 1, 1, verbose=False)

    input = [dict(name='a', sequence='AAA'),
             dict(name='b', sequence='TTT'), ]
    outfp = StringIO()

    tsp.start(input, outfp)

    x = load_records_d(outfp)
    assert len(x) == 2, x
    assert x['a'] == 'AAA'
    assert x['b'] == 'TTT'


def test_basic_fastq_like():
    tsp = ThreadedSequenceProcessor(idem, 1, 1, verbose=False)

    input = [dict(name='a', sequence='AAA', quality='###'),
             dict(name='b', sequence='TTT', quality='###'), ]
    outfp = StringIO()

    tsp.start(input, outfp)

    x = load_records_fastq(outfp)
    for i in x:
        assert i['quality'] == '###'


def test_odd():
    tsp = ThreadedSequenceProcessor(every_other, 1, 1, verbose=False)

    input = [dict(name='a', sequence='AAA'),
             dict(name='b', sequence='TTT'), ]
    outfp = StringIO()

    tsp.start(input, outfp)

    x = load_records_d(outfp)
    assert len(x) == 1, x
    assert x['b'] == 'TTT'


def test_basic_2thread():
    tsp = ThreadedSequenceProcessor(idem, 2, 1, verbose=False)

    input = [dict(name='a', sequence='AAA'),
             dict(name='b', sequence='TTT'), ]
    outfp = StringIO()

    tsp.start(input, outfp)

    x = load_records_d(outfp)
    assert len(x) == 2, x
    assert x['a'] == 'AAA'
    assert x['b'] == 'TTT'


def test_paired_2thread():
    class TSP_TestPairedProcess(ThreadedSequenceProcessor):
        # write a new do_process function that ensures paired ends are kept.

        def do_process(self):
            inq = self.inqueue
            outq = self.outqueue

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

    input = [dict(name='a/1', sequence='AAA'),
             dict(name='a/2', sequence='TTT'), ]
    outfp = StringIO()

    tsp.start(input, outfp)

    x = load_records_d(outfp)
    assert len(x) == 2, x
    assert x['a/1'] == 'AAA'
    assert x['a/2'] == 'TTT'


def test_paired_2thread_more_seq():
    class TSP_TestPairedProcess(ThreadedSequenceProcessor):
        # write a new do_process function that ensures paired ends are kept.

        def do_process(self):
            inq = self.inqueue
            outq = self.outqueue

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

    input = [dict(name='b/1', sequence='AAA'),
             dict(name='a/1', sequence='AAA'),
             dict(name='a/2', sequence='TTT'),
             dict(name='c/2', sequence='AAA'), ]
    outfp = StringIO()

    tsp.start(input, outfp)

    x = load_records_d(outfp)
    assert len(x) == 4, x
    assert x['a/1'] == 'AAA'
    assert x['a/2'] == 'TTT'
    assert x['b/1'] == 'AAA'
    assert x['c/2'] == 'AAA'
