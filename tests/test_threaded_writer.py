from khmer.thread_utils import ThreadedWriter, PairThreadedWriter
from cStringIO import StringIO

def fasta_test():
    fp = StringIO()
    tw = ThreadedWriter(fp, fastq=False)

    tw.save('a', 'ATCG')
    tw._looper()

    x = fp.getvalue()
    assert x == '>a\nATCG\n'

def fastq_test():
    fp = StringIO()
    tw = ThreadedWriter(fp, fastq=True)

    tw.save('a', 'ATCG', '####')
    tw._looper()

    x = fp.getvalue()
    assert x == '@a\nATCG\n+\n####\n', x

def fastq_sniff_test():
    fp = StringIO()
    tw = ThreadedWriter(fp)

    tw.save('a', 'ATCG', '####')
    tw._looper()

    x = fp.getvalue()
    assert x == '@a\nATCG\n+\n####\n', x

def fasta_sniff_test():
    fp = StringIO()
    tw = ThreadedWriter(fp)

    tw.save('a', 'ATCG')
    tw._looper()

    x = fp.getvalue()
    assert x == '>a\nATCG\n'

def test_several():
    fp = StringIO()
    tw = ThreadedWriter(fp)

    tw.save('a', 'ATCG')
    tw.save('a', 'ATCG')
    tw.save('a', 'ATCG')

    tw.exit()
    tw.run()

    assert len(fp.getvalue().splitlines()) == 6, fp.getvalue()

def pair_fasta_test():
    fp = StringIO()
    tw = PairThreadedWriter(fp, fastq=False)

    tw.save(('a', 'ATCG'), ('b', 'TAGC'))
    tw._looper()

    x = fp.getvalue()
    assert x == '>a\nATCG\n>b\nTAGC\n'

def pair_fastq_test():
    fp = StringIO()
    tw = PairThreadedWriter(fp, fastq=True)

    tw.save(('a', 'ATCG', '####'), ('b', 'TAGC', 'BBBB'))
    tw._looper()

    x = fp.getvalue()
    assert x == '@a\nATCG\n+\n####\n@b\nTAGC\n+\nBBBB\n', x

def pair_fastq_sniff_test():
    fp = StringIO()
    tw = PairThreadedWriter(fp)

    tw.save(('a', 'ATCG', '####'), ('b', 'TAGC', 'BBBB'))
    tw._looper()

    x = fp.getvalue()
    assert x == '@a\nATCG\n+\n####\n@b\nTAGC\n+\nBBBB\n', x

def pair_fasta_sniff_test():
    fp = StringIO()
    tw = PairThreadedWriter(fp)

    tw.save(('a', 'ATCG'), ('b', 'TAGC'))
    tw._looper()

    x = fp.getvalue()
    assert x == '>a\nATCG\n>b\nTAGC\n', x
    

def pair_test_several():
    fp = StringIO()
    tw = PairThreadedWriter(fp)

    tw.save(('a', 'ATCG'), ('b', 'TAGC'))
    tw.save(('a', 'ATCG'), ('b', 'TAGC'))
    tw.save(('a', 'ATCG'), ('b', 'TAGC'))

    tw.exit()
    tw.run()

    assert len(fp.getvalue().splitlines()) == 12, fp.getvalue()
