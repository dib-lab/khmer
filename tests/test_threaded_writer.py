from khmer.thread_utils import ThreadedProcessor, PairThreadedProcessor
from cStringIO import StringIO

def fasta_test():
    fp = StringIO()
    tw = ThreadedProcessor(fp, fastq=False)

    tw.writer.save('a', 'ATCG')
    tw.exit()
    tw.run()

    x = fp.getvalue()
    assert x == '>a\nATCG\n'

def fastq_test():
    fp = StringIO()
    tw = ThreadedProcessor(fp, fastq=True)

    tw.writer.save('a', 'ATCG', '####')
    tw.exit()
    tw.run()

    x = fp.getvalue()
    assert x == '@a\nATCG\n+\n####\n', x

def fastq_sniff_test():
    fp = StringIO()
    tw = ThreadedProcessor(fp)

    tw.writer.save('a', 'ATCG', '####')
    tw.exit()
    tw.run()

    x = fp.getvalue()
    assert x == '@a\nATCG\n+\n####\n', x

def fasta_sniff_test():
    fp = StringIO()
    tw = ThreadedProcessor(fp)

    tw.writer.save('a', 'ATCG')
    tw.exit()
    tw.run()

    x = fp.getvalue()
    assert x == '>a\nATCG\n'

def test_several():
    fp = StringIO()
    tw = ThreadedProcessor(fp)

    tw.writer.save('a', 'ATCG')
    tw.writer.save('a', 'ATCG')
    tw.writer.save('a', 'ATCG')

    tw.exit()
    tw.run()

    assert len(fp.getvalue().splitlines()) == 6, fp.getvalue()

def pair_fasta_test():
    fp = StringIO()
    tw = PairThreadedProcessor(fp, fastq=False)

    tw.writer.save(('a', 'ATCG'), ('b', 'TAGC'))
    tw.exit()
    tw.run()
        
    x = fp.getvalue()
    assert x == '>a\nATCG\n>b\nTAGC\n'

def pair_fastq_test():
    fp = StringIO()
    tw = PairThreadedProcessor(fp, fastq=True)

    tw.writer.save(('a', 'ATCG', '####'), ('b', 'TAGC', 'BBBB'))
    tw.exit()
    tw.run()

    x = fp.getvalue()
    assert x == '@a\nATCG\n+\n####\n@b\nTAGC\n+\nBBBB\n', x

def pair_fastq_sniff_test():
    fp = StringIO()
    tw = PairThreadedProcessor(fp)

    tw.writer.save(('a', 'ATCG', '####'), ('b', 'TAGC', 'BBBB'))
    tw.exit()
    tw.run()
    
    x = fp.getvalue()
    assert x == '@a\nATCG\n+\n####\n@b\nTAGC\n+\nBBBB\n', x

def pair_fasta_sniff_test():
    fp = StringIO()
    tw = PairThreadedProcessor(fp)

    tw.writer.save(('a', 'ATCG'), ('b', 'TAGC'))
    tw.exit()
    tw.run()

    x = fp.getvalue()
    assert x == '>a\nATCG\n>b\nTAGC\n', x
    

def pair_test_several():
    fp = StringIO()
    tw = PairThreadedProcessor(fp)

    tw.writer.save(('a', 'ATCG'), ('b', 'TAGC'))
    tw.writer.save(('a', 'ATCG'), ('b', 'TAGC'))
    tw.writer.save(('a', 'ATCG'), ('b', 'TAGC'))

    tw.exit()
    tw.run()

    assert len(fp.getvalue().splitlines()) == 12, fp.getvalue()
