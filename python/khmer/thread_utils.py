"""
Utilities for dealing with multithreaded processing of short reads.
"""
import threading
import Queue
import sys
import screed
import time

DEFAULT_WORKER_THREADS = 8
DEFAULT_GROUPSIZE = 100


def verbose_loader(filename):
    it = screed.open(filename)
    for n, record in enumerate(it):
        if n % 100000 == 0:
            print >>sys.stderr, '... filtering', n
        yield record

verbose_fasta_iter = verbose_loader


class SequenceGroup(object):
    def __init__(self, order, seqlist):
        self.order = order
        self.seqlist = seqlist


def is_pair(r1, r2):
    a = r1['name'].split('/')[0]
    b = r2['name'].split('/')[0]

    return (a == b)

class FakeQueue(object):
    def __init__(self, size):
        self.q = []

    def put(self, item):
        self.q.append(item)

    def empty(self):
        return not len(self.q)

    def get(self, raise_on_empty, wait):
        if self.q:
            return self.q.pop()
        time.sleep(wait)
        if self.q:
            return self.q.pop()
            
        raise Queue.Empty

def start_threads(n_threads, target, args):
    "Start up n_threads, running the given target with the given args."
    threads = []
    for tnum in xrange(n_threads):
        t = threading.Thread(target=target, args=args)
        threads.append(t)
        t.start()

    return threads

class ThreadedWriter(object):
    """A queue-based threadsafe FASTA/FASTQ output writer.

    Public API:
       x = ThreadedWriter(outfp) - write records to outfp, sniffing format

       x.start()             - start self up in thread.
       x.join(other_threads) - wait on other threads, then flush & exit

    Typical usage:
       x = ThreadedWriter(outfp).start()
       ...
       x.join(other_threads)
    """
    QUEUESIZE = 1000                    # what should this be? @CTB

    def __init__(self, outfp, fastq=None):
        self.fp = outfp
        #self.outq = Queue.Queue(self.QUEUESIZE)
        self.outq = FakeQueue(self.QUEUESIZE)
        self.fastq = fastq              # None => sniff; else T(fastq) / F(fa)
        self._exit = False
        self._thread = None
        self.report_every = 10000
        self.report_lock = threading.Lock()
        self.n_read = self.n_saved = self.bp_read = self.bp_saved = 0
        self.last_report_n_read = 0

    def save_record(self, record):
        """Queue up the given record for output.

        Note: run in non-writer threads.
        """
        if self.fastq is None:
            if record.accuracy:
                self.fastq = True
            else:
                self.fastq = False
                
        if self.fastq:
            self.save(record.name, record.sequence, record.accuracy)
        else:
            self.save(record.name, record.sequence)

    def save(self, name, sequence, quality=None):
        """Queue up the given sequence for output.

        Note: run in non-writer threads.
        """
        if self.fastq is None:
            if quality:
                self.fastq = True
            else:
                self.fastq = False

        if self.fastq and quality is None:
            raise Exception("Error: empty quality given for read %s" % name)

        self.outq.put((name, sequence, quality))

    ###

    def start(self):
        """Start the writer thread, run self.run()."""
        if self._thread:
            raise Exception("writer thread already running!")
        self._thread = threading.Thread(target=self.run)
        self._thread.start()
        return self

    def join(self, other_threads):
        """Wait until reader threads are finished, then flush & exit.

        Note: run in writer thread.
        """
        for t in other_threads:
            t.join()
            
        self.exit()
        self._thread.join()

    def run(self):
        """Run until _exit flag is set & queue is empty.

        Note: run in writer thread.
        """
        while (not self._exit) or (not self.outq.empty()):
            self._looper()

    def exit(self):
        """Set the exit flag - no more 'save's will be called.

        Can be called from any thread.
        """
        self._exit = True

    ## internal functions

    def _looper(self):
        """Loop until the queue is empty and a timeout is received.

        Note: run in writer thread.
        """
        try:
            while 1:
                # wait for the queue to stay empty and then raise Empty
                item = self.outq.get(True, 0.0001)
                self._write(item)
        except Queue.Empty:
            return

    def _write(self, item):
        """Write out a single record in appropriate FASTA/FASTQ format.

        Note: run in writer thread.
        """
        
        if self.fastq:
            self.fp.write('@%s\n%s\n+\n%s\n' % (item[0], item[1], item[2]))
        else:
            self.fp.write('>%s\n%s\n' % (item[0], item[1]))

    def process_fn(self, rparser, filter_fn):
        """\
        Handle single reads one at a time.

        If 'filter_fn' returns anything, it's assumed to be sequence to be
        saved; if quality score is present, it's trimmed appropriately and
        then sent to the queue.

        Note: run in non-writer threads.
        """
        # note, all local variables => threadsafe
        n_read = 0
        n_saved = 0
        bp_read = 0
        bp_saved = 0
        
        for read in rparser:
            n_read += 1
            bp_read += len(read.sequence)
            
            seq = filter_fn(read.name, read.sequence)
            if seq:                     # save?
                n_saved += 1
                bp_saved += len(seq)
                
                if read.accuracy:
                    accuracy = read.accuracy[:len(seq)]
                    self.save(read.name, seq, accuracy)
                else:
                    self.save(read.name, seq)

            if n_read % self.report_every == 0:
                self.report(n_read, bp_read, n_saved, bp_saved)
                n_read = n_saved = bp_read = bp_saved = 0
                
        self.report(n_read, bp_read, n_saved, bp_saved, force=True)

    def report(self, n_read, bp_read, n_saved, bp_saved, force=False):
        with self.report_lock:
            since = self.n_read - self.last_report_n_read
            
            self.n_read += n_read
            self.bp_read += bp_read
            self.n_saved += n_saved
            self.bp_saved += bp_saved

            if since + n_read >= self.report_every or force: # report!
                n = self.n_read
                self.last_report_n_read = n

                print \
                 "... read %d sequences (%d kb); wrote %d sequences (%d kb)" %\
                     (self.n_read, round(self.bp_read/1000.),
                      self.n_saved, round(self.bp_saved/1000.))
            

class PairThreadedWriter(ThreadedWriter):
    """A queue-based threadsafe FASTA/FASTQ output writer.

    Public API:
       x = PairThreadedWriter(outfp) - write records to outfp, sniffing format

       x.start()             - start self up in thread.
       x.join(other_threads) - wait on other threads, then flush & exit

    Typical usage:
       x = PairThreadedWriter(outfp).start()
       ...
       x.join(other_threads)
    """
    QUEUESIZE = 1000                    # what should this be? @CTB

    def save(self, a, b):
        """Add a pair of sequences to be output.

        Note: run in non-writer threads.
        """
        if self.fastq is None:
            if len(a) == 3:
                self.fastq = True
            else:
                self.fastq = False

        if self.fastq and (not a[2] or not b[2]):
            raise Exception("Error: empty quality given for pair %s" % a[0])

        self.outq.put((a, b))

    ###


    def _write(self, item):
        """Write out a pair of records in appropriate FASTA/FASTQ format.

        Note: run in non-writer threads.
        """

        a, b = item
        if self.fastq:
            self.fp.write('@%s\n%s\n+\n%s\n' % (a[0], a[1], a[2]))
            self.fp.write('@%s\n%s\n+\n%s\n' % (b[0], b[1], b[2]))
        else:
            self.fp.write('>%s\n%s\n' % (a[0], a[1]))
            self.fp.write('>%s\n%s\n' % (b[0], b[1]))

    def process_fn(self, rparser, filter_fn):
        """Apply filter function; write out sequences appropriately.

        Note: run in non-writer threads.
        """
        # note, all local variables => threadsafe
        n_read = 0
        n_saved = 0
        bp_read = 0
        bp_saved = 0
        
        for r1, r2 in rparser.iter_read_pairs():
            n_read += 2
            bp_read += len(r1.sequence) + len(r2.sequence)
            
            seq1, seq2 = filter_fn(r1.name, r1.sequence,
                                   r2.name, r2.sequence)
            if seq1 and seq2:
                n_saved += 2
                bp_saved += len(seq1) + len(seq2)
                
                if r1.accuracy:
                    accuracy1 = r1.accuracy[:len(seq1)]
                    accuracy2 = r2.accuracy[:len(seq2)]
                    self.save((r1.name, seq1, accuracy1),
                              (r2.name, seq2, accuracy2))
                else:
                    self.save((r1.name, seq1),
                              (r2.name, seq2))

            if n_read % self.report_every == 0:
                self.report(n_read, bp_read, n_saved, bp_saved)
                n_read = n_saved = bp_read = bp_saved = 0
                
        self.report(n_read, bp_read, n_saved, bp_saved, force=True)

class ThreadedSequenceProcessor(object):
    """
    DEPRECATED
    """
    
    QUEUESIZE = 50

    def __init__(self, process_fn, n_workers=DEFAULT_WORKER_THREADS,
                 group_size=DEFAULT_GROUPSIZE, verbose=True):
        self.process_fn = process_fn
        self.n_workers = n_workers
        self.group_size = group_size

        self.inqueue = Queue.Queue(self.QUEUESIZE)
        self.outqueue = Queue.Queue(self.QUEUESIZE)

        self.worker_count = 0
        self.worker_count_lock = threading.Lock()
        self.done = False
        self.verbose = verbose

        self.n_processed = 0
        self.n_written = 0
        self.bp_processed = 0
        self.bp_written = 0
        self.tallies_lock = threading.Lock()

    def start(self, inputiter, outfp):
        if self.verbose:
            print >>sys.stderr, 'starting threads'

        try:
            for i in range(self.n_workers):
                t = threading.Thread(target=self.do_process)
                self.worker_count += 1
                t.start()

            if self.verbose:
                print >>sys.stderr, 'starting writer'

            w = threading.Thread(target=self.do_write, args=(outfp,))
            w.start()

            if self.verbose:
                print >>sys.stderr, 'loading...'

            self.push_sequences(inputiter)

            if self.verbose:
                print >>sys.stderr, 'done loading in sequences'
            self.done = True

            w.join()
        except:
            self.done = True
            raise

    def push_sequences(self, inputiter):
        batch = []
        last_record = None
        i = 0
        for record in inputiter:
            if i >= self.group_size:
                # keep pairs together in batches, to retain the interleaving.
                if is_pair(record, last_record):
                    batch.append(record)
                    g = SequenceGroup(0, batch)
                    self.inqueue.put(g)

                    batch = []
                else:
                    g = SequenceGroup(0, batch)
                    self.inqueue.put(g)
                    batch = [record]

                i = 0
            else:
                batch.append(record)

            last_record = record
            i += 1

        # submit last set of sequences
        if batch:
            g = SequenceGroup(0, batch)
            self.inqueue.put(g)

    def do_process(self):
        inq = self.inqueue
        outq = self.outqueue

        while not self.done or not inq.empty():
            try:
                g = inq.get(True, 1)
            except Queue.Empty:
                continue

            bp_processed = 0
            bp_written = 0

            keep = []
            for record in g.seqlist:
                name, sequence = self.process_fn(record)
                bp_processed += len(record['sequence'])
                if name:
                    accuracy = record.get('accuracy')
                    if accuracy:
                        accuracy = accuracy[:len(sequence)]
                    bp_written += len(sequence)
                    keep.append((name, sequence, accuracy))

            self.outqueue.put(SequenceGroup(0, keep))

            # the tallies are shared among workers, hence we lock
            with self.tallies_lock:

                self.n_processed += len(g.seqlist)
                self.n_written += len(keep)
                self.bp_processed += bp_processed
                self.bp_written += bp_written

                if self.verbose and self.n_processed % 500000 == 0:
                    print >>sys.stderr, \
                        "processed %d / wrote %d / removed %d" % \
                        (self.n_processed, self.n_written,
                         self.n_processed - self.n_written)
                    print >>sys.stderr, \
                        "processed %d bp / wrote %d bp / removed %d bp" % \
                        (self.bp_processed, self.bp_written,
                         self.bp_processed - self.bp_written)
                    discarded = self.bp_processed - self.bp_written
                    f = float(discarded) / float(self.bp_processed) * 100
                    print >>sys.stderr, "discarded %.1f%%" % f

        # end of thread; exit, decrement worker count.
        with self.worker_count_lock:
            self.worker_count -= 1

    def do_write(self, outfp):
        outq = self.outqueue
        while self.worker_count > 0 or not outq.empty():
            try:
                g = outq.get(True, 1)
            except Queue.Empty:
                continue

            for name, seq, accuracy in g.seqlist:
                if accuracy: # write FASTQ; CTB hack.
                    outfp.write('@%s\n%s\n+\n%s\n' % (name, seq, accuracy))
                else:
                    outfp.write('>%s\n%s\n' % (name, seq,))

        if self.verbose:
            print >>sys.stderr, \
                "DONE writing.\nprocessed %d / wrote %d / removed %d" % \
                (self.n_processed, self.n_written,
                 self.n_processed - self.n_written)
            print >>sys.stderr, \
                "processed %d bp / wrote %d bp / removed %d bp" % \
                (self.bp_processed, self.bp_written,
                 self.bp_processed - self.bp_written)
            discarded = self.bp_processed - self.bp_written
            f = float(discarded) / float(self.bp_processed) * 100
            print >>sys.stderr, "discarded %.1f%%" % f

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
