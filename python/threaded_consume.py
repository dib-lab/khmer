from threading import Thread
import time
import khmer

class ConsumeThread(Thread):
    def __init__(self, wordsize, genome):
        self.wordsize = wordsize
        self.genome = genome
        Thread.__init__(self)

    def run(self):
        self.kt  = khmer.consume_genome(self.wordsize, self.genome)

genome = open('/tmp/all2.dna').read()

nthreads = 2
length = len(genome)
genome1 = genome[:length/2]
genome2 = genome[length/2:]
t1 = ConsumeThread(5, genome1)
t2 = ConsumeThread(5, genome2)

t1.start()
t2.start()
t1.join()
t2.join()

master_kt = khmer.new_ktable(5)
master_kt.update(t1.kt)
master_kt.update(t2.kt)
