# Copyright (c) 2008-2010, Michigan State University

"""
screed is a database tool useful for retrieving arbitrary kinds of sequence
data through a on-disk database that emulates a read-only Python dictionary.

For opening a screed database, the 'ScreedDB' class is used. This class
accepts a string file path to a pre-created screed database. Read-only
dictionary methods are implemented here.

For creating a screed database, the 'create_db' function is used. This
function accepts an iterator as an argument which will yield records
from its respective sequence file. create_db will sequentially pull
records from the iterator, writing them to disk in a screed database
until the iterator is done.

Automatic ways for parsing FASTA and FASTQ files are accessed through
the read_fast*_sequences functions. These parse the given sequence
file into a screed database.

Conversion between sequence file types is provided in the ToFastq and
ToFasta functions
"""
from openscreed import ScreedDB, open_writer
from openscreed import open_reader as open
from conversion import ToFastq
from conversion import ToFasta
from createscreed import create_db
from seqparse import read_fastq_sequences
from seqparse import read_fasta_sequences
from dna import rc

__version__ = '0.7'
