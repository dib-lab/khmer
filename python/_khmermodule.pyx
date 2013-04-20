from libc.stdint cimport uint8_t, uint32_t, uint64_t
from libcpp.string cimport string
from libcpp.map cimport map
from libcpp.pair cimport pair
from libcpp.set cimport set
from libcpp.vector cimport vector
from libcpp cimport bool
from cpython cimport bool as pybool
from libc.string cimport strlen
from cython.operator cimport dereference as deref, preincrement as inc

ctypedef unsigned char WordLength
ctypedef unsigned long long int HashIntoType
ctypedef unsigned long long int ExactCounterType
ctypedef unsigned short int BoundedCounterType
ctypedef unsigned char Byte
ctypedef unsigned int IntersectionID
ctypedef unsigned int PartitionID

cdef extern from *:
   ctypedef char* const_char_ptr "const char*"

ctypedef void (*CallbackFn)(const_char_ptr info, 
                            void* callback_data, 
                            unsigned long long n_reads,
                            unsigned long long other) except *

DEF MAX_BIGCOUNT=65535

cdef extern from "stdlib.h":
   ctypedef unsigned long size_t
   void free(void *ptr)
   void *malloc(size_t size)
   void *realloc(void *ptr, size_t size)
   char *strcpy(char *dest, char *scr)

cdef extern from "../lib/storage.hh" namespace "khmer":
   ctypedef struct MinMaxValue:
      BoundedCounterType min_val
      BoundedCounterType max_val

cdef extern from "string" namespace "std":
   cdef cppclass stringclass "std::string":
      void stringclass()
      char *c_str

cdef extern from "string" namespace "std":
   cdef cppclass string:
      string()
      string(char*) 
      char* c_str()
      unsigned int length()


cdef extern from "../lib/khmer_config.hh"	namespace "khmer":
    
  cdef cppclass Config:
      Config( )

      bool has_extra_sanity_checks( )

      uint32_t get_number_of_threads( )
      void set_number_of_threads( uint32_t )

      uint64_t get_reads_input_buffer_size( )
      void set_reads_input_buffer_size( uint64_t )

      uint8_t get_input_buffer_trace_level( )
      void set_input_buffer_trace_level( uint8_t )
      uint8_t get_reads_parser_trace_level( )
      void set_reads_parser_trace_level( uint8_t )

  cdef Config get_active_config( )

cdef class get_config:

   cdef Config *     thisref

   def __cinit__( self ):
      # NOTE! The following code is wrong.
      #       Changes to the active config in the C++ API will not be 
      # reflected in the Python wrapper.
      #       This is due to Cython's pathetic ability to handle C++ 
      # references.
      self.thisref = new Config( )

   def has_extra_sanity_checks( self ):
      return self.thisref.has_extra_sanity_checks( )

   def get_number_of_threads( self ): 
      return self.thisref.get_number_of_threads( )
   def set_number_of_threads( self, uint32_t number_of_threads ):
      # TODO? Handle exceptions from C++ side.
      self.thisref.set_number_of_threads( number_of_threads )

   def get_reads_input_buffer_size( self ):
      return self.thisref.get_reads_input_buffer_size( )
   def set_reads_input_buffer_size( self, uint64_t reads_input_buffer_size ):
      # TODO? Handle exceptions from C++ code.
      self.thisref.set_reads_input_buffer_size( reads_input_buffer_size )

   def get_input_buffer_trace_level( self ):
      return self.thisref.get_input_buffer_trace_level( )
   def set_input_buffer_trace_level( self, uint8_t trace_level ):
      self.thisref.set_input_buffer_trace_level( trace_level )
   def get_reads_parser_trace_level( self ):
      return self.thisref.get_reads_parser_trace_level( )
   def set_reads_parser_trace_level( self, uint8_t trace_level ):
      self.thisref.set_reads_parser_trace_level( trace_level )

cdef extern from "../lib/read_parsers.hh" namespace "khmer:: read_parsers":
   
   cdef cppclass Read:
      
      Read( ) except +
      
      string name
      string annotations
      string sequence
      string accuracy
      # TODO? Add 'description'.
    
   cdef cppclass IParser:

      #IParser( ) except + # TODO: Fill out parameters.

      bool is_complete( )
      Read get_next_read( )


cdef extern from "../lib/read_parsers.hh" \
namespace "khmer:: read_parsers:: IParser":
   
   cdef IParser * get_parser(
      string ifile_name, uint32_t number_of_threads, 
      uint64_t cache_size, uint8_t trace_level
   )


cdef class _Read:
   
   name        = ""
   annotations = ""
   sequence    = ""
   accuracy    = ""

   def __cinit__(
      self, char * name, char * annotations, char * sequence, char * accuracy
   ):
      self.name         = name
      self.annotations  = annotations
      self.sequence     = sequence
      self.accuracy     = accuracy


cdef class _ReadParser:

   cdef IParser *     thisref

   def __cinit__(
      self,
      char * ifile_name, uint32_t number_of_threads, uint64_t cache_size, uint8_t trace_level
   ): 
      cdef string    ifile_name_STRING = string( ifile_name )

      self.thisref = \
      get_parser(
         ifile_name_STRING, number_of_threads, cache_size, trace_level
      )

   def is_complete( self ):
      return self.thisref.is_complete( )

   def get_next_read( self ):
      cdef Read the_read
      the_read = self.thisref.get_next_read( )
      return _Read(
         the_read.name.c_str( ), the_read.annotations.c_str( ),
         the_read.sequence.c_str( ), the_read.accuracy.c_str( )
      ) 
  

cdef extern from "../lib/ktable.hh" namespace "khmer":
   HashIntoType _hash(char*, WordLength)
   HashIntoType _hash(char*, WordLength, HashIntoType&, HashIntoType&)
   HashIntoType _hash_forward(char*, WordLength)
   string _revhash(HashIntoType, WordLength)

   cdef cppclass KTable:
      KTable(WordLength ksize)
      KTable()

      WordLength ksize()
      HashIntoType max_hash()
      HashIntoType n_entries()
      void count(char *)
      ExactCounterType get_count(char *)
      ExactCounterType get_count(HashIntoType)
      void set_count(char *, ExactCounterType)
      void set_count(HashIntoType, ExactCounterType)
      void consume_string(char *)
      void clear()
      void update(KTable)
      KTable * intersect(KTable)

cdef extern from "../lib/hashtable.hh" namespace "khmer":
   cdef cppclass Hashtable:
      Hashtable(WordLength)

      void _init_bitstuff()
      HashIntoType _next_hash(char, HashIntoType, HashIntoType)
      WordLength ksize()
      BoundedCounterType get_count(char*)
      BoundedCounterType get_count(HashIntoType)
      unsigned int consume_string(char*, HashIntoType, HashIntoType)
      unsigned int check_read(char*)


cdef extern from "../lib/storage.hh" namespace "khmer":

   cdef cppclass MinMaxTable:
      MinMaxTable(unsigned int)      
      MinMaxTable()

      unsigned int get_tablesize()
      BoundedCounterType get_min(unsigned int)
      BoundedCounterType get_max(unsigned int)
      BoundedCounterType add_max(unsigned int, unsigned int)
      BoundedCounterType add_min(unsigned int, unsigned int)
      void clear(unsigned int)
      void merge(MinMaxTable)
      void save(char*)
      void load(char*)

cdef extern from "../lib/counting.hh" namespace "khmer":
   cdef cppclass Hashbits

   cdef cppclass CountingHash:
      CountingHash(WordLength, HashIntoType)
      CountingHash(WordLength, vector[HashIntoType])
      CountingHash()

      vector[HashIntoType] get_tablesizes()
      WordLength ksize()
      void count(char*)
      void save(char*)
      void load(char*)
      HashIntoType n_entries()
      HashIntoType n_occupied(HashIntoType, HashIntoType)
      BoundedCounterType get_count(char*)
      BoundedCounterType get_count(HashIntoType)
      MinMaxTable * fasta_file_to_minmax(char*, unsigned int, 
                                         CallbackFn,
                                         void*) except *
      void output_fasta_kmer_pos_freq(char*, char*)
      BoundedCounterType get_min_count(char*, HashIntoType, HashIntoType)
      BoundedCounterType get_max_count(char*, HashIntoType, HashIntoType)
      HashIntoType * abundance_distribution(char*, Hashbits* tracking,
                                            CallbackFn,
                                            void*) except *
      HashIntoType * fasta_count_kmers_by_position(char*, unsigned int,
                                                   BoundedCounterType,
                                                   CallbackFn,
                                                   void*) except *
      void fasta_dump_kmers_by_abundance(char*, BoundedCounterType,
                                         CallbackFn, 
                                         void*) except *
      unsigned int consume_string(char*, HashIntoType, HashIntoType)
      void consume_fasta(char*, unsigned int&, unsigned long long&,
                         HashIntoType, HashIntoType, 
                         CallbackFn, void*) except *
      void set_use_bigcount(bool)
      bool get_use_bigcount()
      void get_median_count(char*, BoundedCounterType&, 
                            float&, float&)
      void get_kadian_count(char*, BoundedCounterType&, unsigned int)

      void get_kmer_abund_mean(char*, unsigned long long &total, 
                               unsigned long long &count, float &mean)
      void get_kmer_abund_abs_deviation(char*, float mean, float &abs_deviation)
      unsigned int max_hamming1_count(char*)
      unsigned int trim_on_abundance(char*, BoundedCounterType)
      # TODO: Take care of corresponding test when this is modified.
      # void collect_high_abundance_kmers(char*, unsigned int, unsigned int, ) @CTB

cdef extern from "../lib/subset.hh" namespace "khmer":
   cdef cppclass Hashbits

   cdef cppclass pre_partition_info:
      pre_partition_info(HashIntoType)
      HashIntoType kmer
      set[HashIntoType] tagged_kmers

   cdef cppclass SubsetPartition:
      SubsetPartition(Hashbits* ht)
      PartitionID assign_partition_id(HashIntoType, set[HashIntoType])
      void set_partition_id(HashIntoType, PartitionID)
      void set_partition_id(char *, PartitionID)
      PartitionID join_partitions(PartitionID, PartitionID)
      PartitionID get_partition_id(char *)
      PartitionID get_partition_id(HashIntoType)
      PartitionID* get_new_partition()
      void merge(SubsetPartition *)
      void merge_from_disk(char *)
      void save_partitionmap(char *)
      void load_partitionmap(char *)
      void _validate_pmap()
      void find_all_tags(HashIntoType, HashIntoType, set[HashIntoType]&,
                         set[HashIntoType]&)
      void do_partition(HashIntoType, HashIntoType, bool, bool,
                        CallbackFn, void*) except *
      void count_partitions(unsigned int&, unsigned int&)
      unsigned int output_partitioned_file(char*, char*, bool,
                                           CallbackFn, void*) except *
      unsigned int find_unpart(char*, bool, bool, CallbackFn, void*)
      bool is_single_partition(char*)
      void join_partitions_by_path(char*)
      void partition_size_distribution(map[unsigned long long, unsigned long long], unsigned int&)
      unsigned int repartition_largest_partition(unsigned int, unsigned int,
                                                 unsigned int, CountingHash&)
      void repartition_a_partition(set[HashIntoType]&)

cdef extern from "../lib/hashbits.hh" namespace "khmer":
   cdef cppclass Hashbits:
      SubsetPartition * partition

      set[HashIntoType] all_tags

      Hashbits(WordLength, vector[HashIntoType])
      void save(char*)
      void load(char*)
      void save_tagset(char*)
      void load_tagset(char*, bool clear_tags)
      WordLength ksize()
      vector[HashIntoType] get_tablesizes()

      void calc_connected_graph_size(char *, unsigned long long,
                                     set[HashIntoType], unsigned long long,
                                     bool)
      unsigned int n_tags()
      void divide_tags_into_subsets(unsigned int, set[HashIntoType])
      void add_kmer_to_tags(HashIntoType)
      void clear_tags()
      void consume_fasta_and_tag(const_char_ptr, unsigned int&, unsigned long long&,
                                 CallbackFn,
                                 void*) except *
      void consume_partitioned_fasta(char*, unsigned int&, unsigned long long&,
                                     CallbackFn, 
                                     void*) except *
      void consume_fasta_and_traverse(char*, unsigned int, unsigned int, unsigned int,
                                      CountingHash&)
      void traverse_from_reads(char*, unsigned int, unsigned int, unsigned int,
                               CountingHash&)
      void consume_fasta_and_tag_with_stoptags(char* filename, unsigned int &total_reads,
                                               unsigned long long &n_consumed,
                                               CallbackFn, void*) except *
      void tags_to_map(map[HashIntoType, unsigned int])
      void discard_tags(map[HashIntoType, unsigned int], unsigned int)
      HashIntoType n_occupied(HashIntoType, HashIntoType)
      HashIntoType n_kmers(HashIntoType, HashIntoType)
      void count(char*)
      void count(HashIntoType)
      BoundedCounterType get_count(char*)
      BoundedCounterType get_count(HashIntoType)
      void filter_if_present(char*, char*, CallbackFn, void*) except *
      unsigned int count_kmers_within_radius(HashIntoType, HashIntoType,
                                            unsigned int, unsigned int)
      unsigned int consume_string(char*, HashIntoType, HashIntoType)
      void consume_fasta(char*, unsigned int&, unsigned long long&,
                         HashIntoType lower_bound, 
                         HashIntoType upper_bound, 
                         CallbackFn, void*) except *
      void add_tag(HashIntoType)
      void add_stop_tag(HashIntoType)
      void load_stop_tags(char* filename)
      void save_stop_tags(char*)
      void print_stop_tags(char* filename)
      unsigned int kmer_degree(HashIntoType, HashIntoType)
      unsigned int kmer_degree(char*)
      void _set_tag_density(unsigned int)
      unsigned int _get_tag_density()
      unsigned int trim_on_degree(char*, unsigned int)
      unsigned int trim_on_sodd(char*, unsigned int)
      unsigned int trim_on_density_explosion(char*, unsigned int, unsigned int)
      unsigned int trim_on_stoptags(char*)
      void identify_stop_tags_by_position(char*, vector[unsigned int]&)
      void hitraverse_to_stoptags(char*, CountingHash&, unsigned int)
      unsigned int find_radius_for_volume(HashIntoType, HashIntoType, unsigned int, unsigned int)
      unsigned int count_kmers_on_radius(HashIntoType, HashIntoType, unsigned int, unsigned int)
      void extract_unique_paths(char*, unsigned int, float, vector[string]&)

_callback_obj = None

cdef void _report_fn(const_char_ptr info, void* data, 
                     unsigned long long n_reads, unsigned long long other) except *:
   global _callback_obj
   _callback_obj(info, n_reads, other)

def set_reporting_callback(o):
   global _callback_obj
   _callback_obj = o

cdef class _pre_partition_info:
   cdef pre_partition_info *thisptr
   def __cinit__(self, HashIntoType h):
      self.thisptr = new pre_partition_info(h)

cdef class new_ktable:
   cdef KTable *thisptr
   def __cinit__(self, arg1, arg2=0):
      if arg2 == 0:
         self.thisptr = new KTable(<WordLength>arg1)
      else:
         k1 = <new_ktable>arg1
         k2 = <new_ktable>arg2
         self.thisptr = k1.thisptr.intersect(k2.thisptr[0])
   def count(self, kmer):
      self.thisptr.count(kmer)
   def consume(self, char* s):
      self.thisptr.consume_string(s)
      return strlen(s) - self.thisptr.ksize() + 1
   def get(self, arg):
      if isinstance(arg, str):
         return self.thisptr.get_count(<char*>arg)
      else:
         return self.thisptr.get_count(<HashIntoType>arg)
   def __getitem__(self, arg):
      if isinstance(arg, str):
         return self.thisptr.get_count(<char*>arg)
      else:
         return self.thisptr.get_count(<HashIntoType>arg)
   def set(self, HashIntoType i, ExactCounterType c):
      self.thisptr.set_count(i, c)
   def set(self, char* kmer, ExactCounterType c):
      self.thisptr.set_count(kmer, c)
   def forward_hash(self, char* s):
      return _hash(s, self.thisptr.ksize())
   def forward_hash_no_rc(self, char* s):
      return _hash_forward(s, self.thisptr.ksize())
   def reverse_hash(self, HashIntoType h):
      cdef string s = _revhash(h, self.thisptr.ksize())
      cdef char *cstr = <char *>malloc(s.length() + 1)
      strcpy(cstr, s.c_str())
      return cstr
   def max_hash(self):
      return self.thisptr.max_hash()
   def n_entries(self):
      return self.thisptr.n_entries()
   def __len__(self):
      return self.thisptr.n_entries()
   def ksize(self):
      return self.thisptr.ksize()
   def clear(self):
      self.thisptr.clear()
   def update(self, new_ktable ktbl2):
      self.thisptr.update(ktbl2.thisptr[0])
   def intersect(self, new_ktable ktbl2):
      result = new_ktable(self, ktbl2)
      return result
   def __contains__( self, char *kmer ):
      return int( pybool( self.thisptr.get_count( kmer ) ) )

cdef class new_minmax:
   cdef MinMaxTable *thisptr
   def __cinit__(self, unsigned int n=0):
      self.thisptr = new MinMaxTable(n)
   def save(self, char* s):
      self.thisptr.save(s)
   def load(self, char* s):
      self.thisptr.load(s)
   def merge(self, new_minmax tbl):
      self.thisptr.merge(tbl.thisptr[0])
   def clear(self, unsigned int n):
      self.thisptr.clear(n)
   def get_min(self, unsigned int n):
      return self.thisptr.get_min(n)
   def get_max(self, unsigned int n):
      return self.thisptr.get_max(n)
   def add_max(self, unsigned int n1, n2):
      return self.thisptr.add_max(n1, n2)
   def add_min(self, unsigned int n1, n2):
      return self.thisptr.add_min(n1, n2)
   def tablesize(self):
      return self.thisptr.get_tablesize()

cdef class _new_hashbits

cdef class _new_counting_hash:
   cdef CountingHash *thisptr
   # TODO: Handle the 'n_threads' argument.
   def __cinit__(self, WordLength ksize, arg, n_threads = 1):
      cdef vector[HashIntoType] v = vector[HashIntoType]()
      if isinstance(arg, long) or isinstance(arg, int):
         self.thisptr = new CountingHash(ksize, <HashIntoType>arg)
      else:
         primes = <list>arg
         for prime in primes:
            v.push_back(prime)
         self.thisptr = new CountingHash(ksize, v)
   def save(self, char* s):
      self.thisptr.save(s)
   def load(self, char* s):
      self.thisptr.load(s)
   def n_occupied(self, start=0, stop=0):
      return self.thisptr.n_occupied(start, stop)
   def n_entries(self):
      return self.thisptr.n_entries()
   def ksize(self):
      return self.thisptr.ksize()
   def consume(self, char* s, HashIntoType lower = 0, HashIntoType upper = 0):
      return self.thisptr.consume_string(s, lower, upper)
   def count(self, char* s):
      self.thisptr.count(s)
   def get_count(self, char* s):
      return self.thisptr.get_count(s)
   def get(self, arg):
      if isinstance(arg, str):
         return self.thisptr.get_count(<char*>arg)
      else:
         return self.thisptr.get_count(<HashIntoType>arg)
   def __getitem__(self, arg):
      if isinstance(arg, str):
         return self.thisptr.get_count(<char*>arg)
      else:
         return self.thisptr.get_count(<HashIntoType>arg)
   def get_max_count(self, char* s, lower=0, upper=0):
      return self.thisptr.get_max_count(s, lower, upper)
   def get_min_count(self, char* s, lower=0, upper=0):
      return self.thisptr.get_min_count(s, lower, upper)
   def hashsizes(self):
      cdef vector[HashIntoType] res = self.thisptr.get_tablesizes()
      hashes = []

      cdef vector[HashIntoType].iterator it = res.begin()
      while it != res.end():
         hashes.append(deref(it))
         inc(it)
      return hashes
   def output_fasta_kmer_pos_freq(self, char* input, char* output):
      self.thisptr.output_fasta_kmer_pos_freq(input, output)
   def abundance_distribution(self, char* filename, _new_hashbits hb, callback_obj=None):
      cdef Hashbits* hbits = <Hashbits*>hb.thisptr

      global _callback_obj
      if callback_obj is not None:
         _callback_obj = callback_obj

      dist = []
      cdist = self.thisptr.abundance_distribution(filename, hbits, _report_fn, <void*>_callback_obj)
      for i in range(MAX_BIGCOUNT+1):
         dist.append(cdist[i])
      return dist
   def fasta_file_to_minmax(self, char* inputfile, unsigned int total_reads, 
                            callback_obj=None):
      global _callback_obj
      if callback_obj is not None:
         _callback_obj = callback_obj

      cdef MinMaxTable * ret_minmax = self.thisptr.fasta_file_to_minmax(inputfile,
                                                                        total_reads,
                                                                        _report_fn,
                                                                        <void*>_callback_obj)
      return_minmax = new_minmax()
      return_minmax.thisptr = ret_minmax
      return return_minmax
   def fasta_count_kmers_by_position(self, char* inputfile, unsigned int max_read_len, 
                                     BoundedCounterType limit_by_count=0, 
                                     callback_obj=None):
      global _callback_obj
      if callback_obj is not None:
         _callback_obj = callback_obj

      positions = []
      cpositions = self.thisptr.fasta_count_kmers_by_position(inputfile, max_read_len,
                                                              limit_by_count,
                                                              _report_fn, <void*>_callback_obj)
      for i in range(max_read_len):
         positions.append(cpositions[i])
      return positions
   def consume_fasta(self, char* filename, HashIntoType lower_bound=0,
                     HashIntoType upper_bound=0, callback_obj=None):
      global _callback_obj
      cdef unsigned long long n_consumed = 0
      cdef unsigned int total_reads = 0
      if callback_obj is not None:
         _callback_obj = callback_obj
      self.thisptr.consume_fasta(filename, total_reads, n_consumed, lower_bound,
                                 upper_bound, _report_fn, <void*>_callback_obj)
      return total_reads, n_consumed
   def fasta_dump_kmers_by_abundance(self, char* inputfile, 
                                     BoundedCounterType limit_by_count, callback_obj=None):
      global _callback_obj
      if callback_obj is not None:
         _callback_obj = callback_obj
      self.thisptr.fasta_dump_kmers_by_abundance(inputfile, limit_by_count,
                                                 _report_fn, <void*>_callback_obj)
   def set_use_bigcount(self, bool b):
      self.thisptr.set_use_bigcount(b)
   def get_use_bigcount(self):
      return self.thisptr.get_use_bigcount()
   def max_hamming1_count(self, char* kmer):
      return self.thisptr.max_hamming1_count(kmer)
   def trim_on_abundance(self, seq, unsigned int min_count):
      cdef unsigned int trim_at = self.thisptr.trim_on_abundance(seq, min_count)
      return seq[0:trim_at], trim_at
   def get_kmer_abund_abs_deviation(self, char* filename, float mean):
      cdef float abs_dev = 0.0
      self.thisptr.get_kmer_abund_abs_deviation(filename, mean, abs_dev)
      return abs_dev
   def get_kmer_abund_mean(self, char* filename):
      cdef unsigned long long total = 0
      cdef unsigned long long count = 0
      cdef float mean = 0.0
      self.thisptr.get_kmer_abund_mean(filename, total, count, mean)
      return total, count, mean
   def get_median_count(self, s):
      cdef BoundedCounterType med = 0
      cdef float average = 0
      cdef float stddev = 0
      self.thisptr.get_median_count(s, med, average, stddev)
      return med, average, stddev
   def get_kadian_count(self, s, nk=1):
      cdef BoundedCounterType kad = 0
      self.thisptr.get_kadian_count(s, kad, nk)
      return kad

def forward_hash(char* s, WordLength k):
   return _hash(s, k)
def forward_hash_no_rc(char* s, WordLength k):
   return _hash_forward(s, k)
def reverse_hash(HashIntoType h, WordLength k):
   cdef string s = _revhash(h, k)
   cdef char *cstr = <char *>malloc(s.length() + 1)
   strcpy(cstr, s.c_str())
   return cstr

cdef class _new_subsetpartition

cdef class _new_hashbits:
   cdef Hashbits *thisptr
   def __cinit__(self, WordLength k, list primes):
      cdef vector[HashIntoType] v = vector[HashIntoType]()
      for prime in primes:
         v.push_back(prime) 
      self.thisptr = new Hashbits(k, v)
   def n_occupied(self):
      return self.thisptr.n_occupied(0, 0)
   def n_unique_kmers(self, HashIntoType a=0, HashIntoType b=0):
      return self.thisptr.n_kmers(a, b)
   def ksize(self):
      return self.thisptr.ksize()
   def count(self, char* s):
      self.thisptr.count(s)
   def count(self, HashIntoType n):
      self.thisptr.count(n)
   def get(self, arg):
      if isinstance(arg, str):
         return self.thisptr.get_count(<char*>arg)
      else:
         return self.thisptr.get_count(<HashIntoType>arg)
   def __getitem__(self, arg):
      if isinstance(arg, str):
         return self.thisptr.get_count(<char*>arg)
      else:
         return self.thisptr.get_count(<HashIntoType>arg)
   def n_tags(self):
      return self.thisptr.n_tags()
   def save(self, char* s):
      self.thisptr.save(s)
   def load(self, char* s):
      self.thisptr.load(s)
   def save_tagset(self, char* s):
      self.thisptr.save_tagset(s)
   def load_tagset(self, char* s, bool clear_tags=1):
      self.thisptr.load_tagset(s, clear_tags)
   def count_kmers_within_radius(self, char* kmer,
                                 unsigned int n1=0, unsigned int n2=0):
      cdef HashIntoType h1 = 0
      cdef HashIntoType h2 = 0
      _hash(kmer, self.thisptr.ksize(), h1, h2)
      return self.thisptr.count_kmers_within_radius(h1, h2, n1, n2)
   def calc_connected_graph_size(self, char * kmer, 
                                 unsigned long long max_size=0,
                                 bool break_on_circum=0):
      cdef set[HashIntoType] *keeper = new set[HashIntoType]()
      cdef unsigned long long size = 0
      self.thisptr.calc_connected_graph_size(kmer, size, keeper[0], 
                                             max_size, break_on_circum)
      return size
   def consume_partitioned_fasta(self, char* filename, callback_obj=None):
      global _callback_obj
      if callback_obj is not None:
         _callback_obj = callback_obj
      cdef unsigned int total_reads = 0
      cdef unsigned long long n_consumed = 0
      self.thisptr.consume_partitioned_fasta(filename, total_reads, 
                                             n_consumed, _report_fn, <void*>_callback_obj)

      return total_reads, n_consumed
   def consume_fasta_and_traverse(self, char* filename, unsigned int radius,
                                  unsigned int big_threshold, unsigned int transfer_threshold,
                                  _new_counting_hash ch):
      self.thisptr.consume_fasta_and_traverse(filename, radius, big_threshold, transfer_threshold,
                                              ch.thisptr[0])
   def traverse_from_reads(self, char* filename, unsigned int radius, 
                           unsigned int big_threshold, unsigned int transfer_threshold, 
                           _new_counting_hash ch):
      self.thisptr.traverse_from_reads(filename, radius, big_threshold, transfer_threshold,
                           ch.thisptr[0])
   def consume_fasta_and_tag_with_stoptags(self, char* filename, callback_obj=None):
      global _callback_obj
      if callback_obj is not None:
         _callback_obj = callback_obj
      cdef unsigned int total_reads = 0
      cdef unsigned long long n_consumed = 0
      self.thisptr.consume_fasta_and_tag_with_stoptags(filename, total_reads, n_consumed,
                                                       _report_fn, <void*>_callback_obj)
      return total_reads, n_consumed
   def consume_fasta_and_tag(self, char* filename, callback_obj=None):
      global _callback_obj
      if callback_obj is not None:
         _callback_obj = callback_obj
      cdef unsigned int total_reads = 0
      cdef unsigned long long n_consumed = 0
      self.thisptr.consume_fasta_and_tag(filename, total_reads, n_consumed, 
                                         _report_fn, <void*>_callback_obj)
      return total_reads, n_consumed
   def consume(self, char* s, HashIntoType lower = 0, HashIntoType upper = 0):
      return self.thisptr.consume_string(s, lower, upper)
   def consume_fasta(self, char* filename, HashIntoType lower=0, 
                     HashIntoType upper=0, callback_obj=None):
      global _callback_obj
      cdef unsigned long long n_consumed  = 0
      cdef unsigned int total_reads       = 0

      if callback_obj is not None:
         _callback_obj = callback_obj

      self.thisptr.consume_fasta(filename, total_reads, n_consumed, lower, upper, 
                                 _report_fn, <void*>callback_obj)
      return n_consumed, total_reads
   def do_subset_partition(self, HashIntoType start, HashIntoType end, bool break_stop=0, bool stop_big_traversal=0,
                           callback_obj=None):
      subset = _new_subsetpartition(self)
      subset.do_partition(start, end, break_stop, stop_big_traversal, callback_obj)
      return subset
   def save_partitionmap(self, char* filename):
      self.thisptr.partition.save_partitionmap(filename)
   def save_subset_partitionmap(self, subset, char* filename):
      subset.save_partitionmap(filename)
   def load_subset_partitionmap(self, filename):
      cdef SubsetPartition * subset_p = new SubsetPartition(self.thisptr)
      subset_p.load_partitionmap(filename)
      subset = _new_subsetpartition()
      subset.thisptr = subset_p
      return subset
   def load_partitionmap(self, filename):
      self.thisptr.partition.load_partitionmap(filename)
   def _validate_subset_partitionmap(self, subset):
      subset._validate_pmap()
   def add_tag(self, char* s):
      cdef HashIntoType kmer = _hash(s, self.thisptr.ksize())
      self.thisptr.add_tag(kmer)
   def add_stop_tag(self, char* s):
      cdef HashIntoType kmer = _hash(s, self.thisptr.ksize())
      self.thisptr.add_stop_tag(kmer)
   def load_stop_tags(self, char* filename):
      self.thisptr.load_stop_tags(filename)
   def save_stop_tags(self, char* filename):
      self.thisptr.save_stop_tags(filename)
   def print_stop_tags(self, char* filename):
      self.thisptr.print_stop_tags(filename)
   def output_partitions(self, infile, outfile, bool output_unassigned=0, callback_obj=None):
      global _callback_obj
      if callback_obj is not None:
         _callback_obj = callback_obj
      return self.thisptr.partition.output_partitioned_file(infile, outfile,
                                                            output_unassigned, 
                                                            _report_fn, <void*>_callback_obj)
   def find_unpart(self, infile, bool traverse, bool stop_big_traversals, callback_obj=None):
      global _callback_obj
      if callback_obj is not None:
         _callback_obj = callback_obj
      return self.thisptr.partition.find_unpart(infile, traverse, stop_big_traversals, _report_fn, <void*>_callback_obj)
   def count_partitions(self):
      cdef unsigned int n_partitions = 0
      cdef unsigned int n_unassigned = 0
      self.thisptr.partition.count_partitions(n_partitions, n_unassigned)
      return n_partitions, n_unassigned
   def subset_count_partitions(self, _new_subsetpartition subset_obj):
      cdef unsigned int n_partitions = 0
      cdef unsigned int n_unassigned = 0
      subset_obj.thisptr.count_partitions(n_partitions, n_unassigned)
      return n_partitions, n_unassigned
   def merge_subset_from_disk(self, char* filename):
      self.thisptr.partition.merge_from_disk(filename)
   def merge_subset(self, _new_subsetpartition subset_obj):
      self.thisptr.partition.merge(subset_obj.thisptr)
   def merge2_subset(self, _new_subsetpartition s1, _new_subsetpartition s2):
      s1.thisptr.merge(s2.thisptr)
   def merge2_subset_from_disk(self, subset_obj, char* filename):
      subset_obj.thisptr.merge_from_disk(filename)
   def join_partitions_by_path(self, char* sequence):
      self.thisptr.partition.join_partitions_by_path(sequence)
   def divide_tags_into_subsets(self, unsigned int subset_size):
      subsets = []
      cdef set[HashIntoType] divvy
      self.thisptr.divide_tags_into_subsets(subset_size, divvy)
      cdef set[HashIntoType].iterator it = divvy.begin()
      while it != divvy.end():
         subsets.append(deref(it)) 
         inc(it)
      return subsets
   def kmer_degree(self, char* s):
      return self.thisptr.kmer_degree(s)
   def set_partition_id(self, char* kmer, PartitionID p):
      self.thisptr.partition.set_partition_id(kmer, p)
   def join_partitions(self, PartitionID p1, PartitionID p2):
      self.thisptr.partition.join_partitions(p1, p2)
   def get_partition_id(self, char* s):
      return self.thisptr.partition.get_partition_id(s)
   def is_single_partition(self, char* s):
      return self.thisptr.partition.is_single_partition(s) 
   def _set_tag_density(self, unsigned int d):
      self.thisptr._set_tag_density(d)
   def _get_tag_density(self):
      return self.thisptr._get_tag_density()
   def trim_on_degree(self, seq, unsigned int max_degree):
      cdef unsigned int trim_at = self.thisptr.trim_on_degree(seq, max_degree)
      trim_seq = seq[0:trim_at]
      return trim_seq, trim_at
   def trim_on_sodd(self, seq, unsigned int max_sodd):
      cdef unsigned int trim_at = self.thisptr.trim_on_sodd(seq, max_sodd)
      trim_seq = seq[0:trim_at]
      return trim_seq, trim_at
   def trim_on_stoptags(self, seq):
      cdef unsigned int trim_at = self.thisptr.trim_on_stoptags(seq)
      trim_seq = seq[0:trim_at]
      return trim_seq, trim_at
   def trim_on_density_explosion(self, seq, unsigned long radius, unsigned long max_volume):
      cdef unsigned int trim_at = self.thisptr.trim_on_density_explosion(seq, radius, max_volume)
      trim_seq = seq[0:trim_at]
      return trim_seq, trim_at
   def identify_stoptags_by_position(self, seq):
      cdef vector[unsigned int] posns
      self.thisptr.identify_stop_tags_by_position(seq, posns)
      poslist = []
      cdef vector[unsigned int].iterator it = posns.begin()
      while it != posns.end():
         poslist.append(deref(it))
         inc(it)
      return poslist
   def hashsizes(self):
      cdef vector[HashIntoType] res = self.thisptr.get_tablesizes()
      hashes = []
      
      cdef vector[HashIntoType].iterator it = res.begin()
      while it != res.end():
         hashes.append(deref(it))
         inc(it)
      return hashes
   def extract_unique_paths(self, seq, unsigned int min_length, float min_unique_f):
      cdef vector[string] results
      cdef char* cstr
      cdef string s
      self.thisptr.extract_unique_paths(seq, min_length, min_unique_f, results)
      res = []
      cdef vector[string].iterator it = <vector[string].iterator> results.begin()
      while it != results.end():
         s = deref(it)
         cstr = <char*>malloc(s.length()+1)
         strcpy(cstr, s.c_str())
         res.append(cstr)
         inc(it)
      return res
   def filter_if_present(self, char* infile, char* output, callback_obj=None):
      global _callback_obj

      if callback_obj is not None:
         _callback_obj = callback_obj

      self.thisptr.filter_if_present(infile, output, _report_fn, <void*>_callback_obj)
   def _validate_partitionmap(self):
      self.thisptr.partition._validate_pmap()
   def hitraverse_to_stoptags(self, char* filename, _new_counting_hash ch, unsigned int cutoff):
      self.thisptr.hitraverse_to_stoptags(filename, ch.thisptr[0], cutoff)
   def find_radius_for_volume(self, char* kmer, unsigned long max_count, unsigned long max_radius):
      cdef HashIntoType kmer_f = 0, kmer_r = 0
      _hash(kmer, self.ksize(), kmer_f, kmer_r)
      return self.thisptr.find_radius_for_volume(kmer_f, kmer_r, max_count, max_radius)
   def count_kmers_on_radius(self, char* kmer, unsigned int radius, unsigned int max_volume):
      cdef HashIntoType kmer_f = 0, kmer_r = 0
      _hash(kmer, self.ksize(), kmer_f, kmer_r)
      return self.thisptr.count_kmers_on_radius(kmer_f, kmer_r, radius, max_volume)
   def repartition_largest_partition(self, _new_subsetpartition subset_o, 
                                     _new_counting_hash counting_o, 
                                     unsigned int distance, 
                                     unsigned int threshold, unsigned int frequency):
      if subset_o is None:
         subset_o = _new_subsetpartition()
         subset_o.thisptr = self.thisptr.partition
      return subset_o.thisptr.repartition_largest_partition(distance, threshold, frequency, 
                                                     counting_o.thisptr[0])
   def subset_partition_size_distribution(self, _new_subsetpartition subset_o):
      cdef map[unsigned long long, unsigned long long] *d = new map[unsigned long long, unsigned long long]()
      cdef unsigned int n_unassigned = 0 
      subset_o.thisptr.partition_size_distribution(d[0], n_unassigned)
     
      cdef pair[unsigned long long, unsigned long long] mypair
      cdef map[unsigned long long, unsigned long long].iterator it = d.begin()
      cdef unsigned int tempOne = 0
      cdef unsigned int tempTwo = 0
      dist = []
      while it != d.end():
         tempOne = deref(it).first
         tempTwo = deref(it).second
         dist.append((tempOne, tempTwo))
         inc(it)

      return dist, n_unassigned
   def find_all_tags(self, char* kmer_s):
      cdef HashIntoType kmer
      cdef HashIntoType kmer_f = 0
      cdef HashIntoType kmer_r = 0
      kmer = _hash(kmer_s, self.thisptr.ksize(), kmer_f, kmer_r)
      ppi = _pre_partition_info(kmer)
      self.thisptr.partition.find_all_tags(kmer_f, kmer_r, ppi.thisptr.tagged_kmers, 
                                           self.thisptr.all_tags)

      self.thisptr.add_kmer_to_tags(kmer)

      return ppi
   def assign_partition_id(self, _pre_partition_info ppi_obj):
      cdef PartitionID p = self.thisptr.partition.assign_partition_id(ppi_obj.thisptr.kmer, 
                                                                      ppi_obj.thisptr.tagged_kmers)
      return p

cdef class _new_subsetpartition:
   cdef SubsetPartition *thisptr
   def __cinit__(self, arg=None):
      if arg is not None:
         hb = <_new_hashbits> arg
         self.thisptr = new SubsetPartition(&hb.thisptr[0])
   def merge(self, _new_subsetpartition sp):
      self.thisptr.merge(sp.thisptr)
   def merge_from_disk(self, char* s):
      self.thisptr.merge_from_disk(s)
   def save_partitionmap(self, char* s):
      self.thisptr.save_partitionmap(s)
   def load_partitionmap(self, char* s):
      self.thisptr.load_partitionmap(s)
   def _validate_pmap(self):
      self.thisptr._validate_pmap()

   def do_partition(self, HashIntoType first_kmer, HashIntoType last_kmer, bool break_stop, bool stop_big_traversals, callback_obj=None):
      global _callback_obj
      if callback_obj is not None:
         _callback_obj = callback_obj
      self.thisptr.do_partition(first_kmer, last_kmer, break_stop, stop_big_traversals, _report_fn, <void*>_callback_obj)

   def count_partitions(self):
      cdef unsigned int n_partitions = 0
      cdef unsigned int n_unassigned = 0
      self.thisptr.count_partitions(n_partitions, n_unassigned)
      return n_partitions, n_unassigned

   def output_partitioned_file(self, infile, outfile, 
                               bool output_unassigned=0,
                               callback_obj=None):
      global _callback_obj
      if callback_obj is not None:
         _callback_obj = callback_obj
      return self.thisptr.output_partitioned_file(infile, outfile,
                                                  output_unassigned,
                                                  _report_fn, <void*>_callback_obj)

# vim: set sts=3 sw=3 et:
