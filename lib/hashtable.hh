#ifndef HASHTABLE_HH
#define HASHTABLE_HH

#if (__cplusplus >= 201103L)
#   include <cstdint>
#else
extern "C"
{
#   include <stdint.h>
}
#endif

#include <vector>
#include <iostream>
#include <list>
#include <queue>

#include <fstream>
#include <string>
#include <set>
#include <map>
#include <queue>

#include "khmer.hh"
#include "storage.hh"
#include "read_parsers.hh"

#define CALLBACK_PERIOD 100000

namespace khmer {

  typedef unsigned int PartitionID;
  typedef std::set<HashIntoType> SeenSet;
  typedef std::set<PartitionID> PartitionSet;
  typedef std::map<HashIntoType, PartitionID*> PartitionMap;
  typedef std::map<PartitionID, PartitionID*> PartitionPtrMap;
  typedef std::map<PartitionID, SeenSet*> PartitionsToTagsMap;
  typedef std::set<PartitionID *> PartitionPtrSet;
  typedef std::map<PartitionID, PartitionPtrSet*> ReversePartitionMap;
  typedef std::queue<HashIntoType> NodeQueue;
  typedef std::map<PartitionID, PartitionID*> PartitionToPartitionPMap;
  typedef std::map<HashIntoType, unsigned int> TagCountMap;
  typedef std::map<PartitionID, unsigned int> PartitionCountMap;
  typedef std::map<unsigned long long, unsigned long long> PartitionCountDistribution;

  struct HashTableReadConsumerPerformanceMetrics : public IPerformanceMetrics
  {
	
	enum
	{
	    // TODO: Declare.
	};

	// TODO: Declare.

  };

  //
  // Sequence iterator class, test.  Not really a C++ iterator yet.
  //

  class KMerIterator {
  protected:
    const char * _seq;
    const unsigned char _ksize;
    
    HashIntoType _kmer_f, _kmer_r;
    HashIntoType bitmask;
    unsigned int _nbits_sub_1;
    unsigned int index, length;
    bool initialized;
  public:
    KMerIterator(const char * seq, unsigned char k) : _seq(seq), _ksize(k) {
      bitmask = 0;
      for (unsigned int i = 0; i < _ksize; i++) {
	bitmask = (bitmask << 2) | 3;
      }
      _nbits_sub_1 = (_ksize*2 - 2);

      index = _ksize - 1;
      length = strlen(seq);
      initialized = false;
    }

    HashIntoType first(HashIntoType& f, HashIntoType& r) {
      HashIntoType x;
      x = _hash(_seq, _ksize, _kmer_f, _kmer_r);

      f = _kmer_f;
      r = _kmer_r;

      index = _ksize;

      return x;
    }

    HashIntoType next(HashIntoType& f, HashIntoType& r) {
      if (done()) {
	throw std::exception();
      }

      if (!initialized) {
	initialized = true;
	return first(f, r);
      }

      unsigned char ch = _seq[index];
      index++;
      assert(index <= length);

      // left-shift the previous hash over
      _kmer_f = _kmer_f << 2;

      // 'or' in the current nt
      _kmer_f |= twobit_repr(ch);

      // mask off the 2 bits we shifted over.
      _kmer_f &= bitmask;

      // now handle reverse complement
      _kmer_r = _kmer_r >> 2;
      _kmer_r |= (twobit_comp(ch) << _nbits_sub_1);

      f = _kmer_f;
      r = _kmer_r;

      return uniqify_rc(_kmer_f, _kmer_r);
    }

    HashIntoType first() { return first(_kmer_f, _kmer_r); }
    HashIntoType next() { return next(_kmer_f, _kmer_r); }

    bool done() { return index >= length; }
  }; // class KMerIterator


  class Hashtable {		// Base class implementation of a Bloom ht.

  protected:

    struct Hasher
    {

	uint32_t			thread_id;
	// TODO: Add HasherPerformanceMetrics instance.
	TraceLogger			trace_logger;

	Hasher(
	    uint32_t const  thread_id,
	    uint8_t const   trace_level = TraceLogger:: TLVL_NONE
	);
	~Hasher( );

    }; // struct Hasher
    
    uint8_t	    _trace_level;

    uint32_t	    _number_of_threads;
    ThreadIDMap	    _thread_id_map;
    Hasher **	    _hashers;
    unsigned int    _max_count;
    unsigned int    _max_bigcount;

    WordLength	    _ksize;
    HashIntoType    bitmask;
    unsigned int    _nbits_sub_1;

    Hashtable(
	WordLength	ksize,
	uint32_t const	number_of_threads   = 
	get_active_config( ).get_number_of_threads( ),
	uint8_t const	trace_level	    = TraceLogger:: TLVL_NONE
    )
    :	_trace_level( trace_level ),
	_number_of_threads( number_of_threads ), 
	_thread_id_map( ThreadIDMap( number_of_threads ) ),
	_ksize( ksize )
    {

	_hashers = new Hasher *[ number_of_threads ];
	for (uint32_t i = 0; i < number_of_threads; ++i) _hashers[ i ] = NULL;

	// TODO: Transplant this logic to someplace more common.
	if (1 == number_of_threads)
	{
	    _max_count	    = MAX_COUNT;
	    _max_bigcount   = MAX_BIGCOUNT;
	}
	else
	{
	    _max_count	    = MAX_COUNT - number_of_threads + 1;
	    _max_bigcount   = MAX_BIGCOUNT - number_of_threads + 1;
	}

	_init_bitstuff();

    }

    virtual ~Hashtable( )
    {
	
	for (uint32_t i = 0; i < _number_of_threads; ++i)
	{
	    if (NULL != _hashers[ i ])
	    {
		delete _hashers[ i ];
		_hashers[ i ]	= NULL;
	    }
	}
	delete [ ] _hashers;
	_hashers    = NULL;

    }

    void _init_bitstuff() {
      bitmask = 0;
      for (unsigned int i = 0; i < _ksize; i++) {
	bitmask = (bitmask << 2) | 3;
      }
      _nbits_sub_1 = (_ksize*2 - 2);
    }


    inline Hasher   &_get_hasher( )
    {
	uint32_t	thread_id	= _thread_id_map.get_thread_id( );
	Hasher *	hasher_PTR	= NULL;

	assert( NULL != _hashers );

	hasher_PTR = _hashers[ thread_id ];
	if (NULL == hasher_PTR)
	{
	    _hashers[ thread_id ]   =
	    new Hasher( thread_id, _trace_level );
	    hasher_PTR		    = _hashers[ thread_id ];
	}

	return *hasher_PTR;
    }


    HashIntoType _next_hash(char ch, HashIntoType &h, HashIntoType &r) const {
      // left-shift the previous hash over
      h = h << 2;

      // 'or' in the current nt
      h |= twobit_repr(ch);

      // mask off the 2 bits we shifted over.
      h &= bitmask;

      // now handle reverse complement
      r = r >> 2;
      r |= (twobit_comp(ch) << _nbits_sub_1);

      return uniqify_rc(h, r);
    }

  public:

    // accessor to get 'k'
    const WordLength ksize() const { return _ksize; }

    virtual void count(const char * kmer) = 0;
    virtual void count(HashIntoType khash) = 0;

    // get the count for the given k-mer.
    virtual const BoundedCounterType get_count(const char * kmer) const = 0;
    virtual const BoundedCounterType get_count(HashIntoType khash) const = 0;

    virtual void save(std::string) = 0;
    virtual void load(std::string) = 0;

    // count every k-mer in the string.
    unsigned int consume_string(const std::string &s,
				HashIntoType lower_bound = 0,
				HashIntoType upper_bound = 0);
    
    // checks each read for non-ACGT characters
    bool check_and_normalize_read(std::string &read) const;

    // check each read for non-ACGT characters, and then consume it.
    unsigned int check_and_process_read(std::string &read,
					bool &is_valid,
					HashIntoType lower_bound = 0,
					HashIntoType upper_bound = 0);
    
    // Count every k-mer in a FASTA or FASTQ file.
    // Note: Yes, the name 'comsume_fasta' is a bit misleading, 
    //	     but the FASTA format is effectively a subset of the FASTQ format
    //	     and the FASTA portion is what we care about in this case.
    void consume_fasta(
	std::string const   &filename,
	unsigned int	    &total_reads,
	unsigned long long  &n_consumed,
	HashIntoType	    lower_bound	    = 0,
	HashIntoType	    upper_bound	    = 0,
	CallbackFn	    callback	    = NULL,
	void *		    callback_data   = NULL
    );
    // Count every k-mer from a stream of FASTA or FASTQ reads, 
    // using the supplied parser.
    void consume_fasta(
	read_parsers:: IParser *	    parser,
	unsigned int	    &total_reads,
	unsigned long long  &n_consumed,
	HashIntoType	    lower_bound	    = 0,
	HashIntoType	    upper_bound	    = 0,
	CallbackFn	    callback	    = NULL,
	void *		    callback_data   = NULL
    );

  };
			  

};

#endif // HASHTABLE_HH
