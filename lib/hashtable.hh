#ifndef HASHTABLE_HH
#define HASHTABLE_HH

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

  struct HashTablePerformanceMetrics : public IPerformanceMetrics
  {
	
	enum
	{
	    MKEY_TIME_NORM_READ,
	    MKEY_TIME_HASH_KMER,
	    MKEY_TIME_UPDATE_TALLIES
	};

	uint64_t	clock_nsecs_norm_read;
	uint64_t	cpu_nsecs_norm_read;
	uint64_t	clock_nsecs_hash_kmer;
	uint64_t	cpu_nsecs_hash_kmer;
	uint64_t	clock_nsecs_update_tallies;
	uint64_t	cpu_nsecs_update_tallies;

		HashTablePerformanceMetrics( );
	virtual ~HashTablePerformanceMetrics( );

	virtual void	accumulate_timer_deltas( uint32_t metrics_key );

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

	uint32_t			pool_id;
	uint32_t			thread_id;
	HashTablePerformanceMetrics	pmetrics;
	TraceLogger			trace_logger;

	Hasher(
	    uint32_t const  pool_id,
	    uint32_t const  thread_id,
	    uint8_t const   trace_level = TraceLogger:: TLVL_NONE
	);
	~Hasher( );

    }; // struct Hasher
    
    uint8_t	    _trace_level;

    uint32_t	    _number_of_threads;
    uint32_t	    _tpool_map_spin_lock;
    uint32_t	    _thread_pool_counter;
    std:: map< int, uint32_t >
		    _thread_pool_id_map;
    std:: map< uint32_t, ThreadIDMap * >
		    _thread_id_maps;
    std:: map< uint32_t, Hasher ** >
		    _hashers_map;
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
	_tpool_map_spin_lock( 0 ),
	_thread_pool_counter( 0 ),
	_max_count( MAX_COUNT - number_of_threads + 1 ),
	_max_bigcount( MAX_BIGCOUNT - number_of_threads + 1 ),
	_ksize( ksize )
    { _init_bitstuff( ); }

    virtual ~Hashtable( )
    {
	std:: map< int, uint32_t >:: iterator it;
	uint32_t thread_pool_id;
	Hasher ** hashers = NULL;

	for (it = _thread_pool_id_map.begin( );
	     it != _thread_pool_id_map.end( );
	     ++it)
	{
	    thread_pool_id = it->second;

	    delete _thread_id_maps[ thread_pool_id ];
	    _thread_id_maps[ thread_pool_id ] = NULL;

	    hashers = _hashers_map[ thread_pool_id ];
	    for (uint32_t i = 0; i < _number_of_threads; ++i)
	    {
		if (NULL != hashers[ i ])
		{
		    delete hashers[ i ];
		    hashers[ i ] = NULL;
		}
	    }
	    delete [ ] hashers;
	    _hashers_map[ thread_pool_id ] = NULL;
	}
    }

    void _init_bitstuff() {
      bitmask = 0;
      for (unsigned int i = 0; i < _ksize; i++) {
	bitmask = (bitmask << 2) | 3;
      }
      _nbits_sub_1 = (_ksize*2 - 2);
    }


    inline Hasher   &_get_hasher( int uuid = 0 )
    {
	std:: map< int, uint32_t >:: iterator	match;	
	uint32_t				thread_pool_id;
	ThreadIDMap *				thread_id_map	= NULL;
	uint32_t				thread_id;
	Hasher **				hashers		= NULL;
	Hasher *				hasher_PTR	= NULL;
	
	match = _thread_pool_id_map.find( uuid );
	if (match == _thread_pool_id_map.end( ))
	{
	    
	    while (!__sync_bool_compare_and_swap(
		&_tpool_map_spin_lock, 0, 1
	    ));

	    // TODO: Handle 'std:: bad_alloc' exceptions.
	    thread_pool_id			= _thread_pool_counter++;
	    _thread_pool_id_map[ uuid ]		= thread_pool_id;
	    _thread_id_maps[ thread_pool_id ]	=
	    new ThreadIDMap( _number_of_threads );
	    _hashers_map[ thread_pool_id ]	=
	    new Hasher *[ _number_of_threads ];
	    hashers				=
	    _hashers_map[ thread_pool_id ];
	    for (uint32_t i = 0; i < _number_of_threads; ++i)
		hashers[ i ] = NULL;
	    
	    __sync_bool_compare_and_swap( &_tpool_map_spin_lock, 1, 0 );

	    match = _thread_pool_id_map.find( uuid );
	} // no thread pool for UUID

	thread_pool_id	    = (*match).second;
	thread_id_map	    = _thread_id_maps[ thread_pool_id ];
	thread_id	    = thread_id_map->get_thread_id( );
	hashers		    = _hashers_map[ thread_pool_id ];
	hasher_PTR	    = hashers[ thread_id ];
	if (NULL == hasher_PTR)
	{
	    hashers[ thread_id ]    =
	    new Hasher( thread_pool_id, thread_id, _trace_level );
	    hasher_PTR		    = hashers[ thread_id ];
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
			  
    void get_median_count(const std::string &s,
			  BoundedCounterType &median,
			  float &average,
			  float &stddev);

  };
};

#endif // HASHTABLE_HH
