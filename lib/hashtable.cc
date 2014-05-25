//
// This file is part of khmer, http://github.com/ged-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2013. It is licensed under
// the three-clause BSD license; see doc/LICENSE.txt.
// Contact: khmer-project@idyll.org
//

#include "khmer.hh"
#include "hashtable.hh"
#include "read_parsers.hh"
#include "counting.hh"

#include <algorithm>

using namespace std;
using namespace khmer;
using namespace khmer:: read_parsers;

#ifdef WITH_INTERNAL_METRICS
HashTablePerformanceMetrics::
HashTablePerformanceMetrics( )
    : IPerformanceMetrics( ),
      clock_nsecs_norm_read( 0 ),
      cpu_nsecs_norm_read( 0 ),
      clock_nsecs_hash_kmer( 0 ),
      cpu_nsecs_hash_kmer( 0 ),
      clock_nsecs_update_tallies( 0 ),
      cpu_nsecs_update_tallies( 0 )
{ }


HashTablePerformanceMetrics::
~HashTablePerformanceMetrics( )
{ }


void
HashTablePerformanceMetrics::
accumulate_timer_deltas( uint32_t metrics_key )
{

    switch (metrics_key) {
    case MKEY_TIME_NORM_READ:
        clock_nsecs_norm_read +=
            _timespec_diff_in_nsecs( _temp_clock_start, _temp_clock_stop );
        cpu_nsecs_norm_read   +=
            _timespec_diff_in_nsecs( _temp_cpu_start, _temp_cpu_stop );
        break;
    case MKEY_TIME_HASH_KMER:
        clock_nsecs_hash_kmer +=
            _timespec_diff_in_nsecs( _temp_clock_start, _temp_clock_stop );
        cpu_nsecs_hash_kmer   +=
            _timespec_diff_in_nsecs( _temp_cpu_start, _temp_cpu_stop );
        break;
    case MKEY_TIME_UPDATE_TALLIES:
        clock_nsecs_update_tallies +=
            _timespec_diff_in_nsecs( _temp_clock_start, _temp_clock_stop );
        cpu_nsecs_update_tallies   +=
            _timespec_diff_in_nsecs( _temp_cpu_start, _temp_cpu_stop );
        break;
    default:
        throw InvalidPerformanceMetricsKey( );
    }

}
#endif

Hashtable:: Hasher::
Hasher(
    uint32_t const  pool_id,
    uint32_t const  thread_id,
    uint8_t const	  trace_level
)
    : pool_id( pool_id ),
      thread_id( thread_id ),
#ifdef WITH_INTERNAL_METRICS
      pmetrics( HashTablePerformanceMetrics( ) ),
#endif
      trace_logger(
          TraceLogger(
              trace_level, "hashtable-%lu-%lu.log",
              (unsigned long int)pool_id, (unsigned long int)thread_id
          )
      )
{ }


Hashtable:: Hasher::
~Hasher( )
{ }


//
// check_and_process_read: checks for non-ACGT characters before consuming
//

unsigned int Hashtable::check_and_process_read(std::string &read,
        bool &is_valid)
{
    is_valid = check_and_normalize_read(read);

    if (!is_valid) {
        return 0;
    }

    return consume_string(read);
}

//
// check_and_normalize_read: checks for non-ACGT characters
//			     converts lowercase characters to uppercase one
// Note: Usually it is desirable to keep checks and mutations separate.
//	 However, in the interests of efficiency (we are potentially working
//	 with TB of data), a check and mutation have been placed inside the
//	 same loop. Potentially trillions fewer fetches from memory would
//	 seem to be a worthwhile goal.
//

bool Hashtable::check_and_normalize_read(std::string &read) const
{
    bool rc = true;
#if (0)  // TODO: WITH_INTERNAL_TRACING < some_threshold
    Hasher		  &hasher		= _get_hasher( );
#endif

    if (read.length() < _ksize) {
        return false;
    }

#if (0)   // TODO: WITH_INTERNAL_TRACING < some_threshold
    hasher.pmetrics.start_timers( );
#endif
    for (unsigned int i = 0; i < read.length(); i++)  {
        read[ i ] &= 0xdf; // toupper - knock out the "lowercase bit"
        if (!is_valid_dna( read[ i ] )) {
            rc = false;
            break;
        }
    }
#if (0)  // TODO: WITH_INTERNAL_TRACING < some_threshold
    hasher.pmetrics.stop_timers( );
    hasher.pmetrics.accumulate_timer_deltas(
        (uint32_t)HashTablePerformanceMetrics:: MKEY_TIME_NORM_READ
    );
#endif

    return rc;
}

//
// consume_fasta: consume a FASTA file of reads
//

// TODO? Inline in header.
void
Hashtable::
consume_fasta(
    std:: string const  &filename,
    unsigned int	      &total_reads, unsigned long long	&n_consumed,
    CallbackFn	      callback,	    void *		callback_data
)
{
    khmer:: Config    &the_config	  = khmer:: get_active_config( );

    // Note: Always assume only 1 thread if invoked this way.
    IParser *	  parser =
        IParser::get_parser(
            filename, 1, the_config.get_reads_input_buffer_size( ),
            the_config.get_reads_parser_trace_level( )
        );

    consume_fasta(
        parser,
        total_reads, n_consumed,
        callback, callback_data
    );

    delete parser;
}

void
Hashtable::
consume_fasta(
    read_parsers:: IParser *  parser,
    unsigned int		    &total_reads, unsigned long long  &n_consumed,
    CallbackFn		    callback,	  void *	      callback_data
)
{
    Hasher		  &hasher		=
        _get_hasher( parser->uuid( ) );
#if (0) // Note: Used with callback - currently disabled.
    unsigned long long int  n_consumed_LOCAL	= 0;
#endif
    Read			  read;

    hasher.trace_logger(
        TraceLogger:: TLVL_DEBUG2, "Starting trace of 'consume_fasta'....\n"
    );

    // Iterate through the reads and consume their k-mers.
    while (!parser->is_complete( )) {
        unsigned int  this_n_consumed;
        bool	  is_valid;

        read = parser->get_next_read( );

        this_n_consumed =
            check_and_process_read(read.sequence, is_valid);

#ifdef WITH_INTERNAL_METRICS
        hasher.pmetrics.start_timers( );
#endif
#if (0) // Note: Used with callback - currently disabled.
        n_consumed_LOCAL  = __sync_add_and_fetch( &n_consumed, this_n_consumed );
#else
        __sync_add_and_fetch( &n_consumed, this_n_consumed );
#endif
        unsigned int total_reads_LOCAL = __sync_add_and_fetch( &total_reads, 1 );
#ifdef WITH_INTERNAL_METRICS
        hasher.pmetrics.stop_timers( );
        hasher.pmetrics.accumulate_timer_deltas(
            (uint32_t)HashTablePerformanceMetrics:: MKEY_TIME_UPDATE_TALLIES
        );
#endif

        if (0 == (total_reads_LOCAL % 10000))
            hasher.trace_logger(
                TraceLogger:: TLVL_DEBUG3,
                "Total number of reads processed: %llu\n",
                (unsigned long long int)total_reads_LOCAL
            );

        // TODO: Figure out alternative to callback into Python VM
        //       Cannot use in multi-threaded operation.
#if (0)
        // run callback, if specified
        if (callback && (0 == (total_reads_LOCAL % CALLBACK_PERIOD))) {
            try {
                callback(
                    "consume_fasta", callback_data,
                    total_reads_LOCAL, n_consumed_LOCAL
                );
            } catch (...) {
                throw;
            }
        }
#endif // 0

    } // while reads left for parser

} // consume_fasta

//
// consume_string: run through every k-mer in the given string, & hash it.
//

unsigned int Hashtable::consume_string(const std::string &s)
{
    const char * sp = s.c_str();
    unsigned int n_consumed = 0;

    KMerIterator kmers(sp, _ksize);

    while(!kmers.done()) {
        HashIntoType kmer = kmers.next();

        count(kmer);
        n_consumed++;
    }

    return n_consumed;
}

//
// consume_high_abund_kmers: run through every k-mer in the given string,
//     and if the k-mer is high enough abundance, add it. Return number
//     of consumed k-mers.
//

unsigned int Hashtable::consume_high_abund_kmers(const std::string &s,
        BoundedCounterType min_count)
{
    const char * sp = s.c_str();
    unsigned int n_consumed = 0;

    KMerIterator kmers(sp, _ksize);

    while(!kmers.done()) {
        HashIntoType kmer = kmers.next();

        if (get_count(kmer) >= min_count) {
            count(kmer);
            n_consumed++;
        }
    }

    return n_consumed;
}

// technically, get medioid count... our "median" is always a member of the
// population.

void Hashtable::get_median_count(const std::string &s,
                                 BoundedCounterType &median,
                                 float &average,
                                 float &stddev)
{
    std::vector<BoundedCounterType> counts;
    KMerIterator kmers(s.c_str(), _ksize);

    while(!kmers.done()) {
        HashIntoType kmer = kmers.next();
        BoundedCounterType count = this->get_count(kmer);
        counts.push_back(count);
    }

    if (!counts.size()) {
        throw std::exception();
    }

    if (!counts.size()) {
        median = 0;
        average = 0;
        stddev = 0;

        return;
    }

    average = 0;
    for (std::vector<BoundedCounterType>::const_iterator i = counts.begin();
            i != counts.end(); ++i) {
        average += *i;
    }
    average /= float(counts.size());

    stddev = 0;
    for (std::vector<BoundedCounterType>::const_iterator i = counts.begin();
            i != counts.end(); ++i) {
        stddev += (float(*i) - average) * (float(*i) - average);
    }
    stddev /= float(counts.size());
    stddev = sqrt(stddev);

    sort(counts.begin(), counts.end());
    median = counts[counts.size() / 2]; // rounds down
}

void Hashtable::save_tagset(std::string outfilename)
{
    ofstream outfile(outfilename.c_str(), ios::binary);
    const size_t tagset_size = n_tags();
    unsigned int save_ksize = _ksize;

    HashIntoType * buf = new HashIntoType[tagset_size];

    unsigned char version = SAVED_FORMAT_VERSION;
    outfile.write((const char *) &version, 1);

    unsigned char ht_type = SAVED_TAGS;
    outfile.write((const char *) &ht_type, 1);

    outfile.write((const char *) &save_ksize, sizeof(save_ksize));
    outfile.write((const char *) &tagset_size, sizeof(tagset_size));
    outfile.write((const char *) &_tag_density, sizeof(_tag_density));

    unsigned int i = 0;
    for (SeenSet::iterator pi = all_tags.begin(); pi != all_tags.end();
            ++pi, i++) {
        buf[i] = *pi;
    }

    outfile.write((const char *) buf, sizeof(HashIntoType) * tagset_size);
    outfile.close();

    delete[] buf;
}

void Hashtable::load_tagset(std::string infilename, bool clear_tags)
{
  ifstream infile;

  // configure ifstream to raise exceptions for everything.
  infile.exceptions(std::ifstream::failbit | std::ifstream::badbit |
                    std::ifstream::eofbit);

  try {
    infile.open(infilename.c_str(), ios::binary);
  } catch (std::ifstream::failure e) {
    std::string err;
    if (!(infile.is_open())) {
      err = "Cannot open file: " + infilename;
    }
    else {
      err = "Unknown error in opening file: " + infilename;
    }
    throw khmer_file_exception(err.c_str());
  }

    if (clear_tags) {
        all_tags.clear();
    }

    unsigned char version, ht_type;
    unsigned int save_ksize = 0;

    size_t tagset_size = 0;

    try {
      infile.read((char *) &version, 1);
      infile.read((char *) &ht_type, 1);
      if (!(version == SAVED_FORMAT_VERSION) || !(ht_type == SAVED_TAGS)) {
        std::string err = "File format error while reading tagset: " + \
          infilename;
        throw khmer_file_exception(err.c_str());
      }

      infile.read((char *) &save_ksize, sizeof(save_ksize));
      if (!(save_ksize == _ksize)) {
        std::string err = "Incorrect k-mer size while reading tagset: " + \
          infilename;
        throw khmer_file_exception(err.c_str());
      }

      infile.read((char *) &tagset_size, sizeof(tagset_size));
      infile.read((char *) &_tag_density, sizeof(_tag_density));

      HashIntoType * buf = new HashIntoType[tagset_size];

      infile.read((char *) buf, sizeof(HashIntoType) * tagset_size);

      for (unsigned int i = 0; i < tagset_size; i++) {
        all_tags.insert(buf[i]);
      }

      delete[] buf;
    } catch (std::ifstream::failure e) {
      std::string err = "Error reading data from: " + infilename;
      throw khmer_file_exception(err.c_str());
    }
}

void Hashtable::consume_sequence_and_tag(const std::string& seq,
        unsigned long long& n_consumed,
        SeenSet * found_tags)
{
    bool kmer_tagged;

    KMerIterator kmers(seq.c_str(), _ksize);
    HashIntoType kmer;

    unsigned int since = _tag_density / 2 + 1;

    while(!kmers.done()) {
        kmer = kmers.next();
        bool is_new_kmer;

        // Set the bits for the kmer in the various hashtables,
        // and report on whether or not they had already been set.
        // This is probably better than first testing and then setting the bits,
        // as a failed test essentially results in doing the same amount of work
        // twice.
        if ((is_new_kmer = test_and_set_bits( kmer ))) {
            ++n_consumed;
        }

#if (1)
        if (is_new_kmer) {
            ++since;
        } else {
            ACQUIRE_ALL_TAGS_SPIN_LOCK
            kmer_tagged = set_contains(all_tags, kmer);
            RELEASE_ALL_TAGS_SPIN_LOCK
            if (kmer_tagged) {
                since = 1;
                if (found_tags) {
                    found_tags->insert(kmer);
                }
            } else {
                ++since;
            }
        }
#else
        if (!is_new_kmer && set_contains(all_tags, kmer)) {
            since = 1;
            if (found_tags) {
                found_tags->insert(kmer);
            }
        } else {
            since++;
        }
#endif

        if (since >= _tag_density) {
            ACQUIRE_ALL_TAGS_SPIN_LOCK
            all_tags.insert(kmer);
            RELEASE_ALL_TAGS_SPIN_LOCK
            if (found_tags) {
                found_tags->insert(kmer);
            }
            since = 1;
        }

    } // iteration over kmers

    if (since >= _tag_density/2 - 1) {
        ACQUIRE_ALL_TAGS_SPIN_LOCK
        all_tags.insert(kmer);	// insert the last k-mer, too.
        RELEASE_ALL_TAGS_SPIN_LOCK
        if (found_tags) {
            found_tags->insert(kmer);
        }
    }
}

//
// consume_fasta_and_tag: consume a FASTA file of reads, tagging reads every
//     so often.
//

// TODO? Inline in header.
void
Hashtable::
consume_fasta_and_tag(
    std:: string const  &filename,
    unsigned int	      &total_reads, unsigned long long	&n_consumed,
    CallbackFn	      callback,	    void *		callback_data
)
{
    khmer:: Config    &the_config	  = khmer:: get_active_config( );

    // Note: Always assume only 1 thread if invoked this way.
    IParser *	  parser =
        IParser::get_parser(
            filename, 1, the_config.get_reads_input_buffer_size( ),
            the_config.get_reads_parser_trace_level( )
        );


    consume_fasta_and_tag(
        parser,
        total_reads, n_consumed,
        callback, callback_data
    );

    delete parser;
}

void
Hashtable::
consume_fasta_and_tag(
    read_parsers:: IParser *  parser,
    unsigned int		    &total_reads,   unsigned long long	&n_consumed,
    CallbackFn		    callback,	    void *		callback_data
)
{
    Hasher		  &hasher		=
        _get_hasher( parser->uuid( ) );
    unsigned int		  total_reads_LOCAL	= 0;
#if (0) // Note: Used with callback - currently disabled.
    unsigned long long int  n_consumed_LOCAL	= 0;
#endif
    Read			  read;

    // TODO? Delete the following assignments.
    total_reads = 0;
    n_consumed = 0;

    hasher.trace_logger(
        TraceLogger:: TLVL_DEBUG2,
        "Starting trace of 'consume_fasta_and_tag'....\n"
    );

    // Iterate through the reads and consume their k-mers.
    while (!parser->is_complete( )) {

        read = parser->get_next_read( );

        if (check_and_normalize_read( read.sequence )) {
            unsigned long long this_n_consumed = 0;
            consume_sequence_and_tag( read.sequence, this_n_consumed );

#ifdef WITH_INTERNAL_METRICS
            hasher.pmetrics.start_timers( );
#endif
#if (0) // Note: Used with callback - currently disabled.
            n_consumed_LOCAL  = __sync_add_and_fetch( &n_consumed, this_n_consumed );
#else
            __sync_add_and_fetch( &n_consumed, this_n_consumed );
#endif
            total_reads_LOCAL = __sync_add_and_fetch( &total_reads, 1 );
#ifdef WITH_INTERNAL_METRICS
            hasher.pmetrics.stop_timers( );
            hasher.pmetrics.accumulate_timer_deltas(
                (uint32_t)HashTablePerformanceMetrics:: MKEY_TIME_UPDATE_TALLIES
            );
#endif
        }

        if (0 == (total_reads_LOCAL % 10000))
            hasher.trace_logger(
                TraceLogger:: TLVL_DEBUG3,
                "Total number of reads processed: %llu\n",
                (unsigned long long int)total_reads_LOCAL
            );

        // TODO: Figure out alternative to callback into Python VM
        //       Cannot use in multi-threaded operation.
#if (0)
        // run callback, if specified
        if (total_reads_TL % CALLBACK_PERIOD == 0 && callback) {
            std::cout << "n tags: " << all_tags.size() << "\n";
            try {
                callback("consume_fasta_and_tag", callback_data, total_reads_TL,
                         n_consumed);
            } catch (...) {
                delete parser;
                throw;
            }
        }
#endif // 0

    } // while reads left for parser

}

//
// consume_fasta_and_tag_with_stoptags: consume a FASTA file of reads,
//     tagging reads every so often.  Do not insert matches to stoptags,
//     and join the tags across those gaps.
//

void Hashtable::consume_fasta_and_tag_with_stoptags(const std::string &filename,
        unsigned int &total_reads,
        unsigned long long &n_consumed,
        CallbackFn callback,
        void * callback_data)
{
    total_reads = 0;
    n_consumed = 0;

    IParser* parser = IParser::get_parser(filename.c_str());
    Read read;

    string seq = "";

    SeenSet read_tags;

    //
    // iterate through the FASTA file & consume the reads.
    //

    while(!parser->is_complete())  {
        read = parser->get_next_read();
        seq = read.sequence;

        read_tags.clear();

        // n_consumed += this_n_consumed;

        if (check_and_normalize_read(seq)) {	// process?
            bool is_new_kmer;
            KMerIterator kmers(seq.c_str(), _ksize);

            HashIntoType kmer, last_kmer;
            bool is_first_kmer = true;

            unsigned int since = _tag_density / 2 + 1;
            while (!kmers.done()) {
                kmer = kmers.next();

                if (!set_contains(stop_tags, kmer)) { // NOT a stop tag... ok.
                    is_new_kmer = (bool) !get_count(kmer);
                    if (is_new_kmer) {
                        count(kmer);
                        n_consumed++;
                    }

                    if (!is_new_kmer && set_contains(all_tags, kmer)) {
                        read_tags.insert(kmer);
                        since = 1;
                    } else {
                        since++;
                    }

                    if (since >= _tag_density) {
                        all_tags.insert(kmer);
                        read_tags.insert(kmer);
                        since = 1;
                    }
                } else {		// stop tag!  do not insert, but connect.
                    // before first tag insertion; insert last kmer.
                    if (!is_first_kmer && read_tags.empty()) {
                        read_tags.insert(last_kmer);
                        all_tags.insert(last_kmer);
                    }

                    since = _tag_density - 1; // insert next kmer, too.
                }

                last_kmer = kmer;
                is_first_kmer = false;
            }

            if (!set_contains(stop_tags, kmer)) { // NOT a stop tag... ok.
                is_new_kmer = (bool) !get_count(kmer);
                if (is_new_kmer) {
                    count(kmer);
                    n_consumed++;
                }

                if (since >= _tag_density/2 - 1) {
                    all_tags.insert(kmer);	// insert the last k-mer, too.
                    read_tags.insert(kmer);
                }
            }
        }

        if (read_tags.size() > 1) {
            partition->assign_partition_id(*(read_tags.begin()), read_tags);
        }

        // reset the sequence info, increment read number
        total_reads++;

        // run callback, if specified
        if (total_reads % CALLBACK_PERIOD == 0 && callback) {
            std::cout << "n tags: " << all_tags.size() << "\n";
            try {
                callback("consume_fasta_and_tag", callback_data, total_reads,
                         n_consumed);
            } catch (...) {
                delete parser;
                throw;
            }
        }
    }
    delete parser;
}

//
// divide_tags_into_subsets - take all of the tags in 'all_tags', and
//   divide them into subsets (based on starting tag) of <= given size.
//

void Hashtable::divide_tags_into_subsets(unsigned int subset_size,
        SeenSet& divvy)
{
    unsigned int i = 0;

    for (SeenSet::const_iterator si = all_tags.begin(); si != all_tags.end();
            ++si) {
        if (i % subset_size == 0) {
            divvy.insert(*si);
            i = 0;
        }
        i++;
    }
}

//
// consume_partitioned_fasta: consume a FASTA file of reads
//

void Hashtable::consume_partitioned_fasta(const std::string &filename,
        unsigned int &total_reads,
        unsigned long long &n_consumed,
        CallbackFn callback,
        void * callback_data)
{
    total_reads = 0;
    n_consumed = 0;

    IParser* parser = IParser::get_parser(filename.c_str());
    Read read;

    string seq = "";

    // reset the master subset partition
    delete partition;
    partition = new SubsetPartition(this);

    //
    // iterate through the FASTA file & consume the reads.
    //

    while(!parser->is_complete())  {
        read = parser->get_next_read();
        seq = read.sequence;

        if (check_and_normalize_read(seq)) {
            // First, figure out what the partition is (if non-zero), and save that.
            PartitionID p = _parse_partition_id(read.name);

            // Then consume the sequence
            n_consumed += consume_string(seq); // @CTB why are we doing this?

            // Next, compute the tag & set the partition, if nonzero
            HashIntoType kmer = _hash(seq.c_str(), _ksize);
            all_tags.insert(kmer);
            if (p > 0) {
                partition->set_partition_id(kmer, p);
            }
        }

        // reset the sequence info, increment read number
        total_reads++;

        // run callback, if specified
        if (total_reads % CALLBACK_PERIOD == 0 && callback) {
            try {
                callback("consume_partitioned_fasta", callback_data, total_reads,
                         n_consumed);
            } catch (...) {
                delete parser;
                throw;
            }
        }
    }

    delete parser;
}

//
// consume_fasta: consume a FASTA file of reads
//

void Hashtable::consume_fasta_and_traverse(const std::string &filename,
        unsigned int radius,
        unsigned int big_threshold,
        unsigned int transfer_threshold,
        CountingHash &counting)
{
    unsigned long long total_reads = 0;

    IParser* parser = IParser::get_parser(filename.c_str());
    Read read;

    string seq = "";

    //
    // iterate through the FASTA file & consume the reads.
    //

    while(!parser->is_complete())  {
        read = parser->get_next_read();
        seq = read.sequence;

        if (check_and_normalize_read(seq)) {	// process?
            KMerIterator kmers(seq.c_str(), _ksize);

            HashIntoType kmer = 0;
            bool is_first_kmer = true;
            while (!kmers.done()) {
                kmer = kmers.next();

                if (set_contains(stop_tags, kmer)) {
                    break;
                }
                count(kmer);
                is_first_kmer = false;
            }

            if (!is_first_kmer) {	// traverse
                SeenSet keeper;

                unsigned int n = traverse_from_kmer(kmer, radius, keeper);
                if (n >= big_threshold) {
#if VERBOSE_REPARTITION
                    std::cout << "lmp: " << n << "; added: " << stop_tags.size() << "\n";
#endif // VERBOSE_REPARTITION
                    count_and_transfer_to_stoptags(keeper, transfer_threshold, counting);
                }
            }
        }

        // reset the sequence info, increment read number
        total_reads++;

        // run callback, if specified
        if (total_reads % CALLBACK_PERIOD == 0) {
            std::cout << "n reads: " << total_reads << "\n";
        }
    }
    delete parser;
}

//////////////////////////////////////////////////////////////////////
// graph stuff

void Hashtable::calc_connected_graph_size(const HashIntoType kmer_f,
        const HashIntoType kmer_r,
        unsigned long long& count,
        SeenSet& keeper,
        const unsigned long long threshold,
        bool break_on_circum)
const
{
    HashIntoType kmer = uniqify_rc(kmer_f, kmer_r);
    const BoundedCounterType val = get_count(kmer);

    if (val == 0) {
        return;
    }

    // have we already seen me? don't count; exit.
    if (set_contains(keeper, kmer)) {
        return;
    }

    // is this in stop_tags?
    if (set_contains(stop_tags, kmer)) {
        return;
    }

    // keep track of both seen kmers, and counts.
    keeper.insert(kmer);

    // is this a high-circumference k-mer? if so, don't count it; get outta here!
    if (break_on_circum && \
            kmer_degree(kmer_f, kmer_r) > 4) {
        return;
    }

    count += 1;

    // are we past the threshold? truncate search.
    if (threshold && count >= threshold) {
        return;
    }

    // otherwise, explore in all directions.

    // NEXT.

    HashIntoType f, r;
    const unsigned int rc_left_shift = _ksize*2 - 2;

    f = next_f(kmer_f, 'A');
    r = next_r(kmer_r, 'A');
    calc_connected_graph_size(f, r, count, keeper, threshold, break_on_circum);

    f = next_f(kmer_f, 'C');
    r = next_r(kmer_r, 'C');
    calc_connected_graph_size(f, r, count, keeper, threshold, break_on_circum);

    f = next_f(kmer_f, 'G');
    r = next_r(kmer_r, 'G');
    calc_connected_graph_size(f, r, count, keeper, threshold, break_on_circum);

    f = next_f(kmer_f, 'T');
    r = next_r(kmer_r, 'T');
    calc_connected_graph_size(f, r, count, keeper, threshold, break_on_circum);

    // PREVIOUS.


    r = prev_r(kmer_r, 'A');
    f = prev_f(kmer_f, 'A');
    calc_connected_graph_size(f, r, count, keeper, threshold, break_on_circum);

    r = prev_r(kmer_r, 'C');
    f = prev_f(kmer_f, 'C');
    calc_connected_graph_size(f, r, count, keeper, threshold, break_on_circum);

    r = prev_r(kmer_r, 'G');
    f = prev_f(kmer_f, 'G');
    calc_connected_graph_size(f, r, count, keeper, threshold, break_on_circum);

    r = prev_r(kmer_r, 'T');
    f = prev_f(kmer_f, 'T');
    calc_connected_graph_size(f, r, count, keeper, threshold, break_on_circum);
}

unsigned int Hashtable::kmer_degree(HashIntoType kmer_f, HashIntoType kmer_r)
const
{
    unsigned int neighbors = 0;

    const unsigned int rc_left_shift = _ksize*2 - 2;

    HashIntoType f, r;

    // NEXT.
    f = next_f(kmer_f, 'A');
    r = next_r(kmer_r, 'A');
    if (get_count(uniqify_rc(f, r))) {
        neighbors++;
    }

    f = next_f(kmer_f, 'C');
    r = next_r(kmer_r, 'C');
    if (get_count(uniqify_rc(f, r))) {
        neighbors++;
    }

    f = next_f(kmer_f, 'G');
    r = next_r(kmer_r, 'G');
    if (get_count(uniqify_rc(f, r))) {
        neighbors++;
    }

    f = next_f(kmer_f, 'T');
    r = next_r(kmer_r, 'T');
    if (get_count(uniqify_rc(f, r))) {
        neighbors++;
    }

    // PREVIOUS.
    r = prev_r(kmer_r, 'A');
    f = prev_f(kmer_f, 'A');
    if (get_count(uniqify_rc(f, r))) {
        neighbors++;
    }

    r = prev_r(kmer_r, 'C');
    f = prev_f(kmer_f, 'C');
    if (get_count(uniqify_rc(f, r))) {
        neighbors++;
    }

    r = prev_r(kmer_r, 'G');
    f = prev_f(kmer_f, 'G');
    if (get_count(uniqify_rc(f, r))) {
        neighbors++;
    }

    r = prev_r(kmer_r, 'T');
    f = prev_f(kmer_f, 'T');
    if (get_count(uniqify_rc(f, r))) {
        neighbors++;
    }

    return neighbors;
}

void Hashtable::filter_if_present(const std::string infilename,
                                  const std::string outputfile,
                                  CallbackFn callback,
                                  void * callback_data)
{
    IParser* parser = IParser::get_parser(infilename);
    ofstream outfile(outputfile.c_str());

    unsigned int total_reads = 0;
    unsigned int reads_kept = 0;

    Read read;
    string seq;

    HashIntoType kmer;

    while(!parser->is_complete()) {
        read = parser->get_next_read();
        seq = read.sequence;

        if (check_and_normalize_read(seq)) {
            KMerIterator kmers(seq.c_str(), _ksize);
            bool keep = true;

            while (!kmers.done()) {
                kmer = kmers.next();
                if (get_count(kmer)) {
                    keep = false;
                    break;
                }
            }

            if (keep) {
                outfile << ">" << read.name;
                outfile << "\n" << seq << "\n";
                reads_kept++;
            }

            total_reads++;

            // run callback, if specified
            if (total_reads % CALLBACK_PERIOD == 0 && callback) {
                try {
                    callback("filter_if_present", callback_data,total_reads, reads_kept);
                } catch (...) {
                    delete parser;
                    parser = NULL;
                    outfile.close();
                    throw;
                }
            }
        }
    }

    delete parser;
    parser = NULL;

    return;
}


unsigned int Hashtable::count_kmers_within_radius(HashIntoType kmer_f,
        HashIntoType kmer_r,
        unsigned int radius,
        unsigned int max_count,
        const SeenSet * seen)
const
{
    HashIntoType f, r;
    NodeQueue node_q;
    std::queue<unsigned int> breadth_q;
    unsigned int cur_breadth = 0;

    const unsigned int rc_left_shift = _ksize*2 - 2;
    unsigned int total = 0;

    SeenSet keeper;		// keep track of traversed kmers
    if (seen) {
        keeper = *seen;
    }

    // start breadth-first search.

    node_q.push(kmer_f);
    node_q.push(kmer_r);
    breadth_q.push(0);

    while(!node_q.empty()) {
        kmer_f = node_q.front();
        node_q.pop();
        kmer_r = node_q.front();
        node_q.pop();
        unsigned int breadth = breadth_q.front();
        breadth_q.pop();

        if (breadth > radius) {
            break;
        }

        HashIntoType kmer = uniqify_rc(kmer_f, kmer_r);
        if (set_contains(keeper, kmer)) {
            continue;
        }

        // keep track of seen kmers
        keeper.insert(kmer);
        total++;

        if (max_count && total > max_count) {
            break;
        }

        if (!(breadth >= cur_breadth)) { // keep track of watermark, for debugging.
            throw std::exception();
        }
        if (breadth > cur_breadth) {
            cur_breadth = breadth;
        }

        //
        // Enqueue next set of nodes.
        //

        // NEXT.
        f = next_f(kmer_f, 'A');
        r = next_r(kmer_r, 'A');
        if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f, r))) {
            node_q.push(f);
            node_q.push(r);
            breadth_q.push(breadth + 1);
        }

        f = next_f(kmer_f, 'C');
        r = next_r(kmer_r, 'C');
        if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f, r))) {
            node_q.push(f);
            node_q.push(r);
            breadth_q.push(breadth + 1);
        }

        f = next_f(kmer_f, 'G');
        r = next_r(kmer_r, 'G');
        if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f, r))) {
            node_q.push(f);
            node_q.push(r);
            breadth_q.push(breadth + 1);
        }

        f = next_f(kmer_f, 'T');
        r = next_r(kmer_r, 'T');
        if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f, r))) {
            node_q.push(f);
            node_q.push(r);
            breadth_q.push(breadth + 1);
        }

        // PREVIOUS.
        r = prev_r(kmer_r, 'A');
        f = prev_f(kmer_f, 'A');
        if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f, r))) {
            node_q.push(f);
            node_q.push(r);
            breadth_q.push(breadth + 1);
        }

        r = prev_r(kmer_r, 'C');
        f = prev_f(kmer_f, 'C');
        if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f, r))) {
            node_q.push(f);
            node_q.push(r);
            breadth_q.push(breadth + 1);
        }

        r = prev_r(kmer_r, 'G');
        f = prev_f(kmer_f, 'G');
        if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f, r))) {
            node_q.push(f);
            node_q.push(r);
            breadth_q.push(breadth + 1);
        }

        r = prev_r(kmer_r, 'T');
        f = prev_f(kmer_f, 'T');
        if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f, r))) {
            node_q.push(f);
            node_q.push(r);
            breadth_q.push(breadth + 1);
        }
    }

    return total;
}

unsigned int Hashtable::count_kmers_within_depth(HashIntoType kmer_f,
        HashIntoType kmer_r,
        unsigned int depth,
        unsigned int max_count,
        SeenSet * seen)
const
{
    HashIntoType f, r;
    unsigned int count = 1;

    if (depth == 0) {
        return 0;
    }

    const unsigned int rc_left_shift = _ksize*2 - 2;

    seen->insert(uniqify_rc(kmer_f, kmer_r));

    // NEXT.
    f = next_f(kmer_f, 'A');
    r = next_r(kmer_r, 'A');
    if (get_count(uniqify_rc(f,r)) && !set_contains(*seen, uniqify_rc(f, r))) {
        count += count_kmers_within_depth(f, r, depth - 1, max_count - count,
                                          seen);
        if (count >= max_count) {
            return count;
        }
    }

    f = next_f(kmer_f, 'C');
    r = next_r(kmer_r, 'C');
    if (get_count(uniqify_rc(f,r)) && !set_contains(*seen, uniqify_rc(f, r))) {
        count += count_kmers_within_depth(f, r, depth -1, max_count - count,
                                          seen);
        if (count >= max_count) {
            return count;
        }
        ;
    }

    f = next_f(kmer_f, 'G');
    r = next_r(kmer_r, 'G');
    if (get_count(uniqify_rc(f,r)) && !set_contains(*seen, uniqify_rc(f, r))) {
        count += count_kmers_within_depth(f, r, depth -1, max_count - count,
                                          seen);
        if (count >= max_count) {
            return count;
        }
        ;
    }

    f = next_f(kmer_f, 'T');
    r = next_r(kmer_r, 'T');
    if (get_count(uniqify_rc(f,r)) && !set_contains(*seen, uniqify_rc(f, r))) {
        count += count_kmers_within_depth(f, r, depth -1, max_count - count,
                                          seen);
        if (count >= max_count) {
            return count;
        }
        ;
    }

    // PREVIOUS.
    r = prev_r(kmer_r, 'A');
    f = prev_f(kmer_f, 'A');
    if (get_count(uniqify_rc(f,r)) && !set_contains(*seen, uniqify_rc(f, r))) {
        count += count_kmers_within_depth(f, r, depth -1, max_count - count,
                                          seen);
        if (count >= max_count) {
            return count;
        }
        ;
    }

    r = prev_r(kmer_r, 'C');
    f = prev_f(kmer_f, 'C');
    if (get_count(uniqify_rc(f,r)) && !set_contains(*seen, uniqify_rc(f, r))) {
        count += count_kmers_within_depth(f, r, depth -1, max_count - count,
                                          seen);
        if (count >= max_count) {
            return count;
        }
        ;
    }

    r = prev_r(kmer_r, 'G');
    f = prev_f(kmer_f, 'G');
    if (get_count(uniqify_rc(f,r)) && !set_contains(*seen, uniqify_rc(f, r))) {
        count += count_kmers_within_depth(f, r, depth -1, max_count - count,
                                          seen);
        if (count >= max_count) {
            return count;
        }
        ;
    }

    r = prev_r(kmer_r, 'T');
    f = prev_f(kmer_f, 'T');
    if (get_count(uniqify_rc(f,r)) && !set_contains(*seen, uniqify_rc(f, r))) {
        count += count_kmers_within_depth(f, r, depth -1, max_count - count,
                                          seen);
        if (count >= max_count) {
            return count;
        }
        ;
    }

    return count;
}

unsigned int Hashtable::find_radius_for_volume(HashIntoType kmer_f,
        HashIntoType kmer_r,
        unsigned int max_count,
        unsigned int max_radius)
const
{
    HashIntoType f, r;
    NodeQueue node_q;
    std::queue<unsigned int> breadth_q;
    unsigned int breadth = 0;

    const unsigned int rc_left_shift = _ksize*2 - 2;
    unsigned int total = 0;

    SeenSet keeper;		// keep track of traversed kmers

    // start breadth-first search.

    node_q.push(kmer_f);
    node_q.push(kmer_r);
    breadth_q.push(0);

    while(!node_q.empty()) {
        kmer_f = node_q.front();
        node_q.pop();
        kmer_r = node_q.front();
        node_q.pop();
        breadth = breadth_q.front();
        breadth_q.pop();

        HashIntoType kmer = uniqify_rc(kmer_f, kmer_r);
        if (set_contains(keeper, kmer)) {
            continue;
        }

        // keep track of seen kmers
        keeper.insert(kmer);
        total++;

        if (total >= max_count || breadth >= max_radius) {
            break;
        }

        //
        // Enqueue next set of nodes.
        //

        // NEXT.
        f = next_f(kmer_f, 'A');
        r = next_r(kmer_r, 'A');
        if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f,r))) {
            node_q.push(f);
            node_q.push(r);
            breadth_q.push(breadth + 1);
        }

        f = next_f(kmer_f, 'C');
        r = next_r(kmer_r, 'C');
        if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f,r))) {
            node_q.push(f);
            node_q.push(r);
            breadth_q.push(breadth + 1);
        }

        f = next_f(kmer_f, 'G');
        r = next_r(kmer_r, 'G');
        if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f,r))) {
            node_q.push(f);
            node_q.push(r);
            breadth_q.push(breadth + 1);
        }

        f = next_f(kmer_f, 'T');
        r = next_r(kmer_r, 'T');
        if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f,r))) {
            node_q.push(f);
            node_q.push(r);
            breadth_q.push(breadth + 1);
        }

        // PREVIOUS.
        r = prev_r(kmer_r, 'A');
        f = prev_f(kmer_f, 'A');
        if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f,r))) {
            node_q.push(f);
            node_q.push(r);
            breadth_q.push(breadth + 1);
        }

        r = prev_r(kmer_r, 'C');
        f = prev_f(kmer_f, 'C');
        if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f,r))) {
            node_q.push(f);
            node_q.push(r);
            breadth_q.push(breadth + 1);
        }

        r = prev_r(kmer_r, 'G');
        f = prev_f(kmer_f, 'G');
        if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f,r))) {
            node_q.push(f);
            node_q.push(r);
            breadth_q.push(breadth + 1);
        }

        r = prev_r(kmer_r, 'T');
        f = prev_f(kmer_f, 'T');
        if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f,r))) {
            node_q.push(f);
            node_q.push(r);
            breadth_q.push(breadth + 1);
        }

        if (node_q.empty()) {
            breadth = max_radius;
            break;
        }
    }

    return breadth;
}

unsigned int Hashtable::count_kmers_on_radius(HashIntoType kmer_f,
        HashIntoType kmer_r,
        unsigned int radius,
        unsigned int max_volume)
const
{
    HashIntoType f, r;
    NodeQueue node_q;
    std::queue<unsigned int> breadth_q;
    unsigned int cur_breadth = 0;
    unsigned int count = 0;

    const unsigned int rc_left_shift = _ksize*2 - 2;
    unsigned int total = 0;

    SeenSet keeper;		// keep track of traversed kmers

    // start breadth-first search.

    node_q.push(kmer_f);
    node_q.push(kmer_r);
    breadth_q.push(0);

    while(!node_q.empty()) {
        kmer_f = node_q.front();
        node_q.pop();
        kmer_r = node_q.front();
        node_q.pop();
        unsigned int breadth = breadth_q.front();
        breadth_q.pop();

        if (breadth > radius) {
            break;
        }

        HashIntoType kmer = uniqify_rc(kmer_f, kmer_r);
        if (set_contains(keeper, kmer)) {
            continue;
        }

        if (breadth == radius) {
            count++;
        }

        // keep track of seen kmers
        keeper.insert(kmer);
        total++;

        if (max_volume && total > max_volume) {
            break;
        }

        if (!(breadth >= cur_breadth)) { // keep track of watermark, for debugging.
            throw std::exception();
        }
        if (breadth > cur_breadth) {
            cur_breadth = breadth;
        }

        //
        // Enqueue next set of nodes.
        //

        // NEXT.
        f = next_f(kmer_f, 'A');
        r = next_r(kmer_r, 'A');
        if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f,r))) {
            node_q.push(f);
            node_q.push(r);
            breadth_q.push(breadth + 1);
        }

        f = next_f(kmer_f, 'C');
        r = next_r(kmer_r, 'C');
        if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f,r))) {
            node_q.push(f);
            node_q.push(r);
            breadth_q.push(breadth + 1);
        }

        f = next_f(kmer_f, 'G');
        r = next_r(kmer_r, 'G');
        if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f,r))) {
            node_q.push(f);
            node_q.push(r);
            breadth_q.push(breadth + 1);
        }

        f = next_f(kmer_f, 'T');
        r = next_r(kmer_r, 'T');
        if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f,r))) {
            node_q.push(f);
            node_q.push(r);
            breadth_q.push(breadth + 1);
        }

        // PREVIOUS.
        r = prev_r(kmer_r, 'A');
        f = prev_f(kmer_f, 'A');
        if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f,r))) {
            node_q.push(f);
            node_q.push(r);
            breadth_q.push(breadth + 1);
        }

        r = prev_r(kmer_r, 'C');
        f = prev_f(kmer_f, 'C');
        if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f,r))) {
            node_q.push(f);
            node_q.push(r);
            breadth_q.push(breadth + 1);
        }

        r = prev_r(kmer_r, 'G');
        f = prev_f(kmer_f, 'G');
        if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f,r))) {
            node_q.push(f);
            node_q.push(r);
            breadth_q.push(breadth + 1);
        }

        r = prev_r(kmer_r, 'T');
        f = prev_f(kmer_f, 'T');
        if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f,r))) {
            node_q.push(f);
            node_q.push(r);
            breadth_q.push(breadth + 1);
        }
    }

    return count;
}

size_t Hashtable::trim_on_stoptags(std::string seq) const
{
    if (!check_and_normalize_read(seq)) {
        return 0;
    }

    KMerIterator kmers(seq.c_str(), _ksize);

    size_t i = _ksize - 2;
    while (!kmers.done()) {
        HashIntoType kmer = kmers.next();
        if (set_contains(stop_tags, kmer)) {
            return i;
        }
        i++;
    }

    return seq.length();
}

void Hashtable::traverse_from_tags(unsigned int distance,
                                   unsigned int threshold,
                                   unsigned int frequency,
                                   CountingHash &counting)
{
    unsigned int i = 0;
    unsigned int n = 0;
    unsigned int n_big = 0;
    SeenSet keeper;

#if VERBOSE_REPARTITION
    std::cout << all_tags.size() << " tags...\n";
#endif // 0

    for (SeenSet::const_iterator si = all_tags.begin(); si != all_tags.end();
            ++si, i++) {
        n++;
        unsigned int count = traverse_from_kmer(*si, distance, keeper);

        if (count >= threshold) {
            n_big++;

            SeenSet::const_iterator ti;
            for (ti = keeper.begin(); ti != keeper.end(); ++ti) {
                if (counting.get_count(*ti) > frequency) {
                    stop_tags.insert(*ti);
                } else {
                    counting.count(*ti);
                }
            }
#if VERBOSE_REPARTITION
            std::cout << "traversed from " << n << " tags total; "
                      << n_big << " big; " << keeper.size() << "\n";
#endif // 0
        }
        keeper.clear();

        if (n % 100 == 0) {
#if VERBOSE_REPARTITION
            std::cout << "traversed " << n << " " << n_big << " " <<
                      all_tags.size() << " " << stop_tags.size() << "\n";
#endif // 0
        }
    }
}

unsigned int Hashtable::traverse_from_kmer(HashIntoType start,
        unsigned int radius,
        SeenSet &keeper)
const
{
    std::string kmer_s = _revhash(start, _ksize);
    HashIntoType kmer_f, kmer_r;
    _hash(kmer_s.c_str(), _ksize, kmer_f, kmer_r);

    HashIntoType f, r;
    NodeQueue node_q;
    std::queue<unsigned int> breadth_q;
    unsigned int cur_breadth = 0;
    bool is_first_kmer = true;

    const unsigned int rc_left_shift = _ksize*2 - 2;
    unsigned int total = 0;

    // start breadth-first search.

    node_q.push(kmer_f);
    node_q.push(kmer_r);
    breadth_q.push(0);

    while(!node_q.empty()) {
        kmer_f = node_q.front();
        node_q.pop();
        kmer_r = node_q.front();
        node_q.pop();
        unsigned int breadth = breadth_q.front();
        breadth_q.pop();

        if (breadth > radius) {
            break;
        }

        if (total > MAX_KEEPER_SIZE) {
            break;
        }

        HashIntoType kmer = uniqify_rc(kmer_f, kmer_r);
        if (set_contains(keeper, kmer)) {
            continue;
        }

        if (set_contains(stop_tags, kmer)) {
            continue;
        }

        // keep track of seen kmers
        keeper.insert(kmer);
        total++;

        // QUESTION: Huh? What's up with the following?
        if (false && !is_first_kmer && set_contains(all_tags, kmer)) {
            continue;
        }

        if (!(breadth >= cur_breadth)) { // keep track of watermark, for debugging.
            throw std::exception();
        }
        if (breadth > cur_breadth) {
            cur_breadth = breadth;
        }

        //
        // Enqueue next set of nodes.
        //

        // NEXT.
        f = next_f(kmer_f, 'A');
        //r = next_r(kmer_r, 'A');

        // f = ((kmer_f << 2) & bitmask) | twobit_repr('A');
        r = next_r(kmer_r, 'A');
        if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f,r))) {
            node_q.push(f);
            node_q.push(r);
            breadth_q.push(breadth + 1);
        }

        f = next_f(kmer_f, 'C');
        r = next_r(kmer_r, 'C');
        if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f,r))) {
            node_q.push(f);
            node_q.push(r);
            breadth_q.push(breadth + 1);
        }

        f = next_f(kmer_f, 'G');
        r = next_r(kmer_r, 'G');
        if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f,r))) {
            node_q.push(f);
            node_q.push(r);
            breadth_q.push(breadth + 1);
        }

        f = next_f(kmer_f, 'T');
        r = next_r(kmer_r, 'T');
        if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f,r))) {
            node_q.push(f);
            node_q.push(r);
            breadth_q.push(breadth + 1);
        }

        // PREVIOUS.
        r = prev_r(kmer_r, 'A');
        f = prev_f(kmer_f, 'A');
        if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f,r))) {
            node_q.push(f);
            node_q.push(r);
            breadth_q.push(breadth + 1);
        }

        r = prev_r(kmer_r, 'C');
        f = prev_f(kmer_f, 'C');
        if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f,r))) {
            node_q.push(f);
            node_q.push(r);
            breadth_q.push(breadth + 1);
        }

        r = prev_r(kmer_r, 'G');
        f = prev_f(kmer_f, 'G');
        if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f,r))) {
            node_q.push(f);
            node_q.push(r);
            breadth_q.push(breadth + 1);
        }

        r = prev_r(kmer_r, 'T');
        f = prev_f(kmer_f, 'T');
        if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f,r))) {
            node_q.push(f);
            node_q.push(r);
            breadth_q.push(breadth + 1);
        }

        is_first_kmer = false;
    }

    return total;
}

void Hashtable::load_stop_tags(std::string infilename, bool clear_tags)
{
  ifstream infile;

  // configure ifstream to raise exceptions for everything.
  infile.exceptions(std::ifstream::failbit | std::ifstream::badbit |
                    std::ifstream::eofbit);

  try {
    infile.open(infilename.c_str(), ios::binary);
  } catch (std::ifstream::failure e) {
    std::string err;
    if (!(infile.is_open())) {
      err = "Cannot open file: " + infilename;
    }
    else {
      err = "Unknown error in opening file: " + infilename;
    }
    throw khmer_file_exception(err.c_str());
  }

    if (clear_tags) {
        stop_tags.clear();
    }

    unsigned char version, ht_type;
    unsigned int save_ksize = 0;

    size_t tagset_size = 0;

    try {
      infile.read((char *) &version, 1);
      infile.read((char *) &ht_type, 1);
      if (!(version == SAVED_FORMAT_VERSION) || !(ht_type == SAVED_STOPTAGS)) {
        std::string err = "File format error while reading stoptags: " + \
          infilename;
        throw khmer_file_exception(err.c_str());
      }

      infile.read((char *) &save_ksize, sizeof(save_ksize));
      if (!(save_ksize == _ksize)) {
        std::string err = "Incorrect k-mer size while reading stoptags: " + \
          infilename;
        throw khmer_file_exception(err.c_str());
      }
      infile.read((char *) &tagset_size, sizeof(tagset_size));

      HashIntoType * buf = new HashIntoType[tagset_size];

      infile.read((char *) buf, sizeof(HashIntoType) * tagset_size);

      for (unsigned int i = 0; i < tagset_size; i++) {
        stop_tags.insert(buf[i]);
      }
      delete[] buf;
    }
    catch (std::ifstream::failure e) {
      std::string err = "Error reading data from: " + infilename;
      throw khmer_file_exception(err.c_str());
    }
}

void Hashtable::save_stop_tags(std::string outfilename)
{
    ofstream outfile(outfilename.c_str(), ios::binary);
    size_t tagset_size = stop_tags.size();

    HashIntoType * buf = new HashIntoType[tagset_size];

    unsigned char version = SAVED_FORMAT_VERSION;
    outfile.write((const char *) &version, 1);

    unsigned char ht_type = SAVED_STOPTAGS;
    outfile.write((const char *) &ht_type, 1);

    unsigned int save_ksize = _ksize;
    outfile.write((const char *) &save_ksize, sizeof(save_ksize));
    outfile.write((const char *) &tagset_size, sizeof(tagset_size));

    unsigned int i = 0;
    for (SeenSet::iterator pi = stop_tags.begin(); pi != stop_tags.end();
            ++pi, i++) {
        buf[i] = *pi;
    }

    outfile.write((const char *) buf, sizeof(HashIntoType) * tagset_size);
    outfile.close();

    delete[] buf;
}

void Hashtable::print_stop_tags(std::string infilename)
{
    ofstream printfile(infilename.c_str());

    unsigned int i = 0;
    for (SeenSet::iterator pi = stop_tags.begin(); pi != stop_tags.end();
            ++pi, i++) {
        std::string kmer = _revhash(*pi, _ksize);
        printfile << kmer << "\n";
    }

    printfile.close();
}

void Hashtable::print_tagset(std::string infilename)
{
    ofstream printfile(infilename.c_str());

    unsigned int i = 0;
    for (SeenSet::iterator pi = all_tags.begin(); pi != all_tags.end();
            ++pi, i++) {
        std::string kmer = _revhash(*pi, _ksize);
        printfile << kmer << "\n";
    }

    printfile.close();
}

unsigned int Hashtable::count_and_transfer_to_stoptags(SeenSet &keeper,
        unsigned int threshold,
        CountingHash &counting)
{
    unsigned int n_inserted = 0;

    SeenSet::const_iterator ti;
    for (ti = keeper.begin(); ti != keeper.end(); ++ti) {
        if (counting.get_count(*ti) >= threshold) {
            stop_tags.insert(*ti);
            n_inserted++;
        } else {
            counting.count(*ti);
        }
    }

    return n_inserted;
}

void Hashtable::identify_stop_tags_by_position(std::string seq,
        std::vector<unsigned int> &posns)
const
{
    if (!check_and_normalize_read(seq)) {
        return;
    }

    KMerIterator kmers(seq.c_str(), _ksize);

    unsigned int i = 0;
    while(!kmers.done()) {
        HashIntoType kmer = kmers.next();

        if (set_contains(stop_tags, kmer)) {
            posns.push_back(i);
        }
        i++;
    }

    return;
}

void Hashtable::extract_unique_paths(std::string seq,
                                     unsigned int min_length,
                                     float min_unique_f,
                                     std::vector<std::string> &results)
{
    if (seq.size() < min_length) {
        return;
    }

    float max_seen = 1.0 - min_unique_f;

    min_length = min_length - _ksize + 1; // adjust for k-mer size.

    KMerIterator kmers(seq.c_str(), _ksize);

    std::deque<bool> seen_queue;
    unsigned int n_already_seen = 0;
    unsigned int n_kmers = 0;

    // first, put together an array for presence/absence of the k-mer
    // at each given position.
    while (!kmers.done()) {
        HashIntoType kmer = kmers.next();

        if (get_count(kmer)) {
            seen_queue.push_back(true);
            n_already_seen++;
        } else {
            seen_queue.push_back(false);
        }
        n_kmers++;
    }

    // next, run through this array with 'i'.

    unsigned int i = 0;
    while (i < n_kmers - min_length) {
        unsigned int seen_counter, j;

        // For each starting 'i', count the number of 'seen' k-mers in the
        // given window.

        // yes, inefficient n^2 algorithm.  sue me.
        for (seen_counter = 0, j = 0; j < min_length; j++) {
            if (seen_queue[i + j]) {
                seen_counter++;
            }
        }

        // If the fraction seen is small enough to be interesting, suggesting
        // that this, in fact, a "new" window -- extend until it isn't, and
        // then extract.

        if (!(j == min_length)) {
            throw std::exception();
        }
        if ( ((float)seen_counter / (float) j) <= max_seen) {
            unsigned int start = i;

            // extend the window until the end of the sequence...
            while ((start + min_length) < n_kmers) {
                if (seen_queue[start]) {
                    seen_counter--;
                }
                if (seen_queue[start + min_length]) {
                    seen_counter++;
                }
                start++;

                // ...or until we've seen too many of the k-mers.
                if (((float)seen_counter / (float) min_length) > max_seen) {
                    break;
                }
            }

            // adjust for ending point.
            if (start + min_length == n_kmers) {	// potentially decrement twice at end
                if (((float)seen_counter / (float) min_length) > max_seen) {
                    start--;
                }
                start--;
            } else {
                start -= 2;
            }

            // ...and now extract the relevant portion of the sequence, and adjust
            // starting pos'n.
            results.push_back(seq.substr(i, start + min_length + _ksize - i));

            i = start + min_length + 1;
        } else {
            i++;
        }
    }
}
// vim: set sts=2 sw=2:
