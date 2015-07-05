//
// This file is part of khmer, https://github.com/dib-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2015. It is licensed under
// the three-clause BSD license; see LICENSE.
// Contact: khmer-project@idyll.org
//

#include "khmer.hh"
#include "hashtable.hh"
#include "read_parsers.hh"
#include "counting.hh"

#include <algorithm>
#include <sstream>
#include <errno.h>

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

    if (read.length() < _ksize) {
        return false;
    }

    for (unsigned int i = 0; i < read.length(); i++)  {
        read[ i ] &= 0xdf; // toupper - knock out the "lowercase bit"
        if (!is_valid_dna( read[ i ] )) {
            rc = false;
            break;
        }
    }

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
    unsigned int	      &total_reads, unsigned long long	&n_consumed
)
{
    IParser *	  parser =
        IParser::get_parser( filename );

    consume_fasta(
        parser,
        total_reads, n_consumed
    );

    delete parser;
}

void
Hashtable::
consume_fasta(
    read_parsers:: IParser *  parser,
    unsigned int		    &total_reads, unsigned long long  &n_consumed
)
{
    Read			  read;

    // Iterate through the reads and consume their k-mers.
    while (!parser->is_complete( )) {
        bool is_valid;
        try {
            read = parser->get_next_read( );
        } catch (NoMoreReadsAvailable) {
            break;
        }

        unsigned int this_n_consumed =
            check_and_process_read(read.sequence, is_valid);

        __sync_add_and_fetch( &n_consumed, this_n_consumed );
        __sync_add_and_fetch( &total_reads, 1 );

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

// technically, get medioid count... our "median" is always a member of the
// population.

void Hashtable::get_median_count(const std::string &s,
                                 BoundedCounterType &median,
                                 float &average,
                                 float &stddev)
{
    std::vector<BoundedCounterType> counts;
    this->get_kmer_counts(s, counts);

    if (!counts.size()) {
        throw khmer_exception("no k-mer counts for this string; too short?");
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

//
// Optimized filter function for normalize-by-median
//
bool Hashtable::median_at_least(const std::string &s,
                                unsigned int cutoff)
{
    KMerIterator kmers(s.c_str(), _ksize);
    unsigned int min_req = 0.5 + float(s.size() - _ksize + 1) / 2;
    unsigned int num_cutoff_kmers = 0;

    // first loop:
    // accumulate at least min_req worth of counts before checking to see
    // if we have enough high-abundance k-mers to indicate success.
    for (unsigned int i = 0; i < min_req; ++i) {
        HashIntoType kmer = kmers.next();
        if (this->get_count(kmer) >= cutoff) {
            ++num_cutoff_kmers;
        }
    }

    // second loop: now check to see if we pass the threshold for each k-mer.
    if (num_cutoff_kmers >= min_req) {
        return true;
    }
    while(!kmers.done()) {
        HashIntoType kmer = kmers.next();
        if (this->get_count(kmer) >= cutoff) {
            ++num_cutoff_kmers;
            if (num_cutoff_kmers >= min_req) {
                return true;
            }
        }
    }
    return false;
}

void Hashtable::save_tagset(std::string outfilename)
{
    ofstream outfile(outfilename.c_str(), ios::binary);
    const size_t tagset_size = n_tags();
    unsigned int save_ksize = _ksize;

    HashIntoType * buf = new HashIntoType[tagset_size];

    outfile.write(SAVED_SIGNATURE, 4);
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
    if (outfile.fail()) {
        delete[] buf;
        throw khmer_file_exception(strerror(errno));
    }
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
    } catch (std::ifstream::failure &e) {
        std::string err;
        if (!(infile.is_open())) {
            err = "Cannot open tagset file: " + infilename;
        } else {
            err = "Unknown error in opening file: " + infilename;
        }
        throw khmer_file_exception(err);
    }

    if (clear_tags) {
        all_tags.clear();
    }

    unsigned char version, ht_type;
    char signature[4];
    unsigned int save_ksize = 0;

    size_t tagset_size = 0;
    HashIntoType * buf = NULL;

    try {
        infile.read(signature, 4);
        infile.read((char *) &version, 1);
        infile.read((char *) &ht_type, 1);
        if (!(std::string(signature, 4) == SAVED_SIGNATURE)) {
            std::ostringstream err;
            err << "Incorrect file signature " << signature
                << " while reading tagset from " << infilename
                << "; should be " << SAVED_SIGNATURE;
            throw khmer_file_exception(err.str());
        } else if (!(version == SAVED_FORMAT_VERSION)) {
            std::ostringstream err;
            err << "Incorrect file format version " << (int) version
                << " while reading tagset from " << infilename
                << "; should be " << (int) SAVED_FORMAT_VERSION;
            throw khmer_file_exception(err.str());
        } else if (!(ht_type == SAVED_TAGS)) {
            std::ostringstream err;
            err << "Incorrect file format type " << (int) ht_type
                << " while reading tagset from " << infilename;
            throw khmer_file_exception(err.str());
        }

        infile.read((char *) &save_ksize, sizeof(save_ksize));
        if (!(save_ksize == _ksize)) {
            std::ostringstream err;
            err << "Incorrect k-mer size " << save_ksize
                << " while reading tagset from " << infilename;
            throw khmer_file_exception(err.str());
        }

        infile.read((char *) &tagset_size, sizeof(tagset_size));
        infile.read((char *) &_tag_density, sizeof(_tag_density));

        buf = new HashIntoType[tagset_size];

        infile.read((char *) buf, sizeof(HashIntoType) * tagset_size);

        for (unsigned int i = 0; i < tagset_size; i++) {
            all_tags.insert(buf[i]);
        }

        delete[] buf;
    } catch (std::ifstream::failure &e) {
        std::string err = "Error reading data from: " + infilename;
        if (buf != NULL) {
            delete[] buf;
        }
        throw khmer_file_exception(err);
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
    unsigned int	      &total_reads, unsigned long long	&n_consumed
)
{
    IParser *	  parser =
        IParser::get_parser( filename );

    consume_fasta_and_tag(
        parser,
        total_reads, n_consumed
    );

    delete parser;
}

void
Hashtable::
consume_fasta_and_tag(
    read_parsers:: IParser *  parser,
    unsigned int		    &total_reads,   unsigned long long	&n_consumed
)
{
    Read			  read;

    // TODO? Delete the following assignments.
    total_reads = 0;
    n_consumed = 0;

    // Iterate through the reads and consume their k-mers.
    while (!parser->is_complete( )) {

        try {
            read = parser->get_next_read( );
        } catch (NoMoreReadsAvailable &e) {
            // Bail out if this error is raised
            break;
        }

        if (check_and_normalize_read( read.sequence )) {
            unsigned long long this_n_consumed = 0;
            consume_sequence_and_tag( read.sequence, this_n_consumed );

            __sync_add_and_fetch( &n_consumed, this_n_consumed );
            __sync_add_and_fetch( &total_reads, 1 );
        }
    } // while reads left for parser

}

//
// consume_fasta_and_tag_with_stoptags: consume a FASTA file of reads,
//     tagging reads every so often.  Do not insert matches to stoptags,
//     and join the tags across those gaps.
//

void Hashtable::consume_fasta_and_tag_with_stoptags(const std::string &filename,
        unsigned int &total_reads,
        unsigned long long &n_consumed)
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
        try {
            read = parser->get_next_read();
        } catch (NoMoreReadsAvailable &exc) {
            break;
        }
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
        unsigned long long &n_consumed)
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
        try {
            read = parser->get_next_read();
        } catch (NoMoreReadsAvailable &exc) {
            break;
        }
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
        try {
            read = parser->get_next_read();
        } catch (NoMoreReadsAvailable &exc) {
            break;
        }
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

void Hashtable::filter_if_present(const std::string &infilename,
                                  const std::string &outputfile)
{
    IParser* parser = IParser::get_parser(infilename);
    ofstream outfile(outputfile.c_str());

    unsigned int total_reads = 0;
    unsigned int reads_kept = 0;

    Read read;
    string seq;

    HashIntoType kmer;

    while(!parser->is_complete()) {
        try {
            read = parser->get_next_read();
        } catch (NoMoreReadsAvailable &exc) {
            break;
        }
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
            throw khmer_exception();
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
            throw khmer_exception();
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
    } catch (std::ifstream::failure &e) {
        std::string err;
        if (!(infile.is_open())) {
            err = "Cannot open stoptags file: " + infilename;
        } else {
            err = "Unknown error in opening file: " + infilename;
        }
        throw khmer_file_exception(err);
    }

    if (clear_tags) {
        stop_tags.clear();
    }

    unsigned char version, ht_type;
    char signature[4];
    unsigned int save_ksize = 0;

    size_t tagset_size = 0;

    try {
        infile.read(signature, 4);
        infile.read((char *) &version, 1);
        infile.read((char *) &ht_type, 1);
        if (!(std::string(signature, 4) == SAVED_SIGNATURE)) {
            std::ostringstream err;
            err << "Incorrect file signature " << signature
                << " while reading stoptags from " << infilename
                << "; should be " << SAVED_SIGNATURE;
            throw khmer_file_exception(err.str());
        } else if (!(version == SAVED_FORMAT_VERSION)) {
            std::ostringstream err;
            err << "Incorrect file format version " << (int) version
                << " while reading stoptags from " << infilename
                << "; should be " << (int) SAVED_FORMAT_VERSION;
            throw khmer_file_exception(err.str());
        } else if (!(ht_type == SAVED_STOPTAGS)) {
            std::ostringstream err;
            err << "Incorrect file format type " << (int) ht_type
                << " while reading stoptags from " << infilename;
            throw khmer_file_exception(err.str());
        }

        infile.read((char *) &save_ksize, sizeof(save_ksize));
        if (!(save_ksize == _ksize)) {
            std::ostringstream err;
            err << "Incorrect k-mer size " << save_ksize
                << " while reading stoptags from " << infilename;
            throw khmer_file_exception(err.str());
        }
        infile.read((char *) &tagset_size, sizeof(tagset_size));

        HashIntoType * buf = new HashIntoType[tagset_size];

        infile.read((char *) buf, sizeof(HashIntoType) * tagset_size);

        for (unsigned int i = 0; i < tagset_size; i++) {
            stop_tags.insert(buf[i]);
        }
        delete[] buf;
    } catch (std::ifstream::failure &e) {
        std::string err = "Error reading stoptags from: " + infilename;
        throw khmer_file_exception(err);
    }
}

void Hashtable::save_stop_tags(std::string outfilename)
{
    ofstream outfile(outfilename.c_str(), ios::binary);
    size_t tagset_size = stop_tags.size();

    HashIntoType * buf = new HashIntoType[tagset_size];

    outfile.write(SAVED_SIGNATURE, 4);
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
            throw khmer_exception();
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


void Hashtable::get_kmers(const std::string &s,
                          std::vector<std::string> &kmers_vec) const
{
    if (s.length() < _ksize) {
        return;
    }
    for (unsigned int i = 0; i < s.length() - _ksize + 1; i++) {
        std::string sub = s.substr(i, i + _ksize);
        kmers_vec.push_back(sub);
    }
}


void Hashtable::get_kmer_hashes(const std::string &s,
                                std::vector<HashIntoType> &kmers_vec) const
{
    KMerIterator kmers(s.c_str(), _ksize);

    while(!kmers.done()) {
        HashIntoType kmer = kmers.next();
        kmers_vec.push_back(kmer);
    }
}


void Hashtable::get_kmer_counts(const std::string &s,
                                std::vector<BoundedCounterType> &counts) const
{
    KMerIterator kmers(s.c_str(), _ksize);

    while(!kmers.done()) {
        HashIntoType kmer = kmers.next();
        BoundedCounterType c = this->get_count(kmer);
        counts.push_back(c);
    }
}

// vim: set sts=2 sw=2:
