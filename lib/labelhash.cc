//
// This file is part of khmer, https://github.com/dib-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2015. It is licensed under
// the three-clause BSD license; see LICENSE.
// Contact: khmer-project@idyll.org
//

#include "labelhash.hh"

#include <sstream>
#include <errno.h>

#define IO_BUF_SIZE 250*1000*1000

#define LABEL_DBG 0
#define printdbg(m) if(LABEL_DBG) std::cout << #m << std::endl;

using namespace std;
using namespace khmer;
using namespace khmer:: read_parsers;

/*
 * @camillescott
 * Might be time for a refactor: could do a general consume_fasta
 * function which accepts a consume_sequence function pointer as a parameter
 */

void
LabelHash::consume_fasta_and_tag_with_labels(
    std:: string const  &filename,
    unsigned int	      &total_reads, unsigned long long	&n_consumed,
    CallbackFn	      callback,	    void *		callback_data
)
{
    IParser *	  parser =
        IParser::get_parser( filename );

    consume_fasta_and_tag_with_labels(
        parser,
        total_reads, n_consumed,
        callback, callback_data
    );

    delete parser;
}

void
LabelHash::consume_fasta_and_tag_with_labels(
    read_parsers:: IParser *  parser,
    unsigned int		    &total_reads,   unsigned long long	&n_consumed,
    CallbackFn		    callback,	    void *		callback_data
)
{
#if (0) // Note: Used with callback - currently disabled.
    unsigned long long int  n_consumed_LOCAL	= 0;
#endif
    Read			  read;

    // TODO? Delete the following assignments.
    total_reads = 0;
    n_consumed = 0;

    Label _tag_label = 0;

    Label * the_label;
    // Iterate through the reads and consume their k-mers.
    while (!parser->is_complete( )) {
        read = parser->get_next_read( );

        if (graph->check_and_normalize_read( read.sequence )) {
            // TODO: make threadsafe!
            unsigned long long this_n_consumed = 0;
            the_label = check_and_allocate_label(_tag_label);
            consume_sequence_and_tag_with_labels( read.sequence,
                                                  this_n_consumed,
                                                  *the_label );
            _tag_label++;

#if (0) // Note: Used with callback - currently disabled.
            n_consumed_LOCAL  = __sync_add_and_fetch( &n_consumed, this_n_consumed );
#else
            __sync_add_and_fetch( &n_consumed, this_n_consumed );
#endif
            __sync_add_and_fetch( &total_reads, 1 );
        }

        // TODO: Figure out alternative to callback into Python VM
        //       Cannot use in multi-threaded operation.
#if (0)
        // run callback, if specified
        if (total_reads_TL % CALLBACK_PERIOD == 0 && callback) {
            std::cout << "n tags: " << graph->all_tags.size() << "\n";
            try {
                callback("consume_fasta_and_tag_with_labels", callback_data,
                         total_reads_TL,
                         n_consumed);
            } catch (...) {
                delete parser;
                throw;
            }
        }
#endif // 0

    } // while reads left for parser

}

void LabelHash::consume_partitioned_fasta_and_tag_with_labels(
    const std::string &filename,
    unsigned int &total_reads,
    unsigned long long &n_consumed,
    CallbackFn callback,
    void * callback_data)
{
    total_reads = 0;
    n_consumed = 0;

    IParser* parser = IParser::get_parser(filename.c_str());
    Read read;

    std::string seq = "";

    // reset the master subset partition
    //delete partition;
    //partition = new SubsetPartition(this);

    //
    // iterate through the FASTA file & consume the reads.
    //
    Label * c;
    PartitionID p;
    while(!parser->is_complete())  {
        read = parser->get_next_read();
        seq = read.sequence;

        if (graph->check_and_normalize_read(seq)) {
            // First, figure out what the partition is (if non-zero), and
            // save that.
            printdbg(parsing partition id)
            p = _parse_partition_id(read.name);
            printdbg(checking label and allocating if necessary) {
                c = check_and_allocate_label(p);
            }
            printdbg(consuming sequence and tagging)
            consume_sequence_and_tag_with_labels( seq,
                                                  n_consumed,
                                                  *c );
            printdbg(back in consume_partitioned)
        }

        // reset the sequence info, increment read number
        total_reads++;

        // run callback, if specified
        if (total_reads % CALLBACK_PERIOD == 0 && callback) {
            try {
                callback("consume_partitioned_fasta_and_tag_with_labels", callback_data,
                         total_reads, n_consumed);
            } catch (...) {
                delete parser;
                throw;
            }
        }
    }
    printdbg(done with while loop in consume_partitioned)

        // @cswelcher TODO: check that deallocate LabelPtrMap is correct
    {
        delete parser;
    }
    printdbg(deleted parser and exiting)
}

// @cswelcher: double-check -- is it valid to pull the address from a reference?
void LabelHash::link_tag_and_label(HashIntoType& kmer, Label& kmer_label)
{
    printdbg(linking tag and label)
    tag_labels.insert(TagLabelPtrPair(kmer, &kmer_label));
    label_tag_ptrs.insert(LabelTagPair(kmer_label, kmer));
    printdbg(done linking tag and label)
}

void LabelHash::consume_sequence_and_tag_with_labels(const std::string& seq,
        unsigned long long& n_consumed,
        Label& current_label,
        SeenSet * found_tags)
{

    printdbg(inside low-level labelhash consume sequence function)

    bool kmer_tagged;

    KMerIterator kmers(seq.c_str(), graph->_ksize);
    HashIntoType kmer;

    unsigned int since = graph->_tag_density / 2 + 1;

    printdbg(entering while loop)
        while(!kmers.done()) {
            kmer = kmers.next();
            bool is_new_kmer;

            if ((is_new_kmer = graph->test_and_set_bits( kmer ))) {
                ++n_consumed;
                printdbg(test_and_set_bits)
            }
#if (1)
            if (is_new_kmer) {
                printdbg(new kmer...)
                ++since;
            } else {
                printdbg(entering tag spin lock)
                //ACQUIRE_ALL_TAGS_SPIN_LOCK
                kmer_tagged = set_contains(graph->all_tags, kmer);
                //RELEASE_ALL_TAGS_SPIN_LOCK
                printdbg(released tag spin lock)
                if (kmer_tagged) {
                    since = 1;
                    printdbg(kmer already in all_tags)
                    // Labeling code
                    // TODO: MAKE THREADSAFE!

                    if (!_cmap_contains_label(tag_labels, kmer, current_label)) {
                        printdbg(tag was not labeled: adding to labels...)
                        //ACQUIRE_TAG_COLORS_SPIN_LOCK
                        link_tag_and_label(kmer, current_label);
                        //RELEASE_TAG_COLORS_SPIN_LOCK
                        printdbg(released label spin lock)
                    }
                    if (found_tags) {
                        found_tags->insert(kmer);
                    }
                }  else {
                    printdbg(inc since var)
                    ++since;
                }
            }
#else
            if (!is_new_kmer && set_contains(graph->all_tags, kmer)) {
                since = 1;
                if (found_tags) {
                    found_tags->insert(kmer);
                }
            } else {
                since++;
            }
#endif
            //
            if (since >= graph->_tag_density) {
                printdbg(exceeded tag density: drop a tag and label --
                         getting tag lock)
                //ACQUIRE_ALL_TAGS_SPIN_LOCK
                printdbg(in tag spin lock)
                graph->all_tags.insert(kmer);
                //RELEASE_ALL_TAGS_SPIN_LOCK
                printdbg(released tag spin lock)

                // Labeling code
                // TODO: MAKE THREADSAFE!
                //ACQUIRE_TAG_COLORS_SPIN_LOCK
                link_tag_and_label(kmer, current_label);
                //RELEASE_TAG_COLORS_SPIN_LOCK

                if (found_tags) {
                    found_tags->insert(kmer);
                }
                since = 1;
            }
            printdbg(moving to next iter)
        } // iteration over kmers
    printdbg(finished iteration: dropping last tag)
    if (since >= graph->_tag_density/2 - 1) {
        //ACQUIRE_ALL_TAGS_SPIN_LOCK
        graph->all_tags.insert(kmer);	// insert the last k-mer, too.
        //RELEASE_ALL_TAGS_SPIN_LOCK

        // Label code: TODO: MAKE THREADSAFE!
        link_tag_and_label(kmer, current_label);

        if (found_tags) {
            found_tags->insert(kmer);
        }
    }
    printdbg(done with low-level consume)
}

unsigned int LabelHash::sweep_label_neighborhood(const std::string& seq,
        LabelPtrSet& found_labels,
        unsigned int range,
        bool break_on_stoptags,
        bool stop_big_traversals)
{

    SeenSet tagged_kmers;
    unsigned int num_traversed;
    num_traversed = graph->partition->sweep_for_tags(seq, tagged_kmers,
                    graph->all_tags,
                    range, break_on_stoptags, stop_big_traversals);
    traverse_labels_and_resolve(tagged_kmers, found_labels);
    //printf("range=%u ", range);
    if (range == 0) {
        if (!(num_traversed == seq.length()-graph->ksize()+1)) {
            throw khmer_exception();
        }
    }
    tagged_kmers.clear();
    return num_traversed;
}

LabelPtrSet LabelHash::get_tag_labels(const HashIntoType& tag)
{
    LabelPtrSet labels;
    //unsigned int num_labels;
    _get_tag_labels(tag, tag_labels, labels);
    return labels;
}

void LabelHash::traverse_labels_and_resolve(const SeenSet& tagged_kmers,
        LabelPtrSet& found_labels)
{

    SeenSet::const_iterator si;
    for (si=tagged_kmers.begin(); si!=tagged_kmers.end(); ++si) {
        HashIntoType tag = *si;
        // get the labels associated with this tag
        unsigned int num_labels = _get_tag_labels(tag, tag_labels, found_labels);
        if (num_labels > 1) {
            // reconcile labels
            // for now do nothing ha
        }
    }
}

LabelHash::~LabelHash()
{
    for (LabelPtrMap::iterator itr=label_ptrs.begin();
            itr!=label_ptrs.end(); ++itr) {
        delete itr->second;
    }
}


// Save a partition map to disk.

void LabelHash::save_labels_and_tags(std::string filename)
{
    ofstream outfile(filename.c_str(), ios::binary);

    unsigned char version = SAVED_FORMAT_VERSION;
    outfile.write((const char *) &version, 1);

    unsigned char ht_type = SAVED_LABELSET;
    outfile.write((const char *) &ht_type, 1);

    unsigned int save_ksize = graph->ksize();
    outfile.write((const char *) &save_ksize, sizeof(save_ksize));

    unsigned long n_labeltags = tag_labels.size();
    outfile.write((const char *) &n_labeltags, sizeof(n_labeltags));

    ///

    char * buf = NULL;
    buf = new char[IO_BUF_SIZE];
    unsigned int n_bytes = 0;

    // For each tag in the partition map, save the tag and the associated
    // partition ID.

    TagLabelPtrMap::const_iterator pi = tag_labels.begin();
    for (; pi != tag_labels.end(); ++pi) {
        HashIntoType *k_p = (HashIntoType *) (buf + n_bytes);
        *k_p = pi->first;
        n_bytes += sizeof(HashIntoType);

        Label * l_p = (Label *) (buf + n_bytes);
        *l_p = *(pi->second);
        n_bytes += sizeof(Label);

        // flush to disk
        if (n_bytes >= IO_BUF_SIZE - sizeof(HashIntoType) - sizeof(Label)) {
            outfile.write(buf, n_bytes);
            n_bytes = 0;
        }
    }
    // save remainder.
    if (n_bytes) {
        outfile.write(buf, n_bytes);
    }

    if (outfile.fail()) {
        delete[] buf;
        throw khmer_file_exception(strerror(errno));
    }
    outfile.close();

    delete[] buf;
}

void LabelHash::load_labels_and_tags(std::string filename)
{
    ifstream infile;

    // configure ifstream to raise exceptions for everything.
    infile.exceptions(std::ifstream::failbit | std::ifstream::badbit);

    try {
        infile.open(filename.c_str(), ios::binary);
    }  catch (std::ifstream::failure &e) {
        std::string err;
        if (!infile.is_open()) {
            err = "Cannot open labels/tags file: " + filename;
        } else {
            err = "Unknown error in opening file: " + filename;
        }
        throw khmer_file_exception(err);
    }

    unsigned long n_labeltags = 1;
    try {
        unsigned int save_ksize = 0;
        unsigned char version = 0, ht_type = 0;

        infile.read((char *) &version, 1);
        infile.read((char *) &ht_type, 1);
        if (!(version == SAVED_FORMAT_VERSION)) {
            std::ostringstream err;
            err << "Incorrect file format version " << (int) version
                << " while reading labels/tags from " << filename;
            throw khmer_file_exception(err.str());
        } else if (!(ht_type == SAVED_LABELSET)) {
            std::ostringstream err;
            err << "Incorrect file format type " << (int) ht_type
                << " while reading labels/tags from " << filename;
            throw khmer_file_exception(err.str());
        }

        infile.read((char *) &save_ksize, sizeof(save_ksize));
        if (!(save_ksize == graph->ksize())) {
            std::ostringstream err;
            err << "Incorrect k-mer size " << save_ksize
                << " while reading labels/tags from " << filename;
            throw khmer_file_exception(err.str());
        }

        infile.read((char *) &n_labeltags, sizeof(n_labeltags));
    } catch (std::ifstream::failure &e) {
        std::string err;
        err = "Unknown error reading header info from: " + filename;
        throw khmer_file_exception(err);
    }

    char * buf = new char[IO_BUF_SIZE];

    unsigned long loaded = 0;
    long remainder;


    HashIntoType * kmer_p = NULL;
    Label * labelp = NULL;

    remainder = 0;
    unsigned int iteration = 0;
    while (!infile.eof()) {
        unsigned int i;

        try {
            infile.read(buf + remainder, IO_BUF_SIZE - remainder);
        } catch (std::ifstream::failure &e) {

            // We may get an exception here if we fail to read all the
            // expected bytes due to EOF -- only pass it up if we read
            // _nothing_.  Note that the while loop exits on EOF.

            if (infile.gcount() == 0) {
                delete[] buf;

                std::string err;
                err = "Unknown error reading data from: " + filename;
                throw khmer_file_exception(err);
            }
        }

        long n_bytes = infile.gcount() + remainder;
        remainder = n_bytes % (sizeof(Label) + sizeof(HashIntoType));
        n_bytes -= remainder;

        iteration++;

        for (i = 0; i < n_bytes;) {
            kmer_p = (HashIntoType *) (buf + i);
            i += sizeof(HashIntoType);

            labelp = (Label *) (buf + i);
            i += sizeof(Label);

            Label * labelp2;

            graph->all_tags.insert(*kmer_p);
            labelp2 = check_and_allocate_label(*labelp);
            link_tag_and_label(*kmer_p, *labelp2);

            loaded++;
        }
        if (!(i == n_bytes)) {
            delete[] buf;
            throw khmer_file_exception("unknown error reading labels and tags");
        }
        memcpy(buf, buf + n_bytes, remainder);
    }

    if (remainder != 0) {
        delete[] buf;
        throw khmer_file_exception("unknown error reading labels and tags");
    }

    if (loaded != n_labeltags) {
        delete[] buf;
        throw khmer_file_exception("error loading labels: too few loaded");
    }

    delete[] buf;
}
