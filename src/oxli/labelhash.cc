/*
This file is part of khmer, https://github.com/dib-lab/khmer/, and is
Copyright (C) 2013-2015, Michigan State University.
Copyright (C) 2015-2016, The Regents of the University of California.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.

    * Neither the name of the Michigan State University nor the names
      of its contributors may be used to endorse or promote products
      derived from this software without specific prior written
      permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
LICENSE (END)

Contact: khmer-project@idyll.org
*/
#include <errno.h>
#include <string.h>
#include <iostream>
#include <sstream> // IWYU pragma: keep
#include <set>

#include "oxli/hashgraph.hh"
#include "oxli/oxli_exception.hh"
#include "oxli/labelhash.hh"
#include "oxli/read_parsers.hh"
#include "oxli/subset.hh"

#define IO_BUF_SIZE 250*1000*1000

#define LABEL_DBG 0
#define printdbg(m) if(LABEL_DBG) std::cout << #m << std::endl;

#define DEBUG 0

using namespace std;
using namespace oxli;
using namespace oxli:: read_parsers;


/*
 * @camillescott
 * Might be time for a refactor: could do a general consume_seqfile
 * function which accepts a consume_sequence function pointer as a parameter
 */

template<typename SeqIO>
void LabelHash::consume_seqfile_and_tag_with_labels(
    std:: string const  &filename,
    unsigned int	      &total_reads, unsigned long long	&n_consumed,
    CallbackFn	      callback,	    void *		callback_data
)
{
    ReadParserPtr<SeqIO> parser = get_parser<SeqIO>(filename);
    consume_seqfile_and_tag_with_labels<SeqIO>(
        parser,
        total_reads, n_consumed,
        callback, callback_data
    );
}

template<typename SeqIO>
void LabelHash::consume_seqfile_and_tag_with_labels(
    ReadParserPtr<SeqIO>& parser,
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

    Label the_label = 0;

    // Iterate through the reads and consume their k-mers.
    while (!parser->is_complete( )) {
        try {
            read = parser->get_next_read( );
        } catch (NoMoreReadsAvailable &exc) {
            break;
        }

        read.set_clean_seq();

        // TODO: make threadsafe!
        unsigned long long this_n_consumed = 0;
        consume_sequence_and_tag_with_labels( read.cleaned_seq,
                                              this_n_consumed,
                                              the_label );
        the_label++;

#if (0) // Note: Used with callback - currently disabled.
        n_consumed_LOCAL  = __sync_add_and_fetch( &n_consumed, this_n_consumed );
#else
        __sync_add_and_fetch( &n_consumed, this_n_consumed );
#endif
        __sync_add_and_fetch( &total_reads, 1 );

        // TODO: Figure out alternative to callback into Python VM
        //       Cannot use in multi-threaded operation.
#if (0)
        // run callback, if specified
        if (total_reads_TL % CALLBACK_PERIOD == 0 && callback) {
            std::cout << "n tags: " << graph->all_tags.size() << "\n";
            try {
                callback("consume_seqfile_and_tag_with_labels", callback_data,
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

template<typename SeqIO>
void LabelHash::consume_partitioned_fasta_and_tag_with_labels(
    const std::string &filename,
    unsigned int &total_reads,
    unsigned long long &n_consumed,
    CallbackFn callback,
    void * callback_data)
{
    total_reads = 0;
    n_consumed = 0;

    ReadParserPtr<SeqIO> parser = get_parser<SeqIO>(filename);
    Read read;

    std::string seq = "";

    // reset the master subset partition
    //delete partition;
    //partition = new SubsetPartition(this);

    //
    // iterate through the FASTA file & consume the reads.
    //
    PartitionID p;
    while(!parser->is_complete())  {
        read = parser->get_next_read();

        read.set_clean_seq();
        seq = read.cleaned_seq;

        // First, figure out what the partition is (if non-zero), and
        // save that.
        printdbg(parsing partition id)
        p = _parse_partition_id(read.name);
        printdbg(consuming sequence and tagging)
        consume_sequence_and_tag_with_labels( seq,
                                              n_consumed,
                                              p );
        printdbg(back in consume_partitioned)

        // reset the sequence info, increment read number
        total_reads++;

        // run callback, if specified
        if (total_reads % CALLBACK_PERIOD == 0 && callback) {
            try {
                callback("consume_partitioned_fasta_and_tag_with_labels", callback_data,
                         total_reads, n_consumed);
            } catch (...) {
                throw;
            }
        }
    }
    printdbg(done with while loop in consume_partitioned)
    printdbg(deleted parser and exiting)
}

// Note: this function assumes that 'kmer' is already in graph->all_tags;
// see usage elsewhere in this code.

void LabelHash::link_tag_and_label(const HashIntoType kmer,
                                   const Label kmer_label)
{
    printdbg(linking tag and label)
    tag_labels.insert(TagLabelPair(kmer, kmer_label));
    label_tag.insert(LabelTagPair(kmer_label, kmer));
    all_labels.insert(kmer_label);
    printdbg(done linking tag and label)
}

void LabelHash::consume_sequence_and_tag_with_labels(const std::string& seq,
        unsigned long long& n_consumed,
        Label current_label,
        SeenSet * found_tags)
{

    printdbg(inside low-level labelhash consume sequence function)

    bool kmer_tagged;

    KmerIterator kmers(seq.c_str(), graph->_ksize);
    HashIntoType kmer;

    unsigned int since = graph->_tag_density / 2 + 1;

    printdbg(entering while loop)
        while(!kmers.done()) {
            kmer = kmers.next();
            bool is_new_kmer;

            if ((is_new_kmer = graph->store->test_and_set_bits( kmer ))) {
                ++n_consumed;
                printdbg(test_and_set_bits)
            }
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
        LabelSet& found_labels,
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
            throw oxli_exception();
        }
    }
    tagged_kmers.clear();
    return num_traversed;
}

void LabelHash::get_tag_labels(const HashIntoType tag,
                               LabelSet& labels) const
{
    if (set_contains(graph->all_tags, tag)) {
        _get_tag_labels(tag, tag_labels, labels);
    }
}

// get_labels_for_sequence: return labels present in the given sequence.

void LabelHash::get_labels_for_sequence(const std::string& seq,
                                        LabelSet& found_labels)
const
{
    bool kmer_tagged;
    TagSet tags;

    KmerIterator kmers(seq.c_str(), graph->_ksize);
    HashIntoType kmer;

    while(!kmers.done()) {
        kmer = kmers.next();

        kmer_tagged = set_contains(graph->all_tags, kmer);

        if (kmer_tagged) {
            tags.insert(kmer);
        }
    }

    SeenSet::const_iterator si;
    for (si = tags.begin(); si != tags.end(); ++si) {
        _get_tag_labels(*si, tag_labels, found_labels);
    }
}

void LabelHash::get_tags_from_label(const Label label,
                                    TagSet& tags) const
{
    if(set_contains(all_labels, label)) {
        _get_tags_from_label(label, label_tag, tags);
    }
}

void LabelHash::traverse_labels_and_resolve(const SeenSet tagged_kmers,
        LabelSet& found_labels)
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
    ;
}


// Save a partition map to disk.

void LabelHash::save_labels_and_tags(std::string filename)
{
    ofstream outfile(filename.c_str(), ios::binary);

    outfile.write(SAVED_SIGNATURE, 4);
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

    TagLabelMap::const_iterator pi = tag_labels.begin();
    for (; pi != tag_labels.end(); ++pi) {
        HashIntoType *k_p = (HashIntoType *) (buf + n_bytes);
        *k_p = pi->first;
        n_bytes += sizeof(HashIntoType);

        Label * l_p = (Label *) (buf + n_bytes);
        *l_p = pi->second;
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
        throw oxli_file_exception(strerror(errno));
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
        throw oxli_file_exception(err);
    } catch (const std::exception &e) {
        // Catching std::exception is a stopgap for
        // https://gcc.gnu.org/bugzilla/show_bug.cgi?id=66145
        std::string err = "Unknown error opening file: " + filename + " "
                          + strerror(errno);
        throw oxli_file_exception(err);
    }

    unsigned long n_labeltags = 1;
    try {
        unsigned int save_ksize = 0;
        char signature[4];
        unsigned char version = 0, ht_type = 0;

        infile.read(signature, 4);
        infile.read((char *) &version, 1);
        infile.read((char *) &ht_type, 1);
        if (!(std::string(signature, 4) == SAVED_SIGNATURE)) {
            std::ostringstream err;
            err << "Incorrect file signature 0x";
            for(size_t i=0; i < 4; ++i) {
                err << std::hex << (int) signature[i];
            }
            err << " while reading labels/tags from " << filename
                << " Should be: " << SAVED_SIGNATURE;
            throw oxli_file_exception(err.str());
        } else if (!(version == SAVED_FORMAT_VERSION)) {
            std::ostringstream err;
            err << "Incorrect file format version " << (int) version
                << " while reading labels/tags from " << filename;
            throw oxli_file_exception(err.str());
        } else if (!(ht_type == SAVED_LABELSET)) {
            std::ostringstream err;
            err << "Incorrect file format type " << (int) ht_type
                << " while reading labels/tags from " << filename;
            throw oxli_file_exception(err.str());
        }

        infile.read((char *) &save_ksize, sizeof(save_ksize));
        if (!(save_ksize == graph->ksize())) {
            std::ostringstream err;
            err << "Incorrect k-mer size " << save_ksize
                << " while reading labels/tags from " << filename;
            throw oxli_file_exception(err.str());
        }

        infile.read((char *) &n_labeltags, sizeof(n_labeltags));
    } catch (std::ifstream::failure &e) {
        std::string err;
        err = "Unknown error reading header info from: " + filename;
        throw oxli_file_exception(err);
    } catch (oxli_file_exception &e) {
        throw e;
    } catch (const std::exception &e) {
        std::string err = "Unknown error opening file: " + filename + " "
                          + strerror(errno);
        throw oxli_file_exception(err);
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
        } catch (std::exception &e) {

            // We may get an exception here if we fail to read all the
            // expected bytes due to EOF -- only pass it up if we read
            // _nothing_.  Note that the while loop exits on EOF.

            if (infile.gcount() == 0) {
                delete[] buf;

                std::string err;
                err = "Unknown error reading data from: " + filename;
                throw oxli_file_exception(err);
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

            graph->all_tags.insert(*kmer_p);
            all_labels.insert(*labelp);
            link_tag_and_label(*kmer_p, *labelp);

            loaded++;
        }
        if (!(i == n_bytes)) {
            delete[] buf;
            throw oxli_file_exception("unknown error reading labels and tags");
        }
        memcpy(buf, buf + n_bytes, remainder);
    }

    if (remainder != 0) {
        delete[] buf;
        throw oxli_file_exception("unknown error reading labels and tags");
    }

    if (loaded != n_labeltags) {
        delete[] buf;
        throw oxli_file_exception("error loading labels: too few loaded");
    }

    delete[] buf;
}

// tag & label k-mers on either side of an HDN.

void LabelHash::label_across_high_degree_nodes(const char * s,
        SeenSet& high_degree_nodes,
        const Label label)
{
    KmerIterator kmers(s, graph->_ksize);

    unsigned long n = 0;

    Kmer prev_kmer = kmers.next();
    if (kmers.done()) {
        return;
    }
    Kmer kmer = kmers.next();
    if (kmers.done()) {
        return;
    }
    Kmer next_kmer = kmers.next();

    // ignore any situation where HDN is at beginning or end of sequence
    // @CTB testme :)
    while(!kmers.done()) {
        n++;
        if (n % 10000 == 0) {
            std::cout << "... label_across_hdn: " << n << "\n";
        }
        if (set_contains(high_degree_nodes, kmer)) {
            graph->add_tag(prev_kmer);
            graph->add_tag(kmer);
            graph->add_tag(next_kmer);
            link_tag_and_label(prev_kmer, label);
            link_tag_and_label(kmer, label);
            link_tag_and_label(next_kmer, label);
        }
        prev_kmer = kmer;
        kmer = next_kmer;
        next_kmer = kmers.next();
    }
}

template void LabelHash::consume_seqfile_and_tag_with_labels<FastxReader>(
    std:: string const &filename,
    unsigned int &total_reads,
    unsigned long long &n_consumed,
    CallbackFn callback,
    void * callback_data
);
template void LabelHash::consume_seqfile_and_tag_with_labels<FastxReader>(
    ReadParserPtr<FastxReader>& parser,
    unsigned int &total_reads,
    unsigned long long &n_consumed,
    CallbackFn callback,
    void * callback_data
);
template void LabelHash::consume_partitioned_fasta_and_tag_with_labels<FastxReader>(
    const std::string &filename,
    unsigned int &total_reads,
    unsigned long long &n_consumed,
    CallbackFn callback = NULL,
    void * callback_datac = NULL
);

