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

#include "hashbits.hh"
#include "hashtable.hh"
#include "khmer_exception.hh"
#include "labelhash.hh"
#include "read_parsers.hh"
#include "subset.hh"

#define IO_BUF_SIZE 250*1000*1000

#define LABEL_DBG 0
#define printdbg(m) if(LABEL_DBG) std::cout << #m << std::endl;

#define DEBUG 1

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

    Label the_label = 0;

    // Iterate through the reads and consume their k-mers.
    while (!parser->is_complete( )) {
        try {
            read = parser->get_next_read( );
        } catch (NoMoreReadsAvailable &exc) {
            break;
        }

        if (graph->check_and_normalize_read( read.sequence )) {
            // TODO: make threadsafe!
            unsigned long long this_n_consumed = 0;
            consume_sequence_and_tag_with_labels( read.sequence,
                                                  this_n_consumed,
                                                  the_label );
            the_label++;

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
    PartitionID p;
    while(!parser->is_complete())  {
        read = parser->get_next_read();
        seq = read.sequence;

        if (graph->check_and_normalize_read(seq)) {
            // First, figure out what the partition is (if non-zero), and
            // save that.
            printdbg(parsing partition id)
            p = _parse_partition_id(read.name);
            printdbg(consuming sequence and tagging)
            consume_sequence_and_tag_with_labels( seq,
                                                  n_consumed,
                                                  p );
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

            if ((is_new_kmer = graph->test_and_set_bits( kmer ))) {
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
            throw khmer_exception();
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
            throw khmer_file_exception(err.str());
        } else if (!(version == SAVED_FORMAT_VERSION)) {
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

            graph->all_tags.insert(*kmer_p);
            all_labels.insert(*labelp);
            link_tag_and_label(*kmer_p, *labelp);

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

// tag & label k-mers on either side of an HDN.

void LabelHash::label_across_high_degree_nodes(const char * s,
                                               SeenSet& high_degree_nodes,
                                               const Label label)
{
    KmerIterator kmers(s, graph->_ksize);

    unsigned long n = 0;

    Kmer prev_kmer = kmers.next();
    if (kmers.done()) { return; }
    Kmer kmer = kmers.next();
    if (kmers.done()) { return; }
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


// Starting from the given seed k-mer, assemble all maximal linear paths in
// both directions, using labels to skip over tricky bits.

std::vector<std::string> LabelHash::assemble_labeled_path(const Kmer seed_kmer)
    const
{
    std::string start_kmer = seed_kmer.get_string_rep(graph->_ksize);

#if DEBUG
    std::cout << "assemble right: " << start_kmer << std::endl;
#endif
    std::vector<std::string> fwd_paths;
    _assemble_labeled_right(start_kmer.c_str(), fwd_paths);

#if DEBUG
    std::cout << "assemble left: " << start_kmer << std::endl;
#endif
    start_kmer = _revcomp(start_kmer);
    std::vector<std::string> rev_paths;
    _assemble_labeled_right(start_kmer.c_str(), rev_paths);

#if DEBUG
    std::cout << "join right and left contigs: " << rev_paths.size() << std::endl;
#endif
    std::vector<std::string> paths;
    for (unsigned int i = 0; i < rev_paths.size(); i++) {
        for (unsigned int j = 0; j < fwd_paths.size(); j++) {
            std::string left = rev_paths[i];
            left = left.substr(graph->_ksize);
            std::string contig = _revcomp(left) + fwd_paths[j];
            paths.push_back(contig);
        }
    }

    return paths;
}

void LabelHash::_assemble_labeled_right(const char * start_kmer, std::vector<std::string>& paths)
    const
{
    const char bases[] = "ACGT";
    std::string kmer = start_kmer;
    std::string contig = kmer;
    bool found2 = false;
    SeenSet visited;

    while (1) {
        const char * base = &bases[0];
        bool found = false;
        char found_base;

        while(*base != 0) {
            std::string try_kmer = kmer.substr(1) + (char) *base;

            // a hit!
            if (graph->get_count(try_kmer.c_str())) {
                if (set_contains(visited, _hash(try_kmer.c_str(), graph->_ksize))) {
#if DEBUG
                    std::cout << "loop.\n";
#endif // DEBUG
                    base++;
                    continue;
                }
                if (found) {
                    found2 = true;
                    break;
                }
                found_base = (char) *base;
                found = true;
            }
            base++;
        }
        if (!found or found2) {
            break;
        } else {
            contig += found_base;
            kmer = kmer.substr(1) + found_base;
            found = true;
            visited.insert(_hash(kmer.c_str(), graph->_ksize));
#if DEBUG
            std::cout << "extending.\n";
#endif // DEBUG
        }
    }
    visited.clear();

    if (found2) {               // hit a HDN
#if DEBUG
        std::cout << "HDN: " << kmer.length() << "\n";
#endif // DEBUG

        Kmer path_begin(kmer.c_str(), kmer.length());

        LabelSet labels = get_tag_labels(path_begin);

#if DEBUG
        std::cout << "n labels: " << labels.size() << "\n";
#endif // DEBUG
        LabelSet::const_iterator li;
        std::vector<std::string> xpaths;
        for (li = labels.begin(); li != labels.end(); li++) {
            Label label = *li;

#if DEBUG
            std::cout << "working with " << label << "\n";
#endif // DEBUG
            xpaths.push_back(_assemble_linear_labels(kmer.c_str(),
                                                     label));
        }
#if DEBUG
        std::cout << "xpaths is: " << xpaths.size() << "\n";
#endif // DEBUG

        if (xpaths.size() == 0) {
            paths.push_back(contig);
            return;
        }
#if DEBUG
        for (unsigned int j = 0; j < xpaths.size(); j++) {
            std::string this_contig = contig;
            this_contig += xpaths[j];
            std::cout << "xpath " << j << ": " << this_contig << std::endl;
        }
#endif

        for (unsigned int j = 0; j < xpaths.size(); j++) {
            std::string this_contig = contig;
            this_contig += xpaths[j];
#if DEBUG
            std::cout << "recurse " << xpaths[0] << "\n";
#endif // DEBUG
            std::string last_kmer = this_contig.substr(this_contig.length() - kmer.length());
            const char * start_again = last_kmer.c_str();
#if DEBUG
            std::cout << "starting from " << start_again << "\n";
#endif // DEBUG
            std::vector<std::string> newpaths;

            _assemble_labeled_right(start_again, newpaths);

            if (newpaths.size() == 0) {
                paths.push_back(this_contig);
            }

            for (unsigned int i = 0; i < newpaths.size(); i++) {
                std::string xxx = newpaths[i];
                this_contig += xxx.substr(kmer.length());
                paths.push_back(this_contig);
            }
        }
    } else {
        paths.push_back(contig);
    }
}

std::string LabelHash::_assemble_linear_labels(const std::string start_kmer,
                                               const Label label)
    const
{
    const char bases[] = "ACGT";
    std::string kmer = start_kmer;
    std::string contig = "";
    bool found2 = false;

    while (1) {
        const char * base = &bases[0];
        bool found = false;
        char found_base;

#if DEBUG
        std::cout << "now at kmer " << kmer << "\n";
#endif // DEBUG

        while(*base != 0) {
            std::string try_kmer = kmer.substr(1) + (char) *base;

#if DEBUG
            std::cout << "trying " << (char) *base << "\n";
#endif // DEBUG

            // a hit!
            if (graph->get_count(try_kmer.c_str())) {
                Kmer tag(try_kmer.c_str(), try_kmer.length());
                LabelSet ls = get_tag_labels(tag);

#if DEBUG
                std::cout << "got count; now ls: " << ls.size() << "\n";
#endif // DEBUG
                if (set_contains(ls, label)) {
                    if (found) {
#if DEBUG
                        std::cout << "found 2..." << (char) *base << "\n";
#endif // DEBUG
                        found2 = true;
                        break;
                    }
#if DEBUG
                    std::cout << "found 1..." << (char) *base << "\n";
#endif // DEBUG
                    found_base = (char) *base;
                    found = true;
                }
            }
            base++;
        }
        if (!found || found2) {
            if (!found) {
#if DEBUG
                std::cout << "ending.\n";
#endif
            }
            break;
        } else {
            contig += found_base;
            kmer = kmer.substr(1) + found_base;
        }
    }
    return contig;
}
// vim: set sts=2 sw=2:
