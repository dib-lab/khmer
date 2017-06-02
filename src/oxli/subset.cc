/*
This file is part of khmer, https://github.com/dib-lab/khmer/, and is
Copyright (C) 2010-2015, Michigan State University.
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
#include <assert.h>
#include <errno.h>
#include <string.h>
#include <iostream>
#include <sstream> // IWYU pragma: keep
#include <map>
#include <set>
#include <utility>

#include "oxli/hashgraph.hh"
#include "oxli/oxli_exception.hh"
#include "oxli/read_parsers.hh"
#include "oxli/subset.hh"
#include "oxli/traversal.hh"

#define IO_BUF_SIZE 250*1000*1000
#define BIG_TRAVERSALS_ARE 200

// #define VALIDATE_PARTITIONS

using namespace oxli;
using namespace oxli:: read_parsers;
using namespace std;

#if 0

static void print_partition_set(PartitionSet& p)
{
    cout << "\tpartition set: ";
    for (PartitionSet::iterator pi = p.begin(); pi != p.end(); pi++) {
        cout << *pi << ", ";
    }
    cout << "\n";
}

static void print_tag_set(SeenSet& p)
{
    cout << "\ttag set: ";
    for (SeenSet::iterator si = p.begin(); si != p.end(); si++) {
        cout << *si << ", ";
    }
    cout << "\n";
}

#endif //0

SubsetPartition::SubsetPartition(Hashgraph * ht) :
    next_partition_id(2), _ht(ht)
{
}

void SubsetPartition::count_partitions(
    size_t& n_partitions,
    size_t& n_unassigned)
{
    n_partitions = 0;
    n_unassigned = 0;

    PartitionSet partitions;

    //
    // go through all the tagged kmers and count partitions/orphan.
    //

    for (SeenSet::iterator ti = _ht->all_tags.begin();
            ti != _ht->all_tags.end(); ++ti) {
        PartitionID * partition_p = partition_map[*ti];
        if (partition_p) {
            partitions.insert(*partition_p);
        } else {
            n_unassigned++;
        }
    }
    n_partitions = partitions.size();
}


size_t SubsetPartition::output_partitioned_file(
    const std::string	&infilename,
    const std::string	&outputfile,
    bool		output_unassigned,
    CallbackFn		callback,
    void *		callback_data)
{
    auto parser = read_parsers::get_parser<FastxReader>(infilename);
    ofstream outfile(outputfile.c_str());

    unsigned int total_reads = 0;
    unsigned int reads_kept = 0;
    size_t n_singletons = 0;

    PartitionSet partitions;

    Read read;
    string seq;

    HashIntoType kmer = 0;

    const unsigned int ksize = _ht->ksize();

    //
    // go through all the reads, and take those with assigned partitions
    // and output them.
    //

    while(!parser->is_complete()) {
        try {
            read = parser->get_next_read();
        } catch (NoMoreReadsAvailable &exc) {
            break;
        }

        read.set_clean_seq();
        seq = read.cleaned_seq;

        bool found_tag = false;
        KmerHashIteratorPtr kmers = _ht->new_kmer_iterator(read.cleaned_seq);
        while (!kmers->done()) {
            kmer = kmers->next();

            // is this a known tag?
            if (set_contains(partition_map, kmer)) {
                found_tag = true;
                break;
            }
        }

        // all sequences should have at least one tag in them.
        // assert(found_tag);  @CTB currently breaks tests.  give fn flag
        // to disable.

        PartitionID partition_id = 0;
        if (found_tag) {
            PartitionID * partition_p = partition_map[kmer];
            if (partition_p == NULL ) {
                partition_id = 0;
                n_singletons++;
            } else {
                partition_id = *partition_p;
                partitions.insert(partition_id);
            }
        }

        if (partition_id > 0 || output_unassigned) {
            if (read.quality.length()) { // FASTQ
                outfile << "@" << read.name << "\t" << partition_id
                        << "\n";
                outfile << seq << "\n+\n";
                outfile << read.quality << "\n";
            } else {		// FASTA
                outfile << ">" << read.name << "\t" << partition_id;
                outfile << "\n" << seq << "\n";
            }
        }

        total_reads++;

        // run callback, if specified
        if (total_reads % CALLBACK_PERIOD == 0 && callback) {
            try {
                callback("output_partitions", callback_data,
                         total_reads, reads_kept);
            } catch (...) {
                outfile.close();
                throw;
            }
        }
    }

    return partitions.size() + n_singletons;
}

// find_all_tags: the core of the partitioning code.  finds all tagged k-mers
//    connected to kmer_f/kmer_r in the graph.

void SubsetPartition::find_all_tags(
    Kmer start_kmer,
    SeenSet&		tagged_kmers,
    const SeenSet&	all_tags,
    bool		break_on_stop_tags,
    bool		stop_big_traversals)
{

    bool first = true;
    KmerQueue node_q;
    std::queue<unsigned int> breadth_q;

    unsigned int cur_breadth = 0;
    const unsigned int max_breadth = (2 * _ht->_tag_density) + 1;

    unsigned int total = 0;
    unsigned int nfound = 0;

    KmerSet keeper;		// keep track of traversed kmers
    KmerFilter filter = [&] (const Kmer& n) -> bool {
        return set_contains(keeper, n);
    };
    Traverser traverser(_ht, filter);

    node_q.push(start_kmer);
    breadth_q.push(0);

    while(!node_q.empty()) {

        if (stop_big_traversals && keeper.size() > BIG_TRAVERSALS_ARE) {
            tagged_kmers.clear();
            break;
        }

        Kmer node = node_q.front();
        node_q.pop();

        unsigned int breadth = breadth_q.front();
        breadth_q.pop();

        if (set_contains(keeper, node)) {
            continue;
        }

        if (break_on_stop_tags && set_contains(_ht->stop_tags, node)) {
            continue;
        }

        // keep track of seen kmers
        keeper.insert(node);
        total++;

        // Is this a kmer-to-tag, and have we put this tag in a partition
        // already? Search no further in this direction.  (This is where we
        // connect partitions.)
        if (!first && set_contains(all_tags, node)) {
            tagged_kmers.insert(node);
            continue;
        }

        if (!(breadth >= cur_breadth)) { // keep track of watermark, for
            // debugging
            throw oxli_exception("Desynchonization between traversal "
                                  "and breadth tracking. Did you forget "
                                  "to pop the node or breadth queue?");
        }
        if (breadth > cur_breadth) {
            cur_breadth = breadth;
        }

        if (breadth >= max_breadth) {
            continue;    // truncate search @CTB exit?
        }

        nfound = traverser.traverse_right(node, node_q);
        for (unsigned int i = 0; i<nfound; ++i) {
            breadth_q.push(breadth + 1);
        }

        nfound = traverser.traverse_left(node, node_q);
        for (unsigned int i = 0; i<nfound; ++i) {
            breadth_q.push(breadth + 1);
        }


        first = false;
    }
}

// Perform a breadth-first search starting from the k-mers in the given
// sequence.

unsigned int SubsetPartition::sweep_for_tags(
    const std::string&	seq,
    SeenSet&		tagged_kmers,
    const SeenSet&	all_tags,
    unsigned int	range,
    bool		break_on_stop_tags,
    bool		stop_big_traversals)
{

    KmerSet traversed_nodes;
    KmerQueue node_q;
    std::queue<unsigned int> breadth_q;

    unsigned int max_breadth = range;
    unsigned int total = 0;
    unsigned int nfound = 0;

    KmerFilter filter = [&] (const Kmer& n) -> bool {
        return set_contains(traversed_nodes, n);
    };
    Traverser traverser(_ht, filter);

    // Queue up all the sequence's k-mers at breadth zero
    // We are searching around the perimeter of the known k-mers
    KmerIterator kmers(seq.c_str(), _ht->ksize());
    while (!kmers.done()) {
        Kmer node = kmers.next();
        traversed_nodes.insert(node);

        node_q.push(node);
        breadth_q.push(0);
    }

    size_t seq_length = node_q.size() / 2;
    size_t BIG_PERIMETER_TRAVERSALS = BIG_TRAVERSALS_ARE * seq_length;

    while(!node_q.empty()) {
        // change this to a better hueristic
        if (stop_big_traversals && traversed_nodes.size() >
                BIG_PERIMETER_TRAVERSALS) {
            tagged_kmers.clear();
            break;
        }

        Kmer node = node_q.front();
        node_q.pop();

        unsigned int breadth = breadth_q.front();
        breadth_q.pop();

        // Do we want to traverse through this k-mer?  If not, skip.
        if (break_on_stop_tags && set_contains(_ht->stop_tags, node)) {
            continue;
        }

        traversed_nodes.insert(node);
        total++;

        if (set_contains(all_tags, node)) {
            tagged_kmers.insert(node);
            // if we find a tag, finish the remaining queued nodes,
            // but don't queue up any more
            // max_breadth = breadth;
            continue;
        }

        if (breadth == max_breadth) {
            continue;
        }

        // finish up nodes on the current level, but if we go beyond, end it
        // immediately; this keeps from having to look at nodes which have
        // already been queued once we lower the limit after finding a tag
        else if (breadth > max_breadth) {
            return total;
        }

        nfound = traverser.traverse_right(node, node_q);
        for (unsigned int i = 0; i<nfound; ++i) {
            breadth_q.push(breadth + 1);
        }

        nfound = traverser.traverse_left(node, node_q);
        for (unsigned int i = 0; i<nfound; ++i) {
            breadth_q.push(breadth + 1);
        }
    }

    return total;
}

// find_all_tags_truncate_on_abundance: experimental code for streaming
// partitioning.

void SubsetPartition::find_all_tags_truncate_on_abundance(
    Kmer start_kmer,
    SeenSet&		tagged_kmers,
    const SeenSet&	all_tags,
    BoundedCounterType	min_count,
    BoundedCounterType	max_count,
    bool		break_on_stop_tags,
    bool		stop_big_traversals)
{

    bool first = true;
    KmerQueue node_q;
    std::queue<unsigned int> breadth_q;

    unsigned int cur_breadth = 0;
    const unsigned int max_breadth = (2 * _ht->_tag_density) + 1;

    unsigned int total = 0;
    unsigned int nfound = 0;

    KmerSet keeper;		// keep track of traversed kmers
    KmerFilter filter = [&] (const Kmer& n) -> bool {
        return set_contains(keeper, n);
    };

    Traverser traverser(_ht, filter);

    node_q.push(start_kmer);
    breadth_q.push(0);

    while(!node_q.empty()) {
        if (stop_big_traversals && keeper.size() > BIG_TRAVERSALS_ARE) {
            tagged_kmers.clear();
            break;
        }

        Kmer node = node_q.front();
        node_q.pop();

        unsigned int breadth = breadth_q.front();
        breadth_q.pop();

        // Have we already seen this k-mer?  If so, skip.
        // NOTE: redundant, move this to before while loop
        if (set_contains(keeper, node)) {
            continue;
        }

        // Do we want to traverse through this k-mer?  If not, skip.
        if (break_on_stop_tags && set_contains(_ht->stop_tags, node)) {
            // @CTB optimize by inserting into keeper set?
            continue;
        }

        BoundedCounterType count = _ht->get_count(node);
        if (count < min_count || count > max_count) {
            continue;
        }

        // keep track of seen kmers
        keeper.insert(node);
        total++;

        // Is this a kmer-to-tag, and have we put this tag in a partition
        // already? Search no further in this direction.  (This is where we
        // connect partitions.)
        if (!first && set_contains(all_tags, node)) {
            tagged_kmers.insert(node);
            continue;
        }

        // @cswelcher Do these lines actually do anything?
        if (!(breadth >= cur_breadth)) { // keep track of watermark, for
            // debugging.
            throw oxli_exception("Desynchonization between traversal "
                                  "and breadth tracking. Did you forget "
                                  "to pop the node or breadth queue?");
        }
        if (breadth > cur_breadth) {
            cur_breadth = breadth;
        }

        if (breadth >= max_breadth) {
            continue;    // truncate search @CTB exit?
        }

        nfound = traverser.traverse_right(node, node_q);
        for (unsigned int i = 0; i<nfound; ++i) {
            breadth_q.push(breadth + 1);
        }

        nfound = traverser.traverse_left(node, node_q);
        for (unsigned int i = 0; i<nfound; ++i) {
            breadth_q.push(breadth + 1);
        }

        first = false;
    }
}

// do_partition: Main partitioning code, to find components in large graphs.

void SubsetPartition::do_partition(
    HashIntoType	first_kmer,
    HashIntoType	last_kmer,
    bool		break_on_stop_tags,
    bool		stop_big_traversals,
    CallbackFn		callback,
    void *		callback_data)
{
    unsigned int total_reads = 0;

    SeenSet tagged_kmers;
    SeenSet::const_iterator si, end;

    if (first_kmer) {
        si = _ht->all_tags.find(first_kmer);
    } else {
        si = _ht->all_tags.begin();
    }
    if (last_kmer) {
        end = _ht->all_tags.find(last_kmer);
    } else {
        end = _ht->all_tags.end();
    }

    for (; si != end; ++si) {
        total_reads++;

        Kmer kmer = _ht->build_kmer(*si);

        // find all tagged kmers within range.
        tagged_kmers.clear();
        find_all_tags(kmer, tagged_kmers, _ht->all_tags,
                      break_on_stop_tags, stop_big_traversals);

        // assign the partition ID
        assign_partition_id(kmer, tagged_kmers);

        // run callback, if specified
        if (total_reads % CALLBACK_PERIOD == 0 && callback) {
            cout << "...subset-part " << first_kmer << "-" << last_kmer << ": "
                 << total_reads << " <- " << next_partition_id << "\n";
#if 0 // @CTB
            try {
                callback("do_subset_partition/read", callback_data, total_reads,
                         next_partition_id);
            } catch (...) {
                delete parser;
                throw;
            }
#endif // 0
        }
    }
}

// do_partition with abundance: Experimental code to do streaming partitioning.

void SubsetPartition::do_partition_with_abundance(
    HashIntoType	first_kmer,
    HashIntoType	last_kmer,
    BoundedCounterType	min_count,
    BoundedCounterType	max_count,
    bool		break_on_stop_tags,
    bool		stop_big_traversals,
    CallbackFn		callback,
    void *		callback_data)
{
    unsigned int total_reads = 0;

    SeenSet tagged_kmers;
    SeenSet::const_iterator si, end;

    if (first_kmer) {
        si = _ht->all_tags.find(first_kmer);
    } else {
        si = _ht->all_tags.begin();
    }
    if (last_kmer) {
        end = _ht->all_tags.find(last_kmer);
    } else {
        end = _ht->all_tags.end();
    }

    for (; si != end; ++si) {
        total_reads++;

        Kmer kmer = _ht->build_kmer(*si);

        // find all tagged kmers within range.
        tagged_kmers.clear();
        find_all_tags_truncate_on_abundance(kmer, tagged_kmers,
                                            _ht->all_tags, min_count,
                                            max_count, break_on_stop_tags,
                                            stop_big_traversals);

        // assign the partition ID
        assign_partition_id(kmer, tagged_kmers);

        // run callback, if specified
        if (total_reads % CALLBACK_PERIOD == 0 && callback) {
            cout << "...subset-part " << first_kmer << "-" << last_kmer << ": "
                 << total_reads << " <- " << next_partition_id << "\n";
#if 0 // @CTB
            try {
                callback("do_subset_partition/read", callback_data, total_reads,
                         next_partition_id);
            } catch (...) {
                delete parser;
                throw;
            }
#endif // 0
        }
    }
}


//

void SubsetPartition::set_partition_id(
    std::string kmer_s,
    PartitionID p)
{
    HashIntoType kmer;
    if (!(kmer_s.length() >= _ht->ksize())) {
        throw oxli_exception();
    }
    kmer = _ht->hash_dna(kmer_s.c_str());

    set_partition_id(kmer, p);
}

void SubsetPartition::set_partition_id(
    HashIntoType	kmer,
    PartitionID		p)
{
    PartitionPtrSet * s = reverse_pmap[p];
    PartitionID * pp = NULL;
    if (s == NULL) {
        s = new PartitionPtrSet();
        pp = new unsigned int(p);
        s->insert(pp);
        reverse_pmap[p] = s;
    } else {
        pp = *(s->begin());
    }
    partition_map[kmer] = pp;

    if (next_partition_id <= p) {
        next_partition_id = p + 1;
    }
}

PartitionID SubsetPartition::assign_partition_id(
    HashIntoType	kmer,
    SeenSet&		tagged_kmers)

{
    PartitionID return_val = 0;

    // did we find a tagged kmer?
    if (!tagged_kmers.empty()) {
        PartitionID * pp = _join_partitions_by_tags(tagged_kmers, kmer);
        return_val = *pp;
    } else {
        partition_map.erase(kmer);
        return_val = 0;
    }

    return return_val;
}

// _join_partitions_by_tags combines the tags in 'tagged_kmers' into a single
// partition, creating or reassigning partitions as necessary.  Low level
// function!

PartitionID * SubsetPartition::_join_partitions_by_tags(
    const SeenSet&	tagged_kmers,
    const HashIntoType	kmer)
{
    SeenSet::const_iterator it = tagged_kmers.begin();
    unsigned int * this_partition_p = NULL;

    // find first assigned partition ID in tagged set
    while (it != tagged_kmers.end()) {
        this_partition_p = partition_map[*it];
        if (this_partition_p != NULL) {
            break;
        }
        ++it;
    }

    // no partition ID? allocate new!
    if (this_partition_p == NULL) {
        this_partition_p = new PartitionID(next_partition_id);
        next_partition_id++;

        PartitionPtrSet * s = new PartitionPtrSet();
        s->insert(this_partition_p);
        reverse_pmap[*this_partition_p] = s;
    }

    // reassign all partitions individually.
    it = tagged_kmers.begin();
    for (; it != tagged_kmers.end(); ++it) {
        PartitionMap::iterator pi = partition_map.find(*it);

        if (pi == partition_map.end()) { // no entry? insert.
            partition_map[*it] = this_partition_p;
        } else {
            PartitionID * pp_id = pi->second;

            if (pp_id == NULL) {	// entry is null? set;
                pi->second = this_partition_p;
            } else if (*pp_id != *this_partition_p) {
                // != entry? join partitions.
                _merge_two_partitions(this_partition_p, pp_id);
            }
        }
    }

    partition_map[kmer] = this_partition_p;

    return this_partition_p;
}

// _merge_two_partitions merges the 'merge_pp' partition into the
// 'the_pp' partition.  It does this by joining the reverse pointer
// map structures for two partitions and resetting each partition
// pointer individually.

PartitionID * SubsetPartition::_merge_two_partitions(
    PartitionID *the_pp,
    PartitionID *merge_pp)
{
    PartitionPtrSet * s = reverse_pmap[*the_pp];
    PartitionPtrSet * t = reverse_pmap[*merge_pp];

    // Choose the smaller of two sets to loop over.
    if (s->size() < t->size()) {
        PartitionPtrSet * tmp = s;
        s = t;
        t = tmp;
        PartitionID * tmp2 = the_pp;
        the_pp = merge_pp;
        merge_pp = tmp2;
    }

    // Get rid of the reverse pointer for the old partition.
    reverse_pmap.erase(*merge_pp);

    // Merge all of the elements in the to-be-replaced PartitionPtrSet
    // into the merged partition.
    for (PartitionPtrSet::iterator pi = t->begin(); pi != t->end(); ++pi) {
        PartitionID * iter_pp;
        iter_pp = *pi;

        *iter_pp = *the_pp;	// reset the partition ID to the new one.
        s->insert(iter_pp);
    }
    delete t;

    return the_pp;
}

PartitionID SubsetPartition::join_partitions(
    PartitionID orig,
    PartitionID join)
{
    if (orig == join) {
        return orig;
    }
    if (orig == 0 || join == 0) {
        return 0;
    }

    if (!set_contains(reverse_pmap, orig) ||
            !set_contains(reverse_pmap, join) ||
            reverse_pmap[orig] == NULL ||
            reverse_pmap[join] == NULL) {
        return 0;
    }

    PartitionID * orig_pp = *(reverse_pmap[orig]->begin());
    PartitionID * join_pp = *(reverse_pmap[join]->begin());

    _merge_two_partitions(orig_pp, join_pp);

    return orig;
}

PartitionID SubsetPartition::get_partition_id(std::string kmer_s)
{
    HashIntoType kmer;
    if (!(kmer_s.length() >= _ht->ksize())) {
        throw oxli_exception();
    }
    kmer = _ht->hash_dna(kmer_s.c_str());

    return get_partition_id(kmer);
}

PartitionID SubsetPartition::get_partition_id(HashIntoType kmer)
{
    if (set_contains(partition_map, kmer)) {
        PartitionID * pp = partition_map[kmer];
        if (pp == NULL) {
            return 0;
        }
        return *pp;
    }
    return 0;
}

void SubsetPartition::merge(SubsetPartition * other)
{
    if (this == other) {
        return;
    }

    PartitionPtrMap other_to_this;

    PartitionMap::const_iterator pi = other->partition_map.begin();
    for (; pi != other->partition_map.end(); ++pi) {
        if (pi->second) {
            _merge_other(pi->first, *(pi->second), other_to_this);
        }
    }
}

// Merge PartitionIDs from another SubsetPartition, based on overlapping
// tags.  Utility function for merge() and merge_from_disk().

void SubsetPartition::_merge_other(
    HashIntoType	tag,
    PartitionID		other_partition,
    PartitionPtrMap&	diskp_to_pp)
{
    if (set_contains(_ht->stop_tags, tag)) { // don't merge if it's a stop_tag
        return;
    }

    // OK.  Does our current partitionmap have this?
    PartitionID * pp_0;
    pp_0 = partition_map[tag];

    if (pp_0 == NULL) {	// No!  OK, map to new 'un.
        PartitionID * existing_pp_0 = diskp_to_pp[other_partition];

        if (existing_pp_0) {	// already seen this other_partition
            partition_map[tag] = existing_pp_0;
        } else {		// new other_partition! create a new partition.
            pp_0 = get_new_partition();

            PartitionPtrSet * pp_set = new PartitionPtrSet();
            pp_set->insert(pp_0);
            reverse_pmap[*pp_0] = pp_set;
            partition_map[tag] = pp_0;

            diskp_to_pp[other_partition] = pp_0;
        }
    } else {			// yes, we've seen this tag before...
        PartitionID * existing_pp_0 = diskp_to_pp[other_partition];

        if (existing_pp_0) {	// mapping exists.  copacetic?
            if (*pp_0 == *existing_pp_0) {
                ;			// yep! nothing to do, yay!
            } else {
                // remapping must be done... we need to merge!
                // the two partitions to merge are *pp_0 and *existing_pp_0.
                // we also need to reset existing_pp_0 in diskp_to_pp to pp_0.

                pp_0 = _merge_two_partitions(pp_0, existing_pp_0);
                diskp_to_pp[other_partition] = pp_0;
            }
        } else {
            // no, does not exist in our mapping yet.  but that's ok,
            // we can fix that.
            diskp_to_pp[other_partition] = pp_0;
        }
    }
}

// Merge an on-disk SubsetPartition into this one.

void SubsetPartition::merge_from_disk(string other_filename)
{
    ifstream infile;
    unsigned long long expected_pmap_size;

    // configure ifstream to raise exceptions for everything.
    infile.exceptions(std::ifstream::failbit | std::ifstream::badbit);

    try {
        infile.open(other_filename.c_str(), ios::binary);
    }  catch (std::ifstream::failure &e) {
        std::string err;
        if (!infile.is_open()) {
            err = "Cannot open subset pmap file: " + other_filename;
        } else {
            err = "Unknown error in opening file: " + other_filename;
        }
        throw oxli_file_exception(err);
    } catch (const std::exception &e) {
        // Catching std::exception is a stopgap for
        // https://gcc.gnu.org/bugzilla/show_bug.cgi?id=66145
        std::string err = "Unknown error opening file: " + other_filename + " "
                          + strerror(errno);
        throw oxli_file_exception(err);
    }
    infile.seekg(0, infile.end);
    const int length = infile.tellg();
    infile.seekg(0, infile.beg);
    if (length == 18) {
        std::string err;
        err = other_filename + " contains only a header and no partition IDs.";
        throw oxli_file_exception(err);
    }


    try {
        unsigned int save_ksize = 0;
        char signature[4];
        unsigned char version, ht_type;

        infile.read(signature, 4);
        infile.read((char *) &version, 1);
        infile.read((char *) &ht_type, 1);
        if (!(std::string(signature, 4) == SAVED_SIGNATURE)) {
            std::ostringstream err;
            err << "Incorrect file signature 0x";
            for(size_t i=0; i < 4; ++i) {
                err << std::hex << (int) signature[i];
            }
            err << " while reading subset pmap from " << other_filename
                << " Should be: " << SAVED_SIGNATURE;
            throw oxli_file_exception(err.str());
        } else if (!(version == SAVED_FORMAT_VERSION)) {
            std::ostringstream err;
            err << "Incorrect file format version " << (int) version
                << " while reading subset pmap from " << other_filename;
            throw oxli_file_exception(err.str());
        } else if (!(ht_type == SAVED_SUBSET)) {
            std::ostringstream err;
            err << "Incorrect file format type " << (int) ht_type
                << " while reading subset pmap from " << other_filename;
            throw oxli_file_exception(err.str());
        }

        infile.read((char *) &save_ksize, sizeof(save_ksize));
        if (!(save_ksize == _ht->ksize())) {
            std::ostringstream err;
            err << "Incorrect k-mer size " << save_ksize
                << " while reading subset pmap from " << other_filename;
            throw oxli_file_exception(err.str());
        }

        infile.read((char *) &expected_pmap_size, sizeof(expected_pmap_size));
    } catch (std::ifstream::failure &e) {
        std::string err;
        err = "Unknown error reading header info from: " + other_filename;
        throw oxli_file_exception(err);
    } catch (oxli_file_exception &e) {
        throw e;
    } catch (const std::exception &e) {
        std::string err = "Unknown error opening file: " + other_filename + " "
                          + strerror(errno);
        throw oxli_file_exception(err);
    }
    char * buf = new char[IO_BUF_SIZE];

    unsigned int loaded = 0;
    long remainder;


    PartitionPtrMap diskp_to_pp;

    HashIntoType * kmer_p = NULL;
    PartitionID * diskp = NULL;

    //
    // Run through the entire partitionmap file, figuring out what partition
    // IDs are present.
    //

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
                err = "Unknown error reading data from: " + other_filename;
                throw oxli_file_exception(err);
            }
        }

        long n_bytes = infile.gcount() + remainder;
        remainder = n_bytes % (sizeof(PartitionID) + sizeof(HashIntoType));
        n_bytes -= remainder;

        iteration++;

        for (i = 0; i < n_bytes;) {
            kmer_p = (HashIntoType *) (buf + i);
            i += sizeof(HashIntoType);
            diskp = (PartitionID *) (buf + i);
            i += sizeof(PartitionID);

            assert((*diskp != 0)); // sanity check!

            _merge_other(*kmer_p, *diskp, diskp_to_pp);

            loaded++;
        }
        assert(i == n_bytes);
        memcpy(buf, buf + n_bytes, remainder);
    }
    delete[] buf;

    if (loaded != expected_pmap_size) {
        throw oxli_file_exception("error loading partitionmap - "
                                   "invalid # of items");
    }
}

// Save a partition map to disk.

void SubsetPartition::save_partitionmap(string pmap_filename)
{
    ofstream outfile(pmap_filename.c_str(), ios::binary);

    unsigned char version = SAVED_FORMAT_VERSION;
    outfile.write(SAVED_SIGNATURE, 4);
    outfile.write((const char *) &version, 1);

    unsigned char ht_type = SAVED_SUBSET;
    outfile.write((const char *) &ht_type, 1);

    unsigned int save_ksize = _ht->ksize();
    outfile.write((const char *) &save_ksize, sizeof(save_ksize));

    unsigned long long pmap_size = partition_map.size();
    outfile.write((const char *) &pmap_size, sizeof(pmap_size));

    ///

    char * buf = NULL;
    buf = new char[IO_BUF_SIZE];
    unsigned int n_bytes = 0;

    // For each tag in the partition map, save the tag and the associated
    // partition ID.

    PartitionMap::const_iterator pi = partition_map.begin();
    for (; pi != partition_map.end(); ++pi) {
        HashIntoType kmer = pi->first;
        if (pi->second != NULL) {	// if a partition ID has been
            /// assigned... save.
            PartitionID p_id = *(pi->second);

            // each record consists of one tag followed by one PartitionID.
            HashIntoType * kmer_p = (HashIntoType *) (buf + n_bytes);
            *kmer_p = kmer;
            n_bytes += sizeof(HashIntoType);

            PartitionID * pp = (PartitionID *) (buf + n_bytes);
            *pp = p_id;
            n_bytes += sizeof(PartitionID);

            // flush to disk
            if (n_bytes >= IO_BUF_SIZE - sizeof(HashIntoType) -
                    sizeof(PartitionID)) {
                outfile.write(buf, n_bytes);
                n_bytes = 0;
            }
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

// Load a partition map from disk.

void SubsetPartition::load_partitionmap(string infilename)
{
    // @CTB make sure this is an empty partition...?
    merge_from_disk(infilename);
}


void SubsetPartition::_validate_pmap()
{
    for (PartitionMap::const_iterator pi = partition_map.begin();
            pi != partition_map.end(); ++pi) {
        //HashIntoType kmer = (*pi).first;
        PartitionID * pp_id = (*pi).second;

        if (pp_id != NULL) {
            if (!(*pp_id >= 1) || !(*pp_id < next_partition_id)) {
                throw oxli_exception();
            }
        }
    }

    for (ReversePartitionMap::const_iterator ri = reverse_pmap.begin();
            ri != reverse_pmap.end(); ++ri) {
        PartitionID p = (*ri).first;
        PartitionPtrSet *s = (*ri).second;

        if (!(s != NULL)) {
            throw oxli_exception();
        }

        for (PartitionPtrSet::const_iterator si = s->begin(); si != s->end();
                ++si) {
            PartitionID * pp;
            pp = *si;

            if (!(p == *pp)) {
                throw oxli_exception();
            }
        }
    }
}

// Get rid of all partitions & partition information.

void SubsetPartition::_clear_all_partitions()
{
    for (ReversePartitionMap::iterator ri = reverse_pmap.begin();
            ri != reverse_pmap.end(); ++ri) {
        PartitionPtrSet * s = (*ri).second;

        for (PartitionPtrSet::iterator pi = s->begin(); pi != s->end(); ++pi) {
            PartitionID * pp = (*pi);
            delete pp;
        }
        delete s;
    }
    partition_map.clear();
    next_partition_id = 1;
}


void SubsetPartition::partition_size_distribution(
    PartitionCountDistribution	&d,
    unsigned int		&n_unassigned)
const
{
    PartitionCountMap cm;

    partition_sizes(cm, n_unassigned);

    for (PartitionCountMap::iterator cmi = cm.begin(); cmi != cm.end();
            ++cmi) {
        d[cmi->second]++;
    }
}

void SubsetPartition::partition_sizes(
    PartitionCountMap	&cm,
    unsigned int	&n_unassigned)
const
{
    n_unassigned = 0;

    // @CTB: should this be all_tags? See count_partitions.
    for (PartitionMap::const_iterator pi = partition_map.begin();
            pi != partition_map.end(); ++pi) {
        if (pi->second) {
            cm[*(pi->second)]++;
        } else {
            n_unassigned++;
        }
    }
}

void SubsetPartition::partition_average_coverages(
    PartitionCountMap	&cm,
    Countgraph * ht) const
{
    PartitionCountMap csum;
    PartitionCountMap cN;

    // CTB: should *only* be members of this partition, so *not* all_tags.
    for (PartitionMap::const_iterator pi = partition_map.begin();
            pi != partition_map.end(); ++pi) {
        if (pi->second) {
            BoundedCounterType count = ht->get_count(pi->first);
            csum[*(pi->second)] += count;
            cN[*(pi->second)]++;
        }
    }

    for (PartitionCountMap::iterator pi = csum.begin();
            pi != csum.end(); ++pi) {
        cm[pi->first] = pi->second / float(cN[pi->first]);
    }
}

unsigned long long SubsetPartition::repartition_largest_partition(
    unsigned int distance,
    unsigned int threshold,
    unsigned int frequency,
    Countgraph &counting)
{
    PartitionCountMap cm;
    unsigned int n_unassigned = 0;
    PartitionID biggest_p = 0;
    unsigned long long next_largest = 0;

#if VERBOSE_REPARTITION
    std::cout << "calculating partition size distribution.\n";
#endif // 0

    // first, count the number of members in each partition.
    for (PartitionMap::const_iterator pi = partition_map.begin();
            pi != partition_map.end(); ++pi) {
        if (pi->second) {
            cm[*(pi->second)]++;
        } else {
            n_unassigned++;
        }
    }

    // then, build the distribution.
    PartitionCountDistribution d;

    for (PartitionCountMap::const_iterator cmi = cm.begin(); cmi != cm.end();
            ++cmi) {
        d[cmi->second]++;
    }

    // find biggest.
    PartitionCountDistribution::const_iterator di = d.end();

    if (d.empty()) {
        throw oxli_exception();
    }
    --di;

    for (PartitionCountMap::const_iterator cmi = cm.begin(); cmi != cm.end();
            ++cmi) {
        if (cmi->second == di->first) {
            biggest_p = cmi->first;	// find PID of largest partition
        }
    }
    if (!(biggest_p != 0)) {
        throw oxli_exception();
    }

#if VERBOSE_REPARTITION
    std::cout << "biggest partition: " << di->first << "\n";
#endif // 0

#if VERBOSE_REPARTITION
    std::cout << "biggest partition ID: " << biggest_p << "\n";
#endif // 0

    if (di != d.begin()) {
        --di;
        next_largest = di->first;
    }

#if VERBOSE_REPARTITION
    std::cout << "next biggest partition: " << di->first << "\n";
#endif // 0

    ///

    SeenSet bigtags;
    _clear_partition(biggest_p, bigtags);
#if VERBOSE_REPARTITION
    std::cout << "gathered/cleared " << bigtags.size() << " tags.\n";
#endif // 0

    /// Now, go through and traverse from all the bigtags, tracking
    // those that lead to well-connected sets.

    unsigned int i = 0;
    unsigned int n = 0;
    unsigned int count;
    unsigned int n_big = 0;
    KmerSet keeper;

    SeenSet::const_iterator si = bigtags.begin();

    for (; si != bigtags.end(); ++si, i++) {
        n++;

#if 1
        if (set_contains(_ht->repart_small_tags, *si)) {
            continue;
        }
#endif //0

        count = _ht->traverse_from_kmer(_ht->build_kmer(*si),
                                        distance, keeper);

        if (count >= threshold) {
            n_big++;

            KmerSet::const_iterator ti;
            for (ti = keeper.begin(); ti != keeper.end(); ++ti) {
                if (counting.get_count(*ti) > frequency) {
                    _ht->stop_tags.insert((*ti).kmer_u);
                } else {
                    counting.count(*ti);
                }
            }
#if VERBOSE_REPARTITION
            std::cout << "traversed from " << n << " tags total, of "
                      << bigtags.size() << "; "
                      << n_big << " big; size is " << keeper.size()
                      << "; " << _ht->repart_small_tags.size() << " small\n";
#endif // 0
        } else {
#if 1
            _ht->repart_small_tags.insert(*si);
#endif //0
        }
        keeper.clear();

#if VERBOSE_REPARTITION
        if (n % 1000 == 0) {
            std::cout << "found big 'un!  traversed " << n << " tags, " <<
                      n_big << " big; " << bigtags.size() << " total tags; " <<
                      _ht->stop_tags.size() << " stop tags\n";
        }
#endif // 0
    }

    // return next_largest;
#if VERBOSE_REPARTITION
    std::cout << "repartitioning...\n";
#endif // 0
    repartition_a_partition(bigtags);

    //

    return next_largest;
}

void SubsetPartition::repartition_a_partition(const SeenSet& partition_tags)
{
    SeenSet tagged_kmers;
    SeenSet::const_iterator si;

    unsigned n = 0;
    for (si = partition_tags.begin(); si != partition_tags.end(); ++si, n++) {
        if (n % 1000 == 0) {
#if VERBOSE_REPARTITION
            std::cout << "repartitioning... on " << n << " of " <<
                      partition_tags.size() << "\n";
#endif // 0
        }

        Kmer kmer = _ht->build_kmer(*si);

        tagged_kmers.clear();
        find_all_tags(kmer, tagged_kmers, _ht->all_tags, true, false);

        // only join things already in bigtags.
        SeenSet::iterator ssi = tagged_kmers.begin();
        while (ssi != tagged_kmers.end()) {
            if (!set_contains(partition_tags, *ssi)) {
                tagged_kmers.erase(ssi++);
            } else {
                ++ssi;
            }
        }
        // std::cout << "joining: " << tagged_kmers.size() << "\n";
        assign_partition_id(kmer, tagged_kmers);
    }
}

// _clear_partition: given a partition ID, identifies all tags that belong
//    to that partition & (a) clears their PID, and (b) adds them to
//    the SeenSet partition_tags.  partition_tags is cleared first.

void SubsetPartition::_clear_partition(
    PartitionID	the_partition,
    SeenSet&	partition_tags)
{
    partition_tags.clear();

    for (PartitionMap::iterator pi = partition_map.begin();
            pi != partition_map.end(); ++pi) {
        if (pi->second && *(pi->second) == the_partition) {
            partition_tags.insert(pi->first);
        }
    }

    for (SeenSet::const_iterator si = partition_tags.begin();
            si != partition_tags.end(); ++si) {
        partition_map.erase(*si);
    }

    // clear out the reverse partition mapping, too.
    PartitionPtrSet * ps = reverse_pmap[the_partition];
    for (PartitionPtrSet::iterator psi = ps->begin(); psi != ps->end();
            ++psi) {
        delete *psi;
    }
    delete ps;

    reverse_pmap.erase(the_partition);
}

void SubsetPartition::report_on_partitions()
{
    std::cout << _ht->all_tags.size() << " tags total\n";
    std::cout << reverse_pmap.size() << " partitions total\n";

    for (SeenSet::iterator ti = _ht->all_tags.begin();
            ti != _ht->all_tags.end(); ++ti) {
        std::cout << "TAG: " << _ht->unhash_dna(*ti) << "\n";
        PartitionID *pid = partition_map[*ti];
        if (pid) {
            std::cout << "partition: " << *(partition_map[*ti]) << "\n";
        } else {
            std::cout << "NULL.\n";
        }
        std::cout << "--\n";
    }
}
