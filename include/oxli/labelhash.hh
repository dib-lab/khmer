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
#ifndef LABELHASH_HH
#define LABELHASH_HH

#include <stddef.h>
#include <stdint.h>
#include <map>
#include <string>
#include <utility>

#include "hashgraph.hh"
#include "oxli.hh"
#include "read_parsers.hh"

namespace oxli
{
    class Hashgraph;

    using read_parsers::ReadParser;
    namespace read_parsers
    {
        template<typename SeqIO> class ReadParser;
        class FastxReader;
    }

class LabelHash
{
protected:
    // Does the given tag already have the given label?
    bool _cmap_contains_label(const TagLabelMap& cmap,
                              const HashIntoType kmer,
                              const Label the_label)
    {
        std::pair<TagLabelMap::const_iterator, TagLabelMap::const_iterator> ret;
        ret = cmap.equal_range(kmer);
        for (TagLabelMap::const_iterator it=ret.first; it!=ret.second; ++it) {
            if (it->second == the_label) {
                return true;
            }
        }
        return false;
    }

    // Does the given label already have a tag associated with it?
    bool _cmap_contains_tag(const LabelTagMap& cmap,
                            Label& the_label,
                            HashIntoType& kmer)
    {
        std::pair<LabelTagMap::const_iterator, LabelTagMap::const_iterator> ret;
        ret = cmap.equal_range(the_label);
        for (LabelTagMap::const_iterator it=ret.first; it!=ret.second; ++it) {
            if(it->second == kmer) {
                return true;
            }
        }
        return false;
    }

    unsigned int _get_tag_labels(const HashIntoType tag,
                                 const TagLabelMap cmap,
                                 LabelSet& found_labels) const
    {
        unsigned int num_labels = 0;
        std::pair<TagLabelMap::const_iterator, TagLabelMap::const_iterator> ret;
        ret = cmap.equal_range(tag);
        for (TagLabelMap::const_iterator it=ret.first; it!=ret.second; ++it) {
            found_labels.insert(it->second);
            ++num_labels;
        }
        return num_labels;
    }

    unsigned int _get_tags_from_label(const Label& label,
                                      const LabelTagMap& cmap,
                                      TagSet& labeled_tags) const
    {
        unsigned int num_tags = 0;
        std::pair<LabelTagMap::const_iterator, LabelTagMap::const_iterator> ret;
        ret = cmap.equal_range(label);
        for (LabelTagMap::const_iterator it=ret.first; it!=ret.second; ++it) {
            labeled_tags.insert(it->second);
            ++num_tags;
        }
        return num_tags;
    }

    uint32_t _tag_labels_spin_lock;

public:
    oxli::Hashgraph * graph;

    explicit LabelHash(Hashgraph * ht) : graph(ht)
    {
        _tag_labels_spin_lock = 0;

    }

    ~LabelHash();

    TagLabelMap tag_labels;
    LabelTagMap label_tag;
    LabelSet all_labels;

    size_t n_labels() const
    {
        return all_labels.size();
    }

    template<typename SeqIO>
    void consume_seqfile_and_tag_with_labels(
        std::string const	  &filename,
        unsigned int	  &total_reads,
        unsigned long long  &n_consumed,
        CallbackFn	  callback	  = NULL,
        void *		  callback_data	  = NULL);

    template<typename SeqIO>
    void consume_seqfile_and_tag_with_labels(
        read_parsers::ReadParserPtr<SeqIO>& parser,
        unsigned int	    &total_reads,
        unsigned long long  &n_consumed,
        CallbackFn	    callback	    = NULL,
        void *		    callback_data   = NULL);

    template<typename SeqIO>
    void consume_partitioned_fasta_and_tag_with_labels(
        const std::string &filename,
        unsigned int &total_reads,
        unsigned long long &n_consumed,
        CallbackFn callback = NULL,
        void * callback_datac = NULL);

    void consume_sequence_and_tag_with_labels(const std::string& seq,
            unsigned long long& n_consumed,
            Label current_label,
            SeenSet * new_tags = 0);

    void get_labels_for_sequence(const std::string& seq,
                                 LabelSet& labels) const;

    void get_tag_labels(const HashIntoType tag,
                        LabelSet& labels) const;
    void get_tags_from_label(const Label label,
                             TagSet& tags) const;

    void link_tag_and_label(const HashIntoType kmer, const Label label);

    unsigned int sweep_label_neighborhood(const std::string & seq,
                                          LabelSet& found_labels,
                                          unsigned int range,
                                          bool break_on_stoptags,
                                          bool stop_big_traversals);

    void traverse_labels_and_resolve(const SeenSet tagged_kmers,
                                     LabelSet& found_labels);

    void save_labels_and_tags(std::string);
    void load_labels_and_tags(std::string);

    void label_across_high_degree_nodes(const char * sequence,
                                        SeenSet& high_degree_nodes,
                                        const Label label);
};

} // namespace oxli

#define ACQUIRE_TAG_COLORS_SPIN_LOCK \
  while(!__sync_bool_compare_and_swap( &_tag_labels_spin_lock, 0, 1));

#define RELEASE_TAG_COLORS_SPIN_LOCK \
  __sync_bool_compare_and_swap( &_tag_labels_spin_lock, 1, 0);

#endif
