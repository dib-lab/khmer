/*
This file is part of khmer, https://github.com/dib-lab/khmer/, and is
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
#include <algorithm>

#include "oxli/oxli.hh"
#include "oxli/hashtable.hh"
#include "oxli/labelhash.hh"
#include "oxli/kmer_filters.hh"


namespace oxli
{

bool apply_kmer_filters(const Kmer& node, const std::list<KmerFilter>& filters)
{
    if (!filters.size()) {
        return false;
    }

    for(auto filter : filters) {
        if (filter(node)) {
            return true;
        }
    }

    return false;
}


KmerFilter get_label_filter(const Label label, const LabelHash * lh)
{
    KmerFilter filter = [=] (const Kmer& node) {
        LabelSet ls;
        lh->get_tag_labels(node, ls);
#if DEBUG_FILTERS
        if (ls.size() == 0) {
            // this should never happen
            std::cout << "no labels to jump to!" << std::endl;
        }
#endif

        return !set_contains(ls, label);

    };

    return filter;
}


KmerFilter get_simple_label_intersect_filter(const LabelSet& src_labels,
        const LabelHash * lh,
        const unsigned int min_cov)
{
    auto src_begin = src_labels.begin();
    auto src_end = src_labels.end();
    unsigned int src_size = src_labels.size();

    KmerFilter filter = [=] (const Kmer& node) {
        LabelSet dst_labels;
        lh->get_tag_labels(node, dst_labels);

        LabelSet intersect;
        std::set_intersection(src_begin, src_end,
                              dst_labels.begin(), dst_labels.end(),
                              std::inserter(intersect, intersect.begin()));

        if ((intersect.size() == 1)
                && (dst_labels.size() == 1)
                && (src_size >= min_cov)) {
#if DEBUG_FILTERS
            std::cout << "TIP: " << intersect.size() << ", " <<
                      dst_labels.size() << ", " << src_size << std::endl;
#endif
            // putative error / tip
            return true;
        } else if (intersect.size() > 0) {
            // there's at least one spanning read
            return false;
        } else {
            return true;
        }
    };

    return filter;
}


KmerFilter get_junction_count_filter(const Kmer& src_node,
                                     Countgraph * junctions,
                                     const unsigned int min_cov)
{
    KmerFilter filter = [=] (const Kmer& dst_node) {
        unsigned int jc = junctions->get_count(src_node.kmer_u ^ dst_node.kmer_u);
#if DEBUG_FILTERS
        std::cout << "Junction Count: " << jc << std::endl;
#endif
        return jc < min_cov;
    };

    return filter;
}


KmerFilter get_stop_bf_filter(const Hashtable * stop_bf)
{
    KmerFilter filter = [=] (const Kmer& n) {
        return stop_bf->get_count(n);
    };
    return filter;
}


KmerFilter get_visited_filter(std::shared_ptr<SeenSet> visited)
{
#if DEBUG_FILTERS
    std::cout << "Create new visited filter with " << visited <<
              " containing " << visited->size() << " nodes" << std::endl;
#endif
    KmerFilter filter = [=] (const Kmer& node) {
        return set_contains(*visited, node);
    };
    return filter;
}

}
