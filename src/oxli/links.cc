
#include <iostream>
#include "oxli/links.hh"

using namespace oxli;

namespace oxli {

inline std::ostream& operator<< (std::ostream& stream, const CompactNode& node) {
    stream << "<CompactNode ID=" << node.node_id << " Kmer=" << node.kmer.kmer_u
           << " Count=" << node.count << " in_degree=" << node.in_degree()
           << " out_degree=" << node.out_degree() << ">";
    return stream;
}

}
