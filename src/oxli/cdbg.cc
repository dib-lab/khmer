
#include <iostream>
#include "oxli/cdbg.hh"

using namespace oxli;

namespace oxli {

void CompactEdgeFactory::write_gml(const std::string filename,
                                   const CompactNodeFactory& nodes) const {

    std::ofstream file;
    file.open(filename);
    pdebug("opened " << filename);
    file << "graph" << std::endl << "[" << std::endl;

    pdebug("writing " << nodes.n_nodes() << " nodes");
    for (auto node : nodes.compact_nodes) {
        file << "  node [" << std::endl;
        file << "    id " << std::to_string(node.node_id) << std::endl;
        file << "    kmer \"" << node.sequence << "\"" << std::endl;
        file << "    count \"" << std::to_string(node.count) << "\"" << std::endl;
        file << "  ]" << std::endl;
    }

    uint32_t edge_offset = INT_MAX / 2;
    pdebug("writing " << compact_edges.size() << " edges");
    for (auto edge_pair : compact_edges) {
        
        id_t edge_id = edge_pair.first + edge_offset;
        CompactEdge* edge = edge_pair.second;
        
        file << "  edge [" << std::endl;
        file << "    id " << std::to_string(edge_id) << std::endl;

        id_t in_id, out_id;
        bool in_null = false, out_null = false;
        if (edge->in_node_id == NULL_ID) {
            in_id = INT_MAX - edge_id;
            in_null = true;
        } else {
            in_id = edge->in_node_id;
        }
        if(edge->out_node_id == NULL_ID) {
            out_id = INT_MAX - edge_id;
            out_null = true;
        } else {
            out_id = edge->out_node_id;
        }

        if (in_null && out_null) {
            std::cerr << "in and out nodes NULL_ID, something weird with "
                << edge->edge_id << std::endl;
        }

        file << "    source " << std::to_string(in_id) << std::endl;
        file << "    target " << std::to_string(out_id) << std::endl;
        file << "    sequence \"" << edge->sequence << "\"" << std::endl;
        file << "    Length " << edge->sequence.length() << std::endl;
        file << "    meta \"" << edge_meta_repr(edge->meta) << "\"" << std::endl;
        file << "  ]" << std::endl;

        // dummy nodes for tips
        /*
        if (in_null) {
            file << "  node [" << std::endl;
            file << "    id " << std::to_string(in_id) << std::endl;
            file << "    label \"null_" << std::to_string(in_id) << "\"" << std::endl;
            file << "  ]" << std::endl;
        }

        if (out_null) {
            file << "  node [" << std::endl;
            file << "    id " << std::to_string(out_id) << std::endl;
            file << "    label \"null_" << std::to_string(out_id) << "\"" << std::endl;
            file << "  ]" << std::endl;
        }
        */
    }

    file << "]";

    file.close();
    pdebug("closed file");
} 


void CompactEdgeFactory::write_fasta(const std::string filename) const {
    std::ofstream file;
    file.open(filename);
    pdebug("opened " << filename);
    for (auto edge_pair : compact_edges) {
        
        id_t edge_id = edge_pair.first;
        CompactEdge* edge = edge_pair.second;
        file << ">" << "edge_id=" << edge_id;
        file << " len=" << edge->sequence.length();
        file << " type=" << edge_meta_repr(edge->meta);
        file << " src=" << edge->in_node_id;
        file << " tgt=" << edge->out_node_id;
        file << std::endl;
        file << edge->sequence;
        file << std::endl;
    }

    file.close();
}

};
