
#include <iostream>
#include "oxli/links.hh"

using namespace oxli;

namespace oxli {

void CompactEdgeFactory::write_gml(const std::string filename,
                                   const CompactNodeFactory& nodes) const {

    std::ofstream file;
    file.open(filename);
    pdebug("opened " << filename);
    file << "graph 1" << std::endl << "[" << std::endl;

    pdebug("writing " << nodes.n_nodes() << " nodes");
    for (auto node : nodes.compact_nodes) {
        file << "  node [" << std::endl;
        file << "    id " << std::to_string(node.node_id) << std::endl;
        file << "    kmer \"" << node.sequence << "\"" << std::endl;
        file << "    count \"" << std::to_string(node.count) << "\"" << std::endl;
        file << "  ]" << std::endl;
    }

    pdebug("writing " << compact_edges.size() << " edges");
    for (auto edge_pair : compact_edges) {
        
        id_t edge_id = edge_pair.first;
        CompactEdge* edge = edge_pair.second;
        
        file << "  edge [" << std::endl;
        file << "    id " << std::to_string(edge->edge_id) << std::endl;

        id_t in_id, out_id;
        bool in_null = false, out_null = false;
        if (edge->in_node_id == NULL_ID) {
            in_id = NULL_ID - edge->edge_id;
            in_null = true;
        } else {
            in_id = edge->in_node_id;
        }
        if(edge->out_node_id == NULL_ID) {
            out_id = NULL_ID - edge->edge_id;
            out_null = true;
        } else {
            out_id = edge->out_node_id;
        }

        file << "    source " << std::to_string(in_id) << std::endl;
        file << "    target " << std::to_string(out_id) << std::endl;
        file << "    sequence \"" << edge->sequence << "\"" << std::endl;
        file << "  ]" << std::endl;

        // dummy nodes for tips
        if (in_null) {
            file << "  node [" << std::endl;
            file << "    id " << std::to_string(in_id) << std::endl;
            file << "    label \"null\"" << std::endl;
            file << "  ]" << std::endl;
        }

        if (out_null) {
            file << "  node [" << std::endl;
            file << "    id " << std::to_string(out_id) << std::endl;
            file << "    label \"null\"" << std::endl;
            file << "  ]" << std::endl;
        }
    }

    file << "]";

    file.close();
    pdebug("closed file");
} 

};
