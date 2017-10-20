
#include <iostream>
#include "oxli/links.hh"

using namespace oxli;

namespace oxli {

uint64_t CompactNode::node_counter = 0;

/*
template <bool direction>
void LinkedAssembler::
_assemble_directed(AssemblerTraverser<direction> cursor,
                   StringVector& paths)
const
{
    std::string contig = linear_asm->assemble_directed<direction>(cursor);

    LinkTraversal link_traversal;

    while(1) {

#if DEBUG_ASSEMBLY
        std::cout << "Pop: " << segments.size() << " segments on stack." << std::endl;
        std::cout << "Segment: " << segment << std::endl;
        std::cout << "Cursor: " << cursor.cursor.repr(_ksize) << std::endl;
        std::cout << "n_filters: " << cursor.n_filters() << std::endl;
#endif

        if (cursor.cursor_degree() > 1) {
            // hit a HDN, push new links
            std::shared_ptr<LinkList> new_links;
            linker->get_links(hdn, new_links);
            // add them to traversal, with age being current contig size
            link_traversal.add_links(new_links, contig.size());

            if (links->size() == 0) {
                paths.push_back(segment);
                continue;
            }

            cursor.push_filter(get_link_filter(cursor.cursor,
                                               links,
                                               ages);
            KmerQueue linked_nodes;
            cursor.neighbors(linked_nodes);
            cursor.pop_filter();

            if (branch_starts.empty()) {
                paths.push_back(segment);
            }


        }
        break;

    }
}
*/
}
