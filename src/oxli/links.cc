
#include <iostream>
#include "oxli/links.hh"

using namespace oxli;

std::ostream& operator<<(std::ostream& stream,
                         const oxli::Junction& j)
{
    stream << "Junction(u=" << j.u << ", v=" << j.v <<
        ", count=" << j.count << ", id=" << j.id() << ")";
    return stream;
}

namespace oxli {

uint64_t Link::n_links = 0;

template <bool direction>
void LinkedAssembler::
_assemble_directed(AssemblerTraverser<direction> cursor,
                   StringVector& paths)
const
{
    std::string root_contig = linear_asm->assemble_directed<direction>(cursor);

    StringVector segments;
    
    std::shared_ptr<LinkList> links = make_shared<LinkList>();
    std::shared_ptr< std::list<uint64_t> > ages = make_shared<std::list< uint64_t> >();
    
    paths.push_back(root_contig);

    while(1) {

#if DEBUG_ASSEMBLY
        std::cout << "Pop: " << segments.size() << " segments on stack." << std::endl;
        std::cout << "Segment: " << segment << std::endl;
        std::cout << "Cursor: " << cursor.cursor.repr(_ksize) << std::endl;
        std::cout << "n_filters: " << cursor.n_filters() << std::endl;
#endif

/*
        if (cursor.cursor_degree() > 1) {
            linker->get_links_copy(hdn, links);

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
        */
        break;

    }
}


}
