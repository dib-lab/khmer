
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
_assemble_directed(AssemblerTraverser<direction> start_cursor,
                   StringVector& paths)
const
{
    std::string root_contig = linear_asm->assemble_directed<direction>(start_cursor);

    StringVector segments;
    std::vector< AssemblerTraverser<direction> > cursors;
    
    segments.push_back(root_contig);
    cursors.push_back(start_cursor);

    while(segments.size() != 0) {
    
        std::string segment = segments.back();
        AssemblerTraverser<direction> cursor = cursors.back();
#if DEBUG_ASSEMBLY
        std::cout << "Pop: " << segments.size() << " segments on stack." << std::endl;
        std::cout << "Segment: " << segment << std::endl;
        std::cout << "Cursor: " << cursor.cursor.repr(_ksize) << std::endl;
        std::cout << "n_filters: " << cursor.n_filters() << std::endl;
#endif
        segments.pop_back();
        cursors.pop_back();


        if (cursor.cursor_degree() > 1) {
            std::shared_ptr<LinkList> links = make_shared<LinkList>();
            linker->get_links(hdn, links);

            if (links->size() == 0) {
                paths.push_back(segment);
                continue;
            } else {
                
            }

        }

    }
}


}
