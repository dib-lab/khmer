
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




}
