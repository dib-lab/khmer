#include "oxli/links.hh"

namespace oxli {

uint64_t Link::n_links = 0;


inline std::ostream& operator<< (std::ostream& stream,
                                 const Junction* j)
{
    stream << "Junction(u=" << j->u << ", v=" << j->v <<
        ", count=" << j->count << ", id=" << j->id() << ")";
}














}
