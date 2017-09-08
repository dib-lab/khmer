#include <string>

namespace oxli {

std::string get_version_cpp()
{
#define _macro_xstr(s) _macro_str(s)
#define _macro_str(s) #s
    std::string dVersion = _macro_xstr(VERSION);
    return dVersion;
}

}
