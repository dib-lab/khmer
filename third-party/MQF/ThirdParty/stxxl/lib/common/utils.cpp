/***************************************************************************
 *  lib/common/utils.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <stxxl/bits/common/types.h>
#include <stxxl/bits/common/utils.h>
#include <stxxl/bits/namespace.h>

#include <sstream>
#include <iomanip>

STXXL_BEGIN_NAMESPACE

//! Parse a string like "343KB" or "  44 GiB  " into the corresponding size in
//! bytes.
bool parse_SI_IEC_size(const std::string& str, uint64& size, char def_unit)
{
    char* endptr;
    size = strtoul(str.c_str(), &endptr, 10);
    if (!endptr) return false;          // parse failed, no number

    while (endptr[0] == ' ') ++endptr;  // skip over spaces

    // multiply with base ^ power
    unsigned int base = 1000;
    unsigned int power = 0;

    if (endptr[0] == 'k' || endptr[0] == 'K')
        power = 1, ++endptr;
    else if (endptr[0] == 'm' || endptr[0] == 'M')
        power = 2, ++endptr;
    else if (endptr[0] == 'g' || endptr[0] == 'G')
        power = 3, ++endptr;
    else if (endptr[0] == 't' || endptr[0] == 'T')
        power = 4, ++endptr;
    else if (endptr[0] == 'p' || endptr[0] == 'P')
        power = 5, ++endptr;

    // switch to power of two (only if power was set above)
    if ((endptr[0] == 'i' || endptr[0] == 'I') && power != 0)
        base = 1024, ++endptr;

    // byte indicator
    if (endptr[0] == 'b' || endptr[0] == 'B') {
        ++endptr;
    }
    else if (power == 0)
    {
        // no explicit power indicator, and no 'b' or 'B' -> apply default unit
        switch (def_unit)
        {
        default: break;
        case 'k': power = 1, base = 1000;
            break;
        case 'm': power = 2, base = 1000;
            break;
        case 'g': power = 3, base = 1000;
            break;
        case 't': power = 4, base = 1000;
            break;
        case 'p': power = 5, base = 1000;
            break;
        case 'K': power = 1, base = 1024;
            break;
        case 'M': power = 2, base = 1024;
            break;
        case 'G': power = 3, base = 1024;
            break;
        case 'T': power = 4, base = 1024;
            break;
        case 'P': power = 5, base = 1024;
            break;
        }
    }

    // skip over spaces
    while (endptr[0] == ' ') ++endptr;

    // multiply size
    for (unsigned int p = 0; p < power; ++p)
        size *= base;

    return (endptr[0] == 0);
}

std::string format_SI_size(uint64 number)
{
    // may not overflow, std::numeric_limits<uint64>::max() == 16 EiB
    double multiplier = 1000.0;
    static const char* SIendings[] = { "", "k", "M", "G", "T", "P", "E" };
    unsigned int scale = 0;
    double number_d = (double)number;
    while (number_d >= multiplier) {
        number_d /= multiplier;
        ++scale;
    }
    std::ostringstream out;
    out << std::fixed << std::setprecision(3) << number_d
        << ' ' << SIendings[scale];
    return out.str();
}

std::string format_IEC_size(uint64 number)
{
    // may not overflow, std::numeric_limits<uint64>::max() == 16 EiB
    double multiplier = 1024.0;
    static const char* IECendings[] = { "", "Ki", "Mi", "Gi", "Ti", "Pi", "Ei" };
    unsigned int scale = 0;
    double number_d = (double)number;
    while (number_d >= multiplier) {
        number_d /= multiplier;
        ++scale;
    }
    std::ostringstream out;
    out << std::fixed << std::setprecision(3) << number_d
        << ' ' << IECendings[scale];
    return out.str();
}

STXXL_END_NAMESPACE
