//
// This file is part of khmer, http://github.com/ged-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2013. It is licensed under
// the three-clause BSD license; see doc/LICENSE.txt.
// Contact: khmer-project@idyll.org
//

#ifndef KHMER_EXCEPTION_HH
#define KHMER_EXCEPTION_HH

#include <exception>

namespace khmer
{

/// 
// A base class for all exceptions.
//
// All exceptions should be derived from this base class.
//
class khmer_exception : public std::exception
{
public:
    khmer_exception(const char * msg) : _msg(msg) { };

    virtual const char* what() const throw() {
        return _msg;
    }
protected:
    const char * _msg;
};

class khmer_file_exception : public khmer_exception
{
public:
    khmer_file_exception(const char * msg) : khmer_exception(msg) { };
};

}

#endif // KHMER_EXCEPTION_HH

// vim: set sts=2 sw=2:
