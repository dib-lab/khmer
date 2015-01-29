//
// This file is part of khmer, http://github.com/ged-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2013. It is licensed under
// the three-clause BSD license; see doc/LICENSE.txt.
// Contact: khmer-project@idyll.org
//

#ifndef KHMER_EXCEPTION_HH
#define KHMER_EXCEPTION_HH

#include <exception>
#include <string>

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
    explicit khmer_exception(const char * msg) : _msg(msg) { }
    explicit khmer_exception(const std::string& msg = "Generic khmer exception")
        : _msg(msg.c_str()) { }

    virtual ~khmer_exception() throw() { }
    virtual const char* what() const throw ()
    {
        return _msg;
    }

protected:
    const char * _msg;
};

///
// A base class for file exceptions.
//
class khmer_file_exception : public khmer_exception
{
public:
    explicit khmer_file_exception(const char * msg) : khmer_exception(msg) { }
    explicit khmer_file_exception(const std::string& msg)
        : khmer_exception(msg) { }
};

struct InvalidStreamBuffer : public khmer_exception {
};

class InvalidStreamHandle : public khmer_file_exception
{
public:
    InvalidStreamHandle()
        : khmer_file_exception("Generic InvalidStreamHandle error") {}
    InvalidStreamHandle(const char * msg) : khmer_file_exception(msg) {}
};

class StreamReadError : public khmer_file_exception
{
public:
    StreamReadError()
        : khmer_file_exception("Generic StreamReadError error") {}
    StreamReadError(const char * msg) : khmer_file_exception(msg) {}
};

}

#endif // KHMER_EXCEPTION_HH

// vim: set sts=2 sw=2:
