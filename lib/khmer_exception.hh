/*
This file is part of khmer, https://github.com/dib-lab/khmer/, and is
Copyright (C) 2014-2015, Michigan State University.
Copyright (C) 2015, The Regents of the University of California.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.

    * Neither the name of the Michigan State University nor the names
      of its contributors may be used to endorse or promote products
      derived from this software without specific prior written
      permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
LICENSE (END)

Contact: khmer-project@idyll.org
*/
#ifndef KHMER_EXCEPTION_HH
#define KHMER_EXCEPTION_HH

#include <exception>
#include <string>

namespace khmer
{

///
// A base class for all exceptions.
//
// All exceptions should be derived from this base class or a sub-class
//
class khmer_exception : public std::exception
{
public:
    explicit khmer_exception(const std::string& msg = "Generic khmer exception")
        : _msg(msg) { }

    virtual ~khmer_exception() throw() { }
    virtual const char* what() const throw ()
    {
        return _msg.c_str();
    }

protected:
    const std::string _msg;
};


/////// Base Exceptions /////

///
// A base class for file exceptions.
//
class khmer_file_exception : public khmer_exception
{
public:
    explicit khmer_file_exception(const std::string& msg)
        : khmer_exception(msg) { }
};

// A base exception for value exceptions
class khmer_value_exception : public khmer_exception
{
public:
    explicit khmer_value_exception(const std::string& msg)
        : khmer_exception(msg) { }
};

/////// Specialised Exceptions /////

class InvalidStream : public khmer_file_exception
{
public:
    InvalidStream()
        : khmer_file_exception("Generic InvalidStream error") {}
    explicit InvalidStream(const std::string& msg)
        : khmer_file_exception(msg) {}
};

class StreamReadError : public khmer_file_exception
{
public:
    StreamReadError()
        : khmer_file_exception("Generic StreamReadError error") {}
    explicit StreamReadError(const std::string& msg)
        : khmer_file_exception(msg) {}
};


///
// An exception for invalid arguments to functions
//

class InvalidValue : public khmer_value_exception
{
public:
    explicit InvalidValue(const std::string& msg)
        : khmer_value_exception(msg) { }
};

///
// An exception for trying to change a read-only attributes
//

class ReadOnlyAttribute : public khmer_exception
{
public:
    explicit ReadOnlyAttribute(const std::string& msg)
        : khmer_exception(msg) { }
};

} // end namespace khmer

#endif // KHMER_EXCEPTION_HH

// vim: set sts=2 sw=2:
